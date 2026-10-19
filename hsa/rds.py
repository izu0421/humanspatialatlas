"""Convert Seurat .rds deposits into the atlas's .h5ad form.

376 GB of the corpus -- 914 files across 123 accessions -- arrives as Seurat objects, which the
loader has always rejected as `no_recognised_matrix`. They are not unreadable, only inconvenient:
an R process extracts the pieces to plain files and this module assembles them here.

The R side (rds_extract.R) deliberately avoids the Seurat and SeuratObject packages. Neither is
installed, and depending on them would tie every conversion to the package version the depositor
happened to use; reading the object's attributes needs base R plus Matrix and works across the v3,
v4 and v5 layouts alike.

Worth knowing: these objects often carry richer cell metadata than the deposit's own files --
donor, age, sex and diagnosis columns that the GEO characteristics never mention -- so the obs
table is kept in full rather than reduced to coordinates.
"""
from __future__ import annotations

import json
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

HERE = Path(__file__).resolve().parent
RSCRIPT = HERE / "rds_extract.R"
TIMEOUT = 3600


def extract(rds: str | Path, out_dir: str | Path) -> dict:
    """Run the R side; returns its info dict."""
    out_dir = Path(out_dir)
    r = subprocess.run(["Rscript", str(RSCRIPT), str(rds), str(out_dir)],
                       capture_output=True, text=True, timeout=TIMEOUT)
    info_f = out_dir / "info.json"
    if r.returncode != 0 or not info_f.exists():
        raise RuntimeError(f"R extraction failed ({r.returncode}): "
                           f"{(r.stderr or r.stdout)[-300:]}")
    return json.loads(info_f.read_text())


def assemble(stage: str | Path, info: dict):
    """Build an AnnData from the staged files, oriented cells x genes."""
    import anndata as ad
    import scipy.io as sio
    import scipy.sparse as sp

    stage = Path(stage)
    X = sio.mmread(stage / "counts.mtx").tocsr()          # genes x cells, as Seurat stores it
    feats = pd.read_csv(stage / "features.tsv", header=None)[0].astype(str).values
    cells = pd.read_csv(stage / "barcodes.tsv", header=None)[0].astype(str).values
    if X.shape == (len(cells), len(feats)) and len(cells) != len(feats):
        pass                                              # already cells x genes
    else:
        X = X.T.tocsr()
    a = ad.AnnData(sp.csr_matrix(X.astype(np.float32)),
                   obs=pd.DataFrame(index=pd.Index(cells, dtype=str)),
                   var=pd.DataFrame(index=pd.Index(feats, dtype=str)))
    a.obs_names_make_unique(); a.var_names_make_unique()

    obs_f = stage / "obs.csv"
    if obs_f.exists():
        md = pd.read_csv(obs_f, index_col=0, low_memory=False)
        md.index = md.index.astype(str)
        # Seurat row order matches the matrix columns; align by name where possible, by
        # position only when the names do not line up at all.
        if md.index.isin(a.obs_names).mean() > 0.5:
            md = md.reindex(a.obs_names)
        elif len(md) == a.n_obs:
            md.index = a.obs_names
        else:
            md = None
        if md is not None:
            md = md.loc[:, [c for c in md.columns if md[c].notna().any()]]
            for c in md.columns:
                if md[c].dtype == object:
                    md[c] = md[c].astype(str)
            a.obs = md

    co_f = stage / "coords.csv"
    if co_f.exists():
        co = pd.read_csv(co_f, index_col=0)
        co.index = co.index.astype(str)
        # column case varies with where the coordinates came from (meta.data keeps the
        # depositor's spelling, e.g. X/Y or sdimx/sdimy)
        co.columns = [str(c).lower() for c in co.columns]
        if co.index.isin(a.obs_names).mean() > 0.5:
            co = co.reindex(a.obs_names)
        elif len(co) != a.n_obs:
            co = None
        if co is not None and {"x", "y"} <= set(co.columns):
            xy = co[["x", "y"]].to_numpy(dtype=float)
            if np.isfinite(xy).all(1).sum() >= 0.5 * len(xy):
                a.obsm["spatial"] = xy
            if "slide" in co.columns and co.slide.nunique() > 1:
                # a merged object: record which slide each cell came from, because its coordinate
                # frames are per-slide and a statistic computed across them is meaningless
                a.obs["seurat_slide"] = co.slide.astype(str).values
                a.uns["merged_slides"] = sorted(co.slide.dropna().astype(str).unique().tolist())
    a.uns["rds_info"] = {k: str(v) for k, v in info.items()}
    return a


def convert(rds: str | Path, keep_stage: bool = False):
    """Full path: .rds -> AnnData. Raises on anything that would give a silent half-sample."""
    stage = Path(tempfile.mkdtemp(prefix="rds_"))
    try:
        info = extract(rds, stage)
        a = assemble(stage, info)
        if a.n_obs == 0 or a.n_vars == 0:
            raise ValueError(f"empty after assembly: {a.shape}")
        if "spatial" not in a.obsm:
            raise ValueError(f"no coordinates recovered (R reported {info.get('coord_source')})")
        return a, info
    finally:
        if not keep_stage:
            shutil.rmtree(stage, ignore_errors=True)
