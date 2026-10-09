"""Make the corpus consistent: unpack stray archives, then load every sample into one contract.

The contract, one .h5ad per sample:
    X                  raw integer counts, cells (or spots/bins) x genes
    var_names          gene symbols; var['ensembl_id'] where resolvable
    obsm['spatial']    float32 (n, 2) x/y
    uns['hsa']         source, accession, sample_id, technology, coord_units, provenance
    obs                harmonised sample-level metadata broadcast per cell

Nothing here guesses: a sample that cannot meet the contract is recorded with the reason and left out
rather than patched into something that looks valid.
"""
from __future__ import annotations

import gzip
import json
import logging
import re
import shutil
import tarfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import db
from .config import DATA

logger = logging.getLogger(__name__)
STD = DATA.parent / "standard"

# archives that were fetched as matrix/coords files and still need opening locally
ARCHIVE_RE = re.compile(r"\.(tar\.gz|tgz|tar|zip)$", re.I)
KEEP_MEMBER = re.compile(
    r"tissue_positions|scalefactors_json|barcodes\.tsv|features\.tsv|genes\.tsv|matrix\.mtx"
    r"|\.h5ad$|\.h5$|cells\.(csv|parquet)|cell_metadata|exprmat|cell_by_gene", re.I)


# ---------------------------------------------------------------- local unpacking
def unpack_local_archives(limit: int | None = None) -> pd.DataFrame:
    """Open already-downloaded .tar.gz/.zip files in place and register their useful members.

    These arrived as role='matrix'/'coords' (e.g. Visium `*_spatial.tar.gz`) so the bundle streamer
    never saw them, leaving samples with a matrix but no reachable coordinates.
    """
    rows = db.query("SELECT * FROM files WHERE dl_status='done' AND local_path IS NOT NULL "
                    "AND role IN ('matrix','coords','matrix+coords')")
    todo = [r for r in rows if ARCHIVE_RE.search(r["local_path"] or "")]
    if limit:
        todo = todo[:limit]
    logger.info("%d downloaded archives to unpack in place", len(todo))
    out = []
    for i, r in enumerate(todo):
        src = Path(r["local_path"])
        if not src.exists():
            continue
        dest = src.parent / (ARCHIVE_RE.sub("", src.name) + "_unpacked")
        n = 0
        try:
            if src.suffix.lower() == ".zip":
                with zipfile.ZipFile(src) as z:
                    members = [m for m in z.namelist() if KEEP_MEMBER.search(m)]
                    for m in members:
                        t = dest / Path(m).name
                        t.parent.mkdir(parents=True, exist_ok=True)
                        if not t.exists():
                            with z.open(m) as fsrc, open(t, "wb") as fdst:
                                shutil.copyfileobj(fsrc, fdst)
                        n += 1
            else:
                with tarfile.open(src, "r:*") as tf:
                    for m in tf:
                        if not m.isfile() or not KEEP_MEMBER.search(m.name):
                            continue
                        t = dest / Path(m.name).name
                        t.parent.mkdir(parents=True, exist_ok=True)
                        if not t.exists():
                            with tf.extractfile(m) as fsrc, open(t, "wb") as fdst:
                                shutil.copyfileobj(fsrc, fdst)
                        n += 1
            status = "ok" if n else "no_useful_members"
        except Exception as e:                       # corrupt or truncated archive
            status, n = f"error: {type(e).__name__}", 0
        for t in (dest.glob("*") if dest.exists() else []):
            db.execute("INSERT OR IGNORE INTO files (source, accession, sample_id, url, role, fmt, "
                       "size_bytes, local_path, dl_status, dl_msg) VALUES (?,?,?,?,?,?,?,?,?,?)",
                       (r["source"], r["accession"], r["sample_id"], f"{r['url']}#local={t.name}",
                        _role_of(t.name), t.suffix.lstrip("."), t.stat().st_size, str(t), "done",
                        f"unpacked locally from {src.name}"))
        out.append({"path": str(src), "members": n, "status": status})
        if i % 100 == 0:
            logger.info("[%d/%d] %s %s", i + 1, len(todo), status, src.name)
    return pd.DataFrame(out)


def _role_of(name: str) -> str:
    from .tools_api import classify
    return classify(name)


# ---------------------------------------------------------------- per-sample loading
def sample_files(source: str, accession: str, sample_id: str) -> list[Path]:
    rows = db.query("SELECT local_path FROM files WHERE source=? AND accession=? AND sample_id=? "
                    "AND dl_status='done' AND local_path IS NOT NULL", (source, accession, sample_id))
    return [Path(r["local_path"]) for r in rows if Path(r["local_path"]).exists()]


def _find(paths: list[Path], pattern: str) -> Path | None:
    rx = re.compile(pattern, re.I)
    hits = [p for p in paths if rx.search(p.name)]
    # prefer filtered over raw, and shallower paths
    hits.sort(key=lambda p: ("raw_feature" in p.name.lower(), len(p.parts)))
    return hits[0] if hits else None


def load_sample(source: str, accession: str, sample_id: str):
    """-> (AnnData | None, info dict). Never raises; failures are reported in info['status']."""
    import anndata as ad
    import scanpy as sc

    paths = sample_files(source, accession, sample_id)
    info = {"source": source, "accession": accession, "sample_id": sample_id,
            "n_files": len(paths), "loader": None, "status": "no_files"}
    if not paths:
        return None, info

    a = None
    try:
        if (p := _find(paths, r"\.h5ad$")):
            a, info["loader"] = ad.read_h5ad(p), "h5ad"
        elif (p := _find(paths, r"(cell_feature_matrix|filtered_feature_bc_matrix|feature_bc_matrix)\.h5$")):
            a, info["loader"] = sc.read_10x_h5(p), "10x_h5"
        elif (p := _find(paths, r"matrix\.mtx(\.gz)?$")):
            a, info["loader"] = _read_mtx_trio(p, paths), "mtx_trio"
        elif (p := _find(paths, r"(exprmat|cell_by_gene|stdata|_counts?\.|expression)")):
            a, info["loader"] = _read_csv_matrix(p), "csv_matrix"
        else:
            info["status"] = "no_recognised_matrix"
            return None, info
    except Exception as e:
        info["status"] = f"matrix_error: {type(e).__name__}: {e}"[:200]
        return None, info

    a.var_names_make_unique()
    a.obs_names_make_unique()
    coords, cinfo = _load_coords(paths, a)
    info.update(cinfo)
    if coords is not None:
        a.obsm["spatial"] = coords
    info["status"] = "ok" if coords is not None else "no_coords"
    info["n_cells"], info["n_genes"] = a.shape
    return a, info


def _read_mtx_trio(mtx: Path, paths: list[Path]):
    import scanpy as sc
    import anndata as ad
    a = ad.AnnData(sp.csr_matrix(sc.read_mtx(mtx).X).T.tocsr())
    sib = [p for p in paths if p.parent == mtx.parent]
    bc = _find(sib, r"barcodes\.tsv")
    ft = _find(sib, r"(features|genes)\.tsv")
    if bc is not None:
        a.obs_names = pd.read_csv(bc, header=None, sep="\t")[0].astype(str).values[: a.n_obs]
    if ft is not None:
        f = pd.read_csv(ft, header=None, sep="\t")
        a.var_names = f[min(1, f.shape[1] - 1)].astype(str).values[: a.n_vars]
        if f.shape[1] > 1:
            a.var["ensembl_id"] = f[0].astype(str).values[: a.n_vars]
    return a


def _read_csv_matrix(p: Path):
    import anndata as ad
    df = pd.read_csv(p, index_col=0, compression="gzip" if p.name.endswith(".gz") else None)
    # orient so rows are cells: gene tables are usually taller than wide
    if df.shape[0] > df.shape[1] * 3 and df.shape[1] < 2000:
        pass
    elif df.shape[1] > df.shape[0] * 3:
        df = df.T
    num = df.select_dtypes(include=[np.number])
    return ad.AnnData(sp.csr_matrix(num.to_numpy(dtype=np.float32)),
                      obs=pd.DataFrame(index=num.index.astype(str)),
                      var=pd.DataFrame(index=num.columns.astype(str)))


def _read_coord_table(p: Path) -> pd.DataFrame:
    """Read a coordinate table, handling the headerless Visium tissue_positions_list.csv.

    That file has no header, so a default read silently eats the first spot and turns barcodes into
    column names. Detect it by testing whether the first field of row 0 looks like a barcode.
    """
    # GEO re-gzips 10x parquet deposits, so the Xenium centroid table usually arrives as
    # cells.parquet.gz. p.suffix is then ".gz" and a plain suffix test sends a binary parquet
    # into read_csv, whose exception was being swallowed as "no coordinates" -- 716 Xenium
    # samples had their centroids on disk the whole time.
    if ".parquet" in p.name:
        if p.name.endswith(".gz"):
            import gzip as _gz
            import io as _io
            with _gz.open(p, "rb") as fh:
                return pd.read_parquet(_io.BytesIO(fh.read()))
        return pd.read_parquet(p)
    comp = "gzip" if p.name.endswith(".gz") else None
    head = pd.read_csv(p, compression=comp, nrows=1, header=None)
    first = str(head.iloc[0, 0])
    headerless = bool(re.fullmatch(r"[ACGT]{8,}-?\d*", first)) or first.replace(".", "").isdigit()
    if headerless:
        df = pd.read_csv(p, compression=comp, header=None)
        if df.shape[1] >= 6:                       # the Space Ranger column order
            df.columns = (["barcode", "in_tissue", "array_row", "array_col",
                           "pxl_row_in_fullres", "pxl_col_in_fullres"] + list(df.columns[6:]))
        return df
    return pd.read_csv(p, compression=comp)


def _load_coords(paths: list[Path], a) -> tuple[np.ndarray | None, dict]:
    """Coordinates from whatever convention the platform used; returns (coords, info)."""
    if "spatial" in getattr(a, "obsm", {}):
        return np.asarray(a.obsm["spatial"], dtype=np.float32)[:, :2], {"coord_source": "h5ad_obsm"}
    for pat, kind in [(r"tissue_positions", "visium"), (r"cells\.(parquet|csv)", "xenium"),
                      (r"(cell_)?metadata_file|cell_metadata", "cosmx"),
                      (r"_coords|bead|position|centroid", "generic")]:
        p = _find(paths, pat)
        if p is None:
            continue
        try:
            df = _read_coord_table(p)
            # index by barcode/cell id FIRST, so the xy slice carries the right index
            if df.index.dtype.kind in "iu" and df.shape[1] > 2:
                df = df.set_index(df.columns[0])
            df.index = df.index.astype(str)
            xy = _pick_xy(df)
            if xy is None:
                continue
            common = a.obs_names.intersection(df.index)
            if len(common) >= 0.5 * a.n_obs:
                a._inplace_subset_obs(a.obs_names.isin(common))
                return xy.loc[a.obs_names].to_numpy(np.float32), {"coord_source": kind}
            if len(df) == a.n_obs:                  # ids do not match but the rows line up
                return xy.to_numpy(np.float32), {"coord_source": f"{kind}_by_order"}
        except Exception as e:
            logger.debug("coord read failed for %s: %s", p.name, e)
            continue
    return None, {"coord_source": None}


def _pick_xy(df: pd.DataFrame):
    cands = [("x_centroid", "y_centroid"), ("pxl_col_in_fullres", "pxl_row_in_fullres"),
             ("CenterX_global_px", "CenterY_global_px"), ("x_global_px", "y_global_px"),
             ("x", "y"), ("X", "Y"), ("xcoord", "ycoord"), ("array_col", "array_row"),
             ("spatial_1", "spatial_2"), ("imagecol", "imagerow")]
    low = {c.lower(): c for c in df.columns}
    for cx, cy in cands:
        if cx.lower() in low and cy.lower() in low:
            return df[[low[cx.lower()], low[cy.lower()]]]
    # Visium tissue_positions without a header
    if df.shape[1] >= 6 and all(pd.api.types.is_numeric_dtype(df[c]) for c in df.columns[-2:]):
        return df[df.columns[-2:]]
    return None


# ---------------------------------------------------------------- QC + standardisation at scale
QC_CSV = DATA.parent / "runs" / "sample_qc.csv"


def qc_one(args) -> dict:
    """Load one sample and report what it is. Never raises."""
    source, accession, sample_id, tech = args
    import numpy as np
    rec = {"source": source, "accession": accession, "sample_id": sample_id, "technology": tech}
    try:
        a, info = load_sample(source, accession, sample_id)
    except Exception as e:
        rec.update(status=f"loader_crash: {type(e).__name__}: {e}"[:160])
        return rec
    rec.update({k: v for k, v in info.items() if k not in ("source", "accession", "sample_id")})
    if a is None:
        return rec
    try:
        X = a.layers["counts"] if "counts" in getattr(a, "layers", {}) else a.X
        tot = np.asarray(X.sum(1)).ravel()
        dat = X.data if sp.issparse(X) else np.asarray(X).ravel()
        sub = dat[:200000]
        rec["median_counts"] = float(np.median(tot)) if len(tot) else 0.0
        rec["median_genes"] = float(np.median(np.asarray((X > 0).sum(1)).ravel())) if a.n_obs else 0.0
        rec["is_integer"] = bool(len(sub) == 0 or np.allclose(sub, np.round(sub)))
        rec["frac_zero_cells"] = float((tot == 0).mean()) if len(tot) else 1.0
        if "spatial" in a.obsm:
            xy = np.asarray(a.obsm["spatial"], dtype=float)
            rec["x_range"] = float(np.ptp(xy[:, 0])) if len(xy) else 0.0
            rec["y_range"] = float(np.ptp(xy[:, 1])) if len(xy) else 0.0
            rec["n_unique_coords"] = int(len(np.unique(xy, axis=0)))
        rec["gene_id_style"] = ("ensembl" if str(a.var_names[0]).startswith("ENS")
                                else "symbol" if a.n_vars else "unknown")
    except Exception as e:
        rec["status"] = f"{rec.get('status','')} | qc_error: {type(e).__name__}"[:160]
    return rec


def qc_all(workers: int = 12, limit: int | None = None, resume: bool = True,
           retry: str | None = None) -> pd.DataFrame:
    """Run the loader over every ready sample; write a per-sample QC table (resumable).

    `retry` re-runs only samples whose recorded status starts with that prefix, which is how a
    loader fix is rolled out: after teaching the reader a new format, `retry="no_coords"` revisits
    exactly the samples that format broke instead of re-reading the whole corpus. The table is
    de-duplicated on (source, accession, sample_id) afterwards, last write winning.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from .export import norm_tech

    ds = {(r["source"], r["accession"]): norm_tech(r["technology"])
          for r in db.query("SELECT source, accession, technology FROM datasets "
                            "WHERE verdict IN ('ready','bundle_only')")}
    rows = db.query("SELECT DISTINCT source, accession, sample_id FROM files "
                    "WHERE dl_status='done' AND local_path IS NOT NULL")
    todo = [(r["source"], r["accession"], r["sample_id"], ds.get((r["source"], r["accession"]), "Unknown"))
            for r in rows if (r["source"], r["accession"]) in ds]

    if retry and QC_CSV.exists():
        prev = pd.read_csv(QC_CSV)
        bad = prev[prev.status.astype(str).str.startswith(retry)]
        keys = set(map(tuple, bad[["source", "accession", "sample_id"]].astype(str).values))
        todo = [t for t in todo if (str(t[0]), str(t[1]), str(t[2])) in keys]
        logger.info("retrying %d samples previously recorded as %s*", len(todo), retry)
    elif resume and QC_CSV.exists():
        prev = pd.read_csv(QC_CSV)
        done_keys = set(map(tuple, prev[["source", "accession", "sample_id"]].astype(str).values))
        todo = [t for t in todo if (str(t[0]), str(t[1]), str(t[2])) not in done_keys]
        logger.info("resuming: %d already done, %d to go", len(done_keys), len(todo))
    if limit:
        todo = todo[:limit]

    out, header = [], not QC_CSV.exists()
    QC_CSV.parent.mkdir(exist_ok=True, parents=True)
    with ProcessPoolExecutor(workers) as ex:
        futs = [ex.submit(qc_one, t) for t in todo]
        for i, fu in enumerate(as_completed(futs)):
            out.append(fu.result())
            if len(out) >= 200 or i == len(futs) - 1:
                pd.DataFrame(out).to_csv(QC_CSV, mode="a", header=header, index=False)
                header, out = False, []
                logger.info("[%d/%d] written", i + 1, len(futs))
    df = pd.read_csv(QC_CSV)
    if retry:                                     # last write wins for a re-run sample
        n0 = len(df)
        df = df.drop_duplicates(subset=["source", "accession", "sample_id"], keep="last")
        df.to_csv(QC_CSV, index=False)
        logger.info("de-duplicated QC table: %d -> %d rows", n0, len(df))
    return df
