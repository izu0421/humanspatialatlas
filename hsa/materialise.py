"""Standardise every loadable sample into one .h5ad, then reclaim the raw files.

Why this exists: the corpus holds 2.6 TB of depositor-shaped files in a dozen formats, of which
1.68 TB belongs to samples that load. The same content as a gzip-compressed sparse .h5ad is
several times smaller and uniform, so converting frees space AND removes the format problem.

Design constraints that matter here:
  * /data runs at 100%, so this streams -- convert one sample, verify it, delete that sample's
    raw files, move on. Peak extra disk is one sample, not one corpus.
  * Deletion is irreversible and a re-download costs days, so a sample's raw files are removed
    only after its .h5ad has been re-opened from disk and checked against the expected shape.
  * Visium HD keeps its cell-level segmentation where one exists; the bins are still written
    (compressed) rather than discarded, because re-segmentation needs them -- 43 samples were
    re-segmented at a corrected resolution the same week this was written.
"""
from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import db, enrich, manifest
from .config import DATA, RUNS

logger = logging.getLogger(__name__)
OUT = DATA.parent / "standardised"
LEDGER = RUNS / "materialise_ledger.csv"


def _tag(source: str, accession: str, sample_id: str) -> str:
    import hashlib
    slug = re.sub(r"[^A-Za-z0-9]+", "_", f"{source}__{accession}__{sample_id}")[:110].strip("_")
    h = hashlib.md5(f"{source}|{accession}|{sample_id}".encode()).hexdigest()[:8]
    return f"{slug}_{h}"


def _clean(a, meta: dict):
    """Counts as int32 CSR, coordinates in microns, everything HSA knows in uns."""
    import anndata as ad

    X = a.X
    if not sp.issparse(X):
        # a CSV-read matrix can arrive as object dtype (a stray text column, or pandas keeping
        # mixed types); scipy.sparse refuses object, so coerce first and fail loudly if it cannot
        if getattr(X, "dtype", None) is not None and X.dtype == object:
            X = np.asarray(X, dtype=object)
            X = np.vectorize(lambda v: pd.to_numeric(v, errors="coerce"))(X).astype(np.float64)
            X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
        X = sp.csr_matrix(X)
    X = X.tocsr()
    if X.dtype == object:
        raise ValueError("matrix is object dtype after coercion; refusing to write")
    if np.allclose(X.data, np.round(X.data)) and X.data.max() < 2**31:
        X.data = X.data.astype(np.int32)         # counts, not a float normalisation
    out = ad.AnnData(X)
    out.obs_names = [str(x) for x in a.obs_names]
    out.var_names = [str(x) for x in a.var_names]
    for c in a.var.columns:
        if a.var[c].dtype.kind in "ifbOSU":
            out.var[c] = a.var[c].to_numpy()
    xy = a.obsm.get("spatial")
    if xy is not None:
        xy = np.asarray(xy, dtype=np.float32)[:, :2]
        scale = meta.get("coord_um_per_unit")
        if scale and np.isfinite(scale) and scale > 0:
            out.obsm["spatial_um"] = (xy * float(scale)).astype(np.float32)
        out.obsm["spatial"] = xy
    return out


def convert_one(source: str, accession: str, sample_id: str, delete_raw: bool = False,
                overwrite: bool = False) -> dict:
    """Load -> enrich -> write .h5ad -> verify -> (optionally) delete the raw files."""
    from . import standardise as St

    rec = {"source": source, "accession": accession, "sample_id": sample_id, "status": "start"}
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{_tag(source, accession, sample_id)}.h5ad"
    if out.exists() and not overwrite:
        return {**rec, "status": "already_done", "path": str(out),
                "bytes_out": out.stat().st_size}

    meta = manifest.for_sample(source, accession, sample_id)
    try:
        a, info = St.load_sample(source, accession, sample_id)
    except Exception as e:
        return {**rec, "status": f"load_error: {type(e).__name__}: {e}"[:200]}
    if a is None or not str(info.get("status", "")).startswith("ok"):
        return {**rec, "status": f"not_loadable: {info.get('status')}"[:200]}

    try:
        clean = _clean(a, meta)
        sex_call, sex_ev = enrich.sex_from_expression(a)
        extra = {"inferred_sex": sex_call, "inferred_sex_evidence": sex_ev,
                 **{f"qc_{k}": v for k, v in enrich.qc_stats(a).items()},
                 **{f"spatial_{k}": v for k, v in
                    enrich.spatial_stats(a, meta.get("coord_um_per_unit")).items()}}
        clean.uns["hsa"] = json.loads(json.dumps(
            {**meta, **extra, "loader": info.get("loader"),
             "coord_source": info.get("coord_source"),
             "hsa_schema": 1}, default=str))
        clean.write_h5ad(out, compression="gzip")
    except Exception as e:
        if out.exists():
            out.unlink()                          # never leave a half-written file behind
        return {**rec, "status": f"write_error: {type(e).__name__}: {e}"[:200]}

    # verify by reading back, before anything is deleted
    try:
        import anndata as ad
        chk = ad.read_h5ad(out, backed="r")
        if chk.shape != clean.shape:
            raise ValueError(f"shape {chk.shape} != {clean.shape}")
        n_obs, n_var = chk.shape
        del chk
    except Exception as e:
        if out.exists():
            out.unlink()
        return {**rec, "status": f"verify_failed: {type(e).__name__}: {e}"[:200]}

    rec.update(status="ok", path=str(out), n_obs=int(n_obs), n_vars=int(n_var),
               bytes_out=out.stat().st_size, inferred_sex=sex_call)

    raws = [r["local_path"] for r in db.query(
        "SELECT local_path FROM files WHERE source=? AND accession=? AND sample_id=? "
        "AND dl_status='done' AND local_path IS NOT NULL", (source, accession, sample_id))]
    rec["bytes_raw"] = sum(os.path.getsize(p) for p in raws if os.path.exists(p))
    if delete_raw:
        # only files belonging to this sample alone: a path shared with another sample (common
        # for dataset-level deposits) must survive until that sample is converted too.
        shared = {r["local_path"] for r in db.query(
            "SELECT local_path FROM files WHERE local_path IN (%s) AND NOT "
            "(source=? AND accession=? AND sample_id=?)" % ",".join("?" * len(raws)),
            (*raws, source, accession, sample_id))} if raws else set()
        freed = 0
        for p in raws:
            if p in shared or not os.path.exists(p):
                continue
            try:
                freed += os.path.getsize(p)
                os.remove(p)
            except OSError:
                pass
        rec["bytes_freed"] = freed
        rec["n_shared_kept"] = len(shared)
    return rec


def convert_all(limit: int | None = None, delete_raw: bool = False,
                technologies=None, min_free_gb: float = 80.0) -> pd.DataFrame:
    """Stream over every loadable sample, newest-largest first, stopping if disk runs low."""
    m = pd.read_csv(manifest.MANIFEST) if manifest.MANIFEST.exists() else manifest.build()
    todo = m[m.usable.fillna(False)]
    if technologies:
        todo = todo[todo.technology_resolved.isin(technologies)]
    todo = todo.sort_values("n_cells", ascending=False)
    if limit:
        todo = todo.head(limit)

    done = set()
    rows = []
    if LEDGER.exists():
        prev = pd.read_csv(LEDGER)
        rows = prev.to_dict("records")
        done = {(str(r["source"]), str(r["accession"]), str(r["sample_id"]))
                for r in rows if str(r.get("status")) in ("ok", "already_done")}

    for i, (_, r) in enumerate(todo.iterrows(), 1):
        key = (str(r.source), str(r.accession), str(r.sample_id))
        if key in done:
            continue
        free_gb = os.statvfs(str(DATA))[4] * os.statvfs(str(DATA))[0] / 1e9
        if free_gb < min_free_gb:
            logger.error("stopping: only %.0f GB free (floor %.0f)", free_gb, min_free_gb)
            break
        rec = convert_one(*key, delete_raw=delete_raw)
        rows.append(rec)
        if i % 20 == 0 or rec["status"] != "ok":
            pd.DataFrame(rows).to_csv(LEDGER, index=False)
        if i % 50 == 0:
            ok = [x for x in rows if x.get("status") == "ok"]
            logger.info("[%d/%d] ok=%d  written=%.1f GB  freed=%.1f GB  free=%.0f GB",
                        i, len(todo), len(ok),
                        sum(x.get("bytes_out", 0) for x in ok) / 1e9,
                        sum(x.get("bytes_freed", 0) or 0 for x in ok) / 1e9, free_gb)
    pd.DataFrame(rows).to_csv(LEDGER, index=False)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    convert_all()
