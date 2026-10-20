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
import time
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


def _ledger_upsert(rec: dict) -> None:
    """Record one sample's outcome in the ledger.

    convert_all writes the ledger in bulk at the end of a batch, so convert_one -- used to redo a
    single sample -- left no trace at all: the .h5ad appeared on disk while the ledger still said
    `deferred_too_large` or `not_loadable`, and manifest.build() reads h5ad_path FROM the ledger,
    so the sample stayed invisible to every downstream count. That bit twice (the 8.2M-bin pancreas
    and 13 recovered Visium sections) before being fixed here rather than patched by hand again.
    """
    k = ["source", "accession", "sample_id"]
    try:
        led = pd.read_csv(LEDGER) if LEDGER.exists() else pd.DataFrame()
    except Exception:
        led = pd.DataFrame()
    if len(led):
        for c in k:
            led[c] = led[c].astype(str)
        m = np.ones(len(led), dtype=bool)
        for c in k:
            m &= led[c].eq(str(rec[c])).to_numpy()
        led = led[~m]
    led = pd.concat([led, pd.DataFrame([rec])], ignore_index=True)
    led.to_csv(LEDGER, index=False)


def convert_one(source: str, accession: str, sample_id: str, delete_raw: bool = False,
                overwrite: bool = False, record: bool = True) -> dict:
    """Load -> enrich -> write .h5ad -> verify -> (optionally) delete the raw files.

    `record` writes the outcome to the ledger. convert_all passes record=False because it writes
    the whole batch itself; a direct caller wants the ledger updated, since nothing else will.
    """
    rec_out = _convert_one_inner(source, accession, sample_id, delete_raw, overwrite)
    if record:
        try:
            _ledger_upsert(rec_out)
        except Exception as e:
            logger.warning("could not update ledger for %s/%s: %s", accession, sample_id, e)
    return rec_out


def _convert_one_inner(source: str, accession: str, sample_id: str, delete_raw: bool = False,
                       overwrite: bool = False) -> dict:
    from . import standardise as St

    rec = {"source": source, "accession": accession, "sample_id": sample_id, "status": "start"}
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{_tag(source, accession, sample_id)}.h5ad"
    if out.exists() and not overwrite:
        # An existing file is not evidence of a finished one. A Xenium sample was recorded
        # already_done against a truncated h5ad (HDF5 claiming data past EOF, only obs written):
        # the shortcut returned before the read-back check that the write path does perform. Open
        # it; if it will not open, treat it as absent and convert again.
        try:
            import anndata as _ad
            _b = _ad.read_h5ad(out, backed="r")
            if _b.n_obs and _b.n_vars:
                return {**rec, "status": "already_done", "path": str(out),
                        "bytes_out": out.stat().st_size}
            raise ValueError(f"empty: {_b.shape}")
        except Exception as e:
            logger.warning("existing h5ad unreadable (%s), reconverting: %s",
                           str(e)[:70], out.name)
            out.unlink()

    meta = manifest.for_sample(source, accession, sample_id)
    try:
        a, info = St.load_sample(source, accession, sample_id)
    except Exception as e:
        return {**rec, "status": f"load_error: {type(e).__name__}: {e}"[:200]}
    if a is None or not str(info.get("status", "")).startswith("ok"):
        return {**rec, "status": f"not_loadable: {info.get('status')}"[:200]}
    # A sample with no gene axis is not a transcriptome. GSE288178 deposits spatial-ATAC h5ads
    # holding fragment counts and obs-level QC (n_fragment, tsse) with shape (N, 0); these reached
    # the writer and failed there as "scipy.sparse does not support dtype object", which named the
    # symptom and hid the cause. Reject them here with the reason.
    if a.n_vars == 0:
        return {**rec, "status": "not_transcriptomic: no gene axis (0 vars) -- "
                                 "source holds a non-expression modality"}
    if a.n_obs == 0:
        return {**rec, "status": "empty: 0 observations"}

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
                technologies=None, min_free_gb: float = 80.0, workers: int = 3,
                max_cells: int | None = 4_000_000,
                retry_failed: bool = False) -> pd.DataFrame:
    """Stream over every loadable sample, smallest first, stopping if disk runs low."""
    m = pd.read_csv(manifest.MANIFEST) if manifest.MANIFEST.exists() else manifest.build()
    todo = m[m.usable.fillna(False)]
    if technologies:
        todo = todo[todo.technology_resolved.isin(technologies)]
    # Smallest first. Largest-first put a multi-million-cell sample at the head of every batch,
    # so a single batch of 25 spent 2.5 h and 150 GB on one whole-transcriptome CosMx table while
    # 1,500 small samples waited behind it. Ascending order clears the bulk and leaves the few
    # giants at the end, where max_cells defers them for a streaming writer.
    todo = todo.sort_values("n_cells", ascending=True)

    done = set()
    rows = []
    if LEDGER.exists():
        prev = pd.read_csv(LEDGER)
        rows = prev.to_dict("records")
        # Every sample already attempted is skipped, not just the successful ones. Counting only
        # ok/already_done left the failures permanently at the head of the ascending queue, so the
        # batch wrapper re-attempted the same 25 unloadable samples on every pass and the
        # conversion count never moved. retry_failed=True is the way to revisit them after a
        # loader fix.
        statuses = ("ok", "already_done") if retry_failed else None
        done = {(str(r["source"]), str(r["accession"]), str(r["sample_id"]))
                for r in rows
                if statuses is None or str(r.get("status")) in statuses}

    # Sequential, but the caller is expected to run this in short batches from a fresh process
    # (see runs/materialise_loop.sh). One long-lived process reached 158 GB RSS: each sample
    # allocates large sparse arrays and glibc does not hand freed arenas back to the OS. A
    # ProcessPoolExecutor was tried first and deadlocked with its worker blocked in pipe_read,
    # so the batch-and-exit approach is used instead -- less clever and it actually runs.
    keys = [(str(r.source), str(r.accession), str(r.sample_id))
            for _, r in todo.iterrows()
            if (str(r.source), str(r.accession), str(r.sample_id)) not in done]
    cells = {(str(r.source), str(r.accession), str(r.sample_id)): r.n_cells
             for _, r in todo.iterrows()}
    if limit:
        keys = keys[:int(limit)]          # limit the REMAINING work, not the whole list
    logger.info("%d samples to convert this batch", len(keys))
    for i, key in enumerate(keys, 1):
        free_gb = os.statvfs(str(DATA))[4] * os.statvfs(str(DATA))[0] / 1e9
        if free_gb < min_free_gb:
            logger.error("stopping: only %.0f GB free (floor %.0f)", free_gb, min_free_gb)
            break
        n_cells = int(cells.get(key) or 0)
        if max_cells and n_cells > max_cells:
            logger.info("[%d/%d] defer %s (%s cells > cap)", i, len(keys), key[2], f"{n_cells:,}")
            rows.append({"source": key[0], "accession": key[1], "sample_id": key[2],
                         "status": f"deferred_too_large: {n_cells:,} cells"})
            continue
        t0 = time.time()
        rec = convert_one(*key, delete_raw=delete_raw, record=False)
        # Per-sample line: the earlier loop only logged every 25th sample, so a batch stuck on
        # sample 1 looked identical to a batch making progress.
        logger.info("[%d/%d] %s %s (%s cells, %.0fs)", i, len(keys),
                    str(rec.get("status"))[:40], key[2][:40], f"{n_cells:,}", time.time() - t0)
        rows.append(rec)
        if i % 10 == 0:
            pd.DataFrame(rows).to_csv(LEDGER, index=False)
        if i % 25 == 0:
            okr = [x for x in rows if x.get("status") == "ok"]
            logger.info("[%d/%d] ok=%d  written=%.1f GB  freed=%.1f GB  free=%.0f GB",
                        i, len(keys), len(okr),
                        sum(x.get("bytes_out", 0) for x in okr) / 1e9,
                        sum(x.get("bytes_freed", 0) or 0 for x in okr) / 1e9, free_gb)
    out = pd.DataFrame(rows)
    if len(out):
        k = ["source", "accession", "sample_id"]
        pref = ~(out.status.astype(str).str.startswith("ok")
                 | out.status.astype(str).eq("already_done"))
        out = (out.assign(_r=pref.astype(int)).sort_values("_r")
                  .drop_duplicates(k, keep="first").drop(columns="_r"))
    out.to_csv(LEDGER, index=False)
    # The return value is the whole accumulated ledger, so its length never reaches zero and the
    # batch wrapper's "nothing left to do" exit never fired. Carry this batch's size separately.
    out.attrs["n_batch"] = len(keys)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    convert_all()


SEX_SUBSAMPLE = 300_000


def refresh_derived(limit: int | None = None, pieces: bool = True,
                    resume: bool = True) -> pd.DataFrame:
    """Recompute inferred_sex, and count tissue pieces, in already-written .h5ad files.

    Sex is rewritten into the file because it is a per-sample fact the file should carry; the
    piece count goes to a sidecar CSV instead. Adding a new uns key to every file would rewrite
    all 310 GB of the atlas to store one integer per sample, so the manifest joins it on instead.

    A fix to the inference does not reach back into files already materialised, and the marker
    list has now changed twice: first to tell "XIST absent from the panel" from "XIST measured as
    zero", and then to drop the pseudoautosomal genes and demote RPS4Y1 after measuring every
    candidate against 800 depositor-labelled samples. Recomputing needs only the stored matrix, so
    this reads each file rather than reloading anything from the raw deposits.
    """
    import anndata as ad

    if not LEDGER.exists():
        return pd.DataFrame()
    L = pd.read_csv(LEDGER)
    todo = L[(L.status == "ok") & L.path.notna()]
    rows = []
    done = set()
    ck = RUNS / "sex_refresh.csv"
    if resume and ck.exists():
        prev = pd.read_csv(ck)
        rows = prev.to_dict("records")
        done = set(prev.path.astype(str))
        todo = todo[~todo.path.astype(str).isin(done)]
        logger.info("resuming: %d already done, %d to go", len(done), len(todo))
    if limit:
        todo = todo.head(limit)
    for i, (_, r) in enumerate(todo.iterrows(), 1):
        if i % 100 == 1 and i > 1:
            logger.info("[%d/%d] %d changed", i - 1, len(todo),
                        sum(1 for x in rows if x.get("changed")))
            pd.DataFrame(rows).to_csv(RUNS / "sex_refresh.csv", index=False)
        f = Path(r.path)
        if not f.exists():
            continue
        try:
            a = ad.read_h5ad(f, backed="r")
            before = a.uns.get("hsa", {}).get("inferred_sex")
            # The sex call is a ratio of two marker families over the whole library, so a random
            # subsample estimates it as well as the full matrix. Reading every cell of a 27-million
            # bin Visium HD section to compute one ratio cost 40 GB of RSS.
            n = a.n_obs
            if n > SEX_SUBSAMPLE:
                idx = np.sort(np.random.default_rng(0).choice(n, SEX_SUBSAMPLE, replace=False))
                sub = a[idx].to_memory()
            else:
                sub = a.to_memory()
            call, ev = enrich.sex_from_expression(sub)
            if n > SEX_SUBSAMPLE:
                ev["subsampled_cells"] = int(SEX_SUBSAMPLE)
            rec = {"path": str(f), "sex": call, "was": before,
                   "reason": str(ev.get("reason", ""))[:90],
                   "core_y_detected": ev.get("core_y_genes_detected"),
                   "xist_cpm": ev.get("xist_cpm"), "y_cpm": ev.get("y_cpm"),
                   "par_excluded": ";".join(ev.get("pseudoautosomal_excluded", []))[:60]}
            if pieces:
                rec.update({k: v for k, v in enrich.tissue_pieces(sub if n <= SEX_SUBSAMPLE
                                                                  else a).items()
                            if k in ("n_pieces", "frac_in_largest", "piece_sizes", "n_fragments",
                                     "frac_in_fragments")})
            if call == before:
                rows.append({**rec, "changed": False})
                continue
            # only now is the whole file needed, and only for the samples whose call changed
            full = ad.read_h5ad(f)
            full.uns.setdefault("hsa", {})
            full.uns["hsa"]["inferred_sex"] = call
            full.uns["hsa"]["inferred_sex_evidence"] = json.loads(json.dumps(ev, default=str))
            full.write_h5ad(f, compression="gzip")
            del full
            rows.append({**rec, "changed": True})
        except Exception as e:
            rows.append({"path": str(f), "error": f"{type(e).__name__}: {e}"[:120]})
    df = pd.DataFrame(rows)
    df.to_csv(RUNS / "sex_refresh.csv", index=False)
    return df


# kept so existing callers and the cycle do not break
refresh_sex = refresh_derived
