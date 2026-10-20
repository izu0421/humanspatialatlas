"""Is each materialised matrix actually counts, or something already processed?

Why this matters: nothing downstream announces it. A normalised matrix passes every call and
returns confident nonsense -- SpatialDM silently drops every pair on centred data, a Moran
statistic on log-CPM is not the statistic anyone reported, and a differential test on scaled
values is meaningless. Six HuBMAP sections turned out to be z-scored (one with X.sum() =
-293,004) and were found only because the library crashed on an empty array, which is luck
rather than a check.

Classification, from the values themselves rather than from any label:

  raw_counts        non-negative integers
  log1p             non-negative, non-integer, max below ~20 -- the range log1p leaves
  normalised        non-negative, non-integer, row sums clustered on a target (1e4, 1e6)
  scaled            negative values present: centred or z-scored, unrecoverable here
  other             non-negative, non-integer, fits none of the above

The row-sum test is what separates log1p from library-size normalisation, since both are
non-negative and non-integer and neither announces itself.
"""
from __future__ import annotations

import logging

import h5py
import numpy as np
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)
OUT = RUNS / "counts_audit.csv"
N_VALS = 500_000
N_ROWS = 2_000


def _x(h):
    X = h.get("X")
    if X is None:
        return None, None
    if isinstance(X, h5py.Group):                 # sparse
        return X["data"], X
    return X, None


def classify_one(path: str) -> dict:
    rec = {"path": str(path), "kind": None, "n_vals": 0}
    with h5py.File(path, "r") as h:
        d, grp = _x(h)
        if d is None:
            return {**rec, "kind": "no_X"}
        v = np.asarray(d[: min(N_VALS, d.shape[0])], dtype=np.float64)
        if v.size == 0:
            return {**rec, "kind": "empty"}
        rec["n_vals"] = int(v.size)
        rec["min"] = float(v.min()); rec["max"] = float(v.max())
        rec["frac_negative"] = float((v < 0).mean())
        is_int = bool(np.allclose(v, np.round(v)))
        rec["is_integer"] = is_int

        if rec["frac_negative"] > 0:
            return {**rec, "kind": "scaled"}
        if is_int:
            return {**rec, "kind": "raw_counts"}

        # non-negative and non-integer: log1p or library-size normalised. Row sums tell them
        # apart, so reconstruct a few rows from the CSR structure rather than loading the matrix.
        rs = None
        try:
            if grp is not None and "indptr" in grp:
                ip = np.asarray(grp["indptr"][: N_ROWS + 1], dtype=np.int64)
                if len(ip) > 2:
                    hi = int(ip[-1])
                    dat = np.asarray(grp["data"][:hi], dtype=np.float64)
                    rs = np.array([dat[ip[i]:ip[i + 1]].sum() for i in range(len(ip) - 1)])
            else:
                rs = np.asarray(h["X"][: min(200, h["X"].shape[0])], dtype=np.float64).sum(1)
        except Exception:
            rs = None
        if rs is not None and len(rs) > 5:
            rs = rs[np.isfinite(rs) & (rs > 0)]
            if len(rs) > 5:
                rec["row_sum_median"] = float(np.median(rs))
                rec["row_sum_cv"] = float(np.std(rs) / max(1e-9, np.mean(rs)))
                # a normalisation target leaves every row on the same total
                if rec["row_sum_cv"] < 0.01:
                    return {**rec, "kind": "normalised"}
        if rec["max"] <= 20:
            return {**rec, "kind": "log1p"}
        return {**rec, "kind": "other"}


def run(limit: int | None = None, log_every: int = 500) -> pd.DataFrame:
    m = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    m = m[m.h5ad_path.notna()]
    if limit:
        m = m.head(limit)
    logger.info("auditing %d materialised matrices", len(m))
    rows = []
    for i, r in enumerate(m.itertuples(index=False), 1):
        base = {"source": r.source, "accession": r.accession, "sample_id": r.sample_id,
                "technology": r.technology_resolved, "n_cells": r.n_cells}
        try:
            rows.append({**base, **classify_one(str(r.h5ad_path))})
        except Exception as e:
            rows.append({**base, "path": r.h5ad_path, "kind": f"error: {type(e).__name__}"})
        if i % log_every == 0:
            logger.info("  %d/%d", i, len(m))
            pd.DataFrame(rows).to_csv(OUT, index=False)
    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    logger.info("kinds:\n%s", d.kind.value_counts().to_string())
    return d


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    run()
