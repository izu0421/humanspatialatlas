"""Re-decide `counts_recoverable` for the matrices the audit could only call `other`.

The audit labels a matrix by its values: raw_counts, log1p, normalised, scaled, other. `other`
is the bucket for everything it cannot name, and 357 sections landed there or in `scaled` and
were marked counts_recoverable=no. That label is load-bearing -- the ligand-receptor scorer
refuses those sections -- so it is worth deciding per file instead of per bucket.

Testing them showed two different things sharing one label:

  * 177 HuBMAP sections are non-negative, ~73% of values are exactly integer, and dividing each
    GENE by its own smallest non-zero value makes them integral (median fraction 1.000). The
    residues are 1/2, 1/3, 1/4, 1/6 -- the signature of a transcript split between overlapping
    segmented cells. These are counts, fractionally allocated. Nothing was destroyed.
  * the `scaled` matrices carry negative values and no lattice. Those genuinely are destroyed:
    a z-score cannot be inverted without the mean and variance that produced it.

So the verdict becomes one of three, and the evidence for each is written out:
  counts        -- integral, usable as counts
  fractional    -- non-negative and on a per-gene integer lattice; counts with split transcripts
  not_recoverable
"""
from __future__ import annotations

import logging

import h5py
import numpy as np
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)
OUT = RUNS / "counts_recover.csv"
N_VALS = 200_000
N_GENES = 400


def probe(path: str) -> dict:
    """Read a corner of the matrix and decide what kind of numbers it holds."""
    with h5py.File(path, "r") as f:
        X = f["X"]
        sparse = isinstance(X, h5py.Group)
        if sparse:
            data = np.asarray(X["data"][:N_VALS], dtype=float)
            cols = np.asarray(X["indices"][:N_VALS])
        else:
            block = np.asarray(X[: min(400, X.shape[0])], dtype=float)
            data = block.ravel()[:N_VALS]
            cols = np.tile(np.arange(block.shape[1]), block.shape[0])[:N_VALS]
    if not len(data):
        return {"verdict": "empty"}
    frac_int = float(np.mean(np.abs(data - np.round(data)) < 1e-6))
    has_neg = bool((data < 0).any())

    # per-gene lattice: within a gene, is every value an integer multiple of its smallest?
    # data and cols must stay aligned -- filtering values without filtering their column
    # indices compares each value against some other gene's, which inverted this verdict.
    keep = data != 0
    data_nz, cols_nz = data[keep], cols[: len(data)][keep]
    lattice = []
    order = pd.Series(cols_nz).value_counts().head(N_GENES).index
    for g in order:
        v = data_nz[cols_nz == g]
        v = v[v > 0]
        if len(v) < 5:
            continue
        u = v.min()
        lattice.append(float(np.mean(np.abs(v / u - np.round(v / u)) < 1e-5)))
    lat = float(np.median(lattice)) if lattice else np.nan

    if has_neg:
        verdict = "not_recoverable"
    elif frac_int > 0.99:
        verdict = "counts"
    elif np.isfinite(lat) and lat > 0.95:
        verdict = "fractional"
    else:
        verdict = "not_recoverable"
    return {"verdict": verdict, "frac_integer": round(frac_int, 4),
            "has_negative": has_neg, "gene_lattice": None if not np.isfinite(lat)
            else round(lat, 4), "min": float(data.min()), "max": float(data.max())}


def run(limit: int | None = None) -> pd.DataFrame:
    h = pd.read_csv(RUNS / "harmonise_log.csv", low_memory=False)
    todo = h[h.counts_recoverable == "no"].path.tolist()[:limit]
    logger.info("probing %d sections currently marked counts_recoverable=no", len(todo))
    rows = []
    for i, p in enumerate(todo, 1):
        try:
            rows.append({"path": p, **probe(p)})
        except Exception as e:
            rows.append({"path": p, "verdict": f"error: {str(e)[:60]}"})
        if i % 100 == 0:
            logger.info("  %d/%d", i, len(todo))
    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    logger.info("verdicts: %s", d.verdict.value_counts().to_dict())
    return d


def apply(d: pd.DataFrame | None = None) -> None:
    """Write the new verdict back to the harmonise log."""
    d = d if d is not None else pd.read_csv(OUT)
    h = pd.read_csv(RUNS / "harmonise_log.csv", low_memory=False)
    lut = d.set_index("path").verdict
    hit = h.path.isin(lut.index)
    h.loc[hit, "counts_recoverable"] = h.loc[hit, "path"].map(
        lambda p: {"counts": "yes", "fractional": "fractional"}.get(lut.get(p), "no"))
    h.to_csv(RUNS / "harmonise_log.csv", index=False)
    logger.info("harmonise log updated: %s", h.counts_recoverable.value_counts().to_dict())
