"""Per-sample facts computed from the matrix itself, rather than from the depositor's notes.

Depositor metadata is thin: sex is stated for 31% of samples, age for 22%. Some of that is
recoverable from the data. Sex in particular is strongly determined by XIST against the
Y-linked genes, and unlike a parsed label it is available for every sample that loads.
"""
from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)

XIST = ["XIST", "TSIX"]
CHRY = ["RPS4Y1", "DDX3Y", "UTY", "USP9Y", "KDM5D", "EIF1AY", "NLGN4Y", "ZFY", "TXLNGY", "TMSB4Y"]


def sex_from_expression(a, min_counts: int = 30, ratio: float = 3.0):
    """-> (call, evidence dict). Returns ('unknown', ...) when the data cannot answer.

    Uses the RATIO of XIST to Y-linked counts, not either one's share of the library. A fraction
    threshold fails on whole-transcriptome platforms: counts spread over 36,601 genes put a clearly
    expressed XIST at 1e-4 of the library, so an absolute floor abstained on most Visium samples
    even though the marker genes were present and detected.

    Abstains when the panel carries no informative gene, or when neither side reaches `min_counts`
    -- a targeted panel without XIST or a Y gene genuinely cannot answer, and saying so is better
    than guessing from noise.
    """
    import scipy.sparse as _sp
    names = {}
    for i, g in enumerate(a.var_names):
        names.setdefault(str(g).upper(), i)
    xi = [names[g] for g in XIST if g in names]
    yi = [names[g] for g in CHRY if g in names]
    ev = {"xist_genes_in_panel": len(xi), "y_genes_in_panel": len(yi)}
    if not xi and not yi:
        return "unknown", {**ev, "reason": "no sex-informative gene in the panel"}
    X = a.X
    tot = lambda idx: float(X[:, idx].sum()) if idx else 0.0
    cx, cy = tot(xi), tot(yi)
    ev.update(xist_counts=int(cx), y_counts=int(cy))
    if max(cx, cy) < min_counts:
        return "unknown", {**ev, "reason": f"fewer than {min_counts} informative counts"}
    if cy >= min_counts and cy > ratio * max(cx, 1.0):
        return "male", ev
    if cx >= min_counts and cx > ratio * max(cy, 1.0):
        return "female", ev
    return "unknown", {**ev, "reason": f"XIST and Y counts within {ratio}x of each other"}


def spatial_stats(a, um_per_unit: float | None = None) -> dict:
    """Physical extent and density, which are comparable across platforms once in microns."""
    out = {}
    xy = a.obsm.get("spatial") if hasattr(a, "obsm") else None
    if xy is None or len(xy) == 0:
        return out
    xy = np.asarray(xy, dtype=float)[:, :2]
    w = float(np.nanmax(xy[:, 0]) - np.nanmin(xy[:, 0]))
    h = float(np.nanmax(xy[:, 1]) - np.nanmin(xy[:, 1]))
    out.update(extent_x=round(w, 1), extent_y=round(h, 1))
    if um_per_unit and np.isfinite(um_per_unit) and um_per_unit > 0:
        wu, hu = w * um_per_unit, h * um_per_unit
        out.update(extent_x_um=round(wu, 1), extent_y_um=round(hu, 1))
        area_mm2 = (wu / 1000.0) * (hu / 1000.0)
        if area_mm2 > 0:
            out.update(area_mm2=round(area_mm2, 3),
                       density_per_mm2=round(a.n_obs / area_mm2, 1))
    return out


def qc_stats(a) -> dict:
    """Depth and sparsity, from counts rather than from whatever normalisation a depositor left."""
    X = a.X
    counts = np.asarray(X.sum(1)).ravel()
    genes = np.asarray((X > 0).sum(1)).ravel()
    return {"total_counts": int(counts.sum()),
            "median_counts_per_unit": float(np.median(counts)) if len(counts) else None,
            "median_genes_per_unit": float(np.median(genes)) if len(genes) else None,
            "mean_counts_per_unit": round(float(counts.mean()), 2) if len(counts) else None,
            "frac_units_under_10_counts": round(float((counts < 10).mean()), 4) if len(counts) else None,
            "n_genes_detected": int((np.asarray((X > 0).sum(0)).ravel() > 0).sum())}
