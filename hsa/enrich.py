"""Per-sample facts computed from the matrix itself, rather than from the depositor's notes.

Depositor metadata is thin: on the materialised atlas, sex is stated for 30% of samples, age for
24% and self-reported ethnicity for 2%. Some of that is recoverable from the data itself, and
unlike a parsed label it is available for every sample that loads.

Two of the facts here exist because the obvious version of them is wrong. Sex looks like a simple
XIST-versus-chrY ratio until you notice that the pseudoautosomal genes sit on both chromosomes and
are expressed in every female. And a "sample" is assumed to be one piece of tissue from one donor
until you look at the coordinates and find several.
"""
from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)

# Sex-informative genes, chosen from a measurement on this atlas rather than from a chromosome
# annotation. 91 samples with a confident XIST call and no Y signal were compared with 112 with
# clear Y signal; the medians below are CPM in those two groups.
#
#   male-specific (MSY)      female CPM   male CPM      pseudoautosomal   female CPM   male CPM
#     RPS4Y1                      0.00      172.5         CD99               201.1       322.5
#     DDX3Y                       0.00       35.6         SLC25A6            300.3       315.0
#     EIF1AY                      0.00       33.7         GTPBP6              26.0        33.0
#     KDM5D                       0.00       21.8         AKAP17A             23.1        30.1
#     USP9Y                       0.00       20.1         ASMTL               29.8        34.8
#     NLGN4Y / ZFY / UTY          0.00     7.4-9.6        ZBED1               21.6        30.6
#
# The pseudoautosomal genes sit on both X and Y and are expressed in every female sample at
# 20-300 CPM -- more than every real Y gene combined. They are listed here only so that a panel
# that annotates them as chrY can be refused rather than silently inverted. SRY, AMELY and TBL1Y
# are omitted entirely: they were detected in no sample of either sex, so they add noise only.
XIST = ["XIST", "TSIX"]

# Split by measured discrimination against 800 depositor-labelled samples from 476 studies
# (400 male / 400 female), AUC for predicting male from that gene's CPM alone:
#
#   core   DDX3Y .958  KDM5D .958  ZFY .958  UTY .958  USP9Y .957  EIF1AY .953  NLGN4Y .951
#   extra  PRKY .939   TMSB4Y .937  RPS4Y1 .875
#   PAR    CD99 .564   DHRSX .609   SLC25A6 .521   AKAP17A .506   XG .400
#
# RPS4Y1 is the trap. It carries the strongest male signal of any Y gene (187 CPM against 36 for
# DDX3Y) and is the one most often reached for, but it is also detected in 33% of labelled
# FEMALES -- cross-mapping from RPS4X, a highly expressed ribosomal gene -- which costs it 8
# points of AUC. It stays as supporting evidence and is kept out of the core so it can never
# carry a call by itself.
CHRY_CORE = ["DDX3Y", "KDM5D", "ZFY", "UTY", "USP9Y", "EIF1AY", "NLGN4Y"]
CHRY_EXTRA = ["PRKY", "TMSB4Y", "RPS4Y1"]
CHRY = CHRY_CORE + CHRY_EXTRA

PSEUDOAUTOSOMAL = {"CD99", "XG", "SLC25A6", "ASMTL", "P2RY8", "AKAP17A", "CSF2RA", "IL3RA",
                   "SHOX", "ZBED1", "DHRSX", "GTPBP6", "PPP2R3B", "PLCXD1", "VAMP7", "SPRY3",
                   "WASH6P", "IL9R"}


def sex_from_expression(a, min_cpm: float = 3.0, ratio: float = 5.0, min_y_genes: int = 2):
    """-> (call, evidence dict). Returns ('unknown', ...) when the data cannot answer.

    Works in CPM, not raw counts, so the thresholds mean the same thing on a 300-gene Xenium panel
    and a 36,601-gene Visium slide. Three rules do the work:

    1. Pseudoautosomal genes are never counted. They are on both sex chromosomes and every female
       sample here expresses them at 20-300 CPM; one of them admitted as "Y-linked" outweighs the
       entire male-specific signal and flips the call.
    2. A male call needs at least `min_y_genes` distinct genes from CHRY_CORE, not a count total
       and not any Y gene. Every Y gene is detected in some labelled females (9-33% of them), so
       one gene alone can be cross-mapping; two core genes agreeing cannot be, at these levels.
       RPS4Y1 contributes to the totals but cannot satisfy this requirement on its own.
    3. A marker family absent from the PANEL is unknown, not zero. Comparing against that implied
       zero previously called five whole-transcriptome samples male on 32-81 ambient Y counts.

    The call describes the dominant signal on the slide. A section carrying biopsies from more
    than one donor can mix sexes, so `tissue_pieces` is reported alongside it; treat a sex call on
    a multi-piece section as provisional.
    """
    import scipy.sparse as _sp
    names = {}
    for i, g in enumerate(a.var_names):
        names.setdefault(str(g).upper(), i)

    bad = sorted(PSEUDOAUTOSOMAL & set(names))
    xi = [names[g] for g in XIST if g in names]
    y_genes = [g for g in CHRY if g in names]
    yi = [names[g] for g in y_genes]
    core_genes = [g for g in CHRY_CORE if g in names]
    ev = {"xist_genes_in_panel": len(xi), "y_genes_in_panel": len(yi),
          "pseudoautosomal_excluded": bad}
    if not xi and not yi:
        return "unknown", {**ev, "reason": "no sex-informative gene in the panel"}

    X = a.X
    total = float(X.sum()) or 1.0
    ev["total_counts"] = int(total)
    percpm = 1e6 / total

    def col(j):
        v = X[:, j]
        return float(v.sum())

    y_counts = {g: col(names[g]) for g in y_genes}
    cx = sum(col(j) for j in xi)
    cy = sum(y_counts.values())
    x_cpm, y_cpm = cx * percpm, cy * percpm
    n_y_detected = sum(1 for g, v in y_counts.items() if v > 0 and g in CHRY_CORE)
    ev.update(xist_cpm=round(x_cpm, 2), y_cpm=round(y_cpm, 2),
              core_y_genes_in_panel=len(core_genes), core_y_genes_detected=n_y_detected,
              y_per_gene_cpm={g: round(v * percpm, 2) for g, v in y_counts.items() if v > 0})

    if max(x_cpm, y_cpm) < min_cpm:
        return "unknown", {**ev, "reason": f"no sex-informative signal above {min_cpm:g} CPM"}

    # Only one family measurable: the other side's zero is "not in the panel", not "absent".
    one_sided = not xi or not yi
    if one_sided:
        if yi and y_cpm >= min_cpm and n_y_detected >= min_y_genes:
            return "male", {**ev, "note": "no XIST in panel; call rests on Y genes alone"}
        if xi and x_cpm >= min_cpm:
            return "female", {**ev, "note": "no Y gene in panel; call rests on XIST alone"}
        return "unknown", {**ev, "reason": "only one marker family in the panel and it is too weak"}

    if y_cpm >= min_cpm and y_cpm > ratio * max(x_cpm, 0.1):
        if n_y_detected < min_y_genes:
            return "unknown", {**ev, "reason": f"Y signal rests on {n_y_detected} core gene(s); "
                                               f"{min_y_genes} required because every Y gene is "
                                               "detected in some labelled females (RPS4Y1 in 33%)"}
        return "male", ev
    if x_cpm >= min_cpm and x_cpm > ratio * max(y_cpm, 0.1):
        return "female", ev
    return "unknown", {**ev, "reason": f"XIST and Y signal within {ratio}x -- possible mixed "
                                       "section, tumour X abnormality, or low depth"}


def tissue_pieces(a, k: int = 6, gap_factor: float = 8.0, min_spots: int = 25,
                  min_frac: float = 0.05, max_points: int = 40_000, seed: int = 0) -> dict:
    """How many physically separate tissue pieces are on this section?

    One capture area often holds more than one biopsy -- paired tumour/normal cores, a TMA, or
    several sections from different donors mounted together. Treated as a single sample, those
    give a donor-level call (sex, age, disease) that belongs to no one donor, and a spatial
    statistic computed across the gap between pieces is measuring empty slide.

    Pieces are found by linking spots closer than `gap_factor` times the typical nearest-neighbour
    spacing, which adapts to the platform's own resolution rather than assuming microns.

    A piece must hold at least `min_frac` of the units to count. An absolute floor cannot work
    across platforms that span 400 to 1.9 million units: on cell-resolved data it counted stray
    cells as biopsies, reporting 643 "pieces" for a MERFISH section whose real content was one
    tissue plus debris, and putting 76-88% of Xenium, CosMx and MERFISH samples above one piece
    against 3% of Visium. Small components are still reported, as `n_fragments`.

    Above `max_points` the coordinates are subsampled first. A radius graph over a million-cell
    Xenium section is not merely slow: one such sample held a refresh pass at 99.9% CPU and 40 GB
    for 42 hours without finishing. Whether tissue pieces are separated by empty slide is a
    property of the layout, not of the sampling density, so it survives thinning -- the linking
    radius is rescaled by the new spacing, which the subsample measures for itself.
    """
    out = {"n_pieces": 1, "piece_sizes": None, "piece_method": None}
    xy = a.obsm.get("spatial") if hasattr(a, "obsm") else None
    if xy is None or len(xy) < min_spots * 2:
        return out
    xy = np.asarray(xy, dtype=float)[:, :2]
    xy = xy[np.isfinite(xy).all(1)]
    if len(xy) < min_spots * 2:
        return out
    n_full = len(xy)
    if n_full > max_points:
        rng = np.random.default_rng(seed)
        xy = xy[rng.choice(n_full, max_points, replace=False)]
        min_spots = max(5, int(min_spots * max_points / n_full))
    try:
        from sklearn.neighbors import NearestNeighbors
        from scipy.sparse.csgraph import connected_components
    except Exception:
        return out

    nn = NearestNeighbors(n_neighbors=min(k + 1, len(xy))).fit(xy)
    d, _ = nn.kneighbors(xy)
    spacing = float(np.median(d[:, 1]))
    if not np.isfinite(spacing) or spacing <= 0:
        return out
    graph = nn.radius_neighbors_graph(xy, radius=spacing * gap_factor, mode="connectivity")
    n, lab = connected_components(graph, directed=False)
    sizes = np.bincount(lab)
    floor = max(min_spots, int(min_frac * sizes.sum()))
    big = sizes[sizes >= floor]
    frag = sizes[(sizes < floor) & (sizes >= 3)]
    scale = n_full / len(xy)
    out.update(n_pieces=int(len(big)),
               piece_sizes=[int(round(v * scale)) for v in sorted(big.tolist(), reverse=True)[:12]],
               n_fragments=int(len(frag)),
               frac_in_fragments=round(float(frag.sum() / sizes.sum()), 4) if len(frag) else 0.0,
               piece_method=f"connected components at {gap_factor:g}x median spacing "
                            f"({spacing:.1f} units), pieces >= {min_frac:.0%} of units"
                            + (f", subsampled {n_full}->{len(xy)}" if scale > 1 else ""),
               frac_in_largest=round(float(big.max() / sizes.sum()), 3) if len(big) else None)
    return out


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
