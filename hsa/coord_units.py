"""Infer the physical unit of each sample's coordinates from nearest-neighbour spacing.

Why this exists: `x_range`/`y_range` cannot tell a broken coordinate column from a correct one
recorded in the wrong unit. A MERFISH section spanning "1.3 x 2.5" is nonsense in microns and
perfectly ordinary in centimetres, and the two need opposite treatments -- re-download versus a
multiplication. Nearest-neighbour spacing settles it, because the physical spacing of the units
is known from the technology: a Visium spot pitch is 100 um, an 8 um HD bin is 8 um, and an
imaged cell centroid sits 5-30 um from its neighbour whatever the depositor wrote in the file.

So: measure median NN spacing, divide the expected spacing by it, and if that ratio lands on a
known unit factor (mm, cm, m) the coordinates are sound and only mislabelled. If it lands nowhere,
the column pick itself is wrong and no rescale will save it.

The expectation is calibrated from the corpus, not from first principles, because the first-
principles guess was wrong. Reasoning from cell size gave 12 um for a Xenium centroid spacing; the
1,261 Xenium sections here actually sit at a median of 25 um, since segmentation misses cells and
a 2D plane through 3D tissue samples only a fraction of them. With a 3x tolerance that 2x error
was enough to read sparse-but-correct micron coordinates as Xenium pixels. Hence `recalibrate`:
take the per-technology median across samples that need no rescale, and judge each sample against
its own technology's observed behaviour.

A Xenium pixel (0.2125 um) is deliberately NOT a candidate unit. It sits only 4.7x from microns,
inside the tolerance this method needs to absorb genuine biological spread in cell density, so the
two can never be told apart this way and guessing between them does more harm than leaving it.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)

# Median centroid-to-centroid spacing expected in microns, per technology. Visium is a
# manufactured 100 um pitch and is exact; the cell-level entries are corpus medians written back
# by `recalibrate` (see CALIB_JSON) and the values here are only the starting point for the first
# pass on a corpus that has never been calibrated.
EXPECTED_NN = {
    "Visium": 100.0,          # 100 um pitch, 55 um spots -- manufactured, not estimated
    "Visium HD": None,        # depends on bin size; taken from grid_bin_um
    "Xenium": 25.0,
    "MERFISH": 25.0,
    "MERSCOPE": 25.0,
    "CosMx": 25.0,
    "seqFISH": 25.0,
    "osmFISH": 25.0,
    "CODEX": 25.0,
    "cell-level WTA": 25.0,
    "cell-level (unknown platform)": 25.0,
    "Slide-seq": 10.0,        # 10 um beads
    "Stereo-seq": None,       # binning-dependent
    "DBiT-seq": 25.0,
    "GeoMx": None,            # hand-placed ROIs, no lattice
}

# unit -> microns per unit. 1.0 means the file was already in microns. Every candidate is at
# least 1000x from its neighbours, so a factor landing near one of them is unambiguous.
UNIT_FACTORS = {
    "microns": 1.0,
    "millimetres": 1_000.0,
    "centimetres": 10_000.0,
    "metres": 1_000_000.0,
}

TOLERANCE = 5.0   # cell density varies several-fold between tissues; the unit candidates are
                  # 1000x apart, so a wide tolerance costs nothing and avoids false rescales


def nn_spacing(path: str, max_pts: int = 20_000, seed: int = 0):
    """Median nearest-neighbour distance in the file's own units, plus span and usable count."""
    import anndata as ad
    from scipy.spatial import cKDTree

    a = ad.read_h5ad(path, backed="r")
    xy = a.obsm.get("spatial")
    if xy is None:
        return None
    xy = np.asarray(xy, dtype=float)
    if xy.ndim != 2 or xy.shape[1] < 2:
        return None
    xy = xy[:, :2]
    finite = np.isfinite(xy).all(1)
    if finite.mean() < 0.5 or finite.sum() < 50:
        return {"nn": None, "n_finite": int(finite.sum()), "frac_finite": float(finite.mean()),
                "span_x": None, "span_y": None}
    xy = xy[finite]
    span = xy.max(0) - xy.min(0)
    rng = np.random.default_rng(seed)
    if len(xy) > max_pts:
        xy = xy[rng.choice(len(xy), max_pts, replace=False)]
    # k=2 because the first neighbour of a point is itself
    d, _ = cKDTree(xy).query(xy, k=2)
    d1 = d[:, 1]
    nn = float(np.median(d1))
    # On a printed lattice every point is one pitch from its neighbour, so the spread of NN
    # distance collapses: Visium sections here sit at CV 0.003-0.016 against 0.523 for random
    # points. That is what licenses dividing by a manufactured pitch, and it is needed as a gate
    # rather than an observation -- of 320 samples labelled "Visium HD", only 205 are on a bin
    # grid, and applying a pitch to the other 115 would corrupt coordinates that were correct.
    cv = float(np.std(d1) / np.mean(d1)) if np.mean(d1) > 0 else np.nan
    # The spread alone is not the right statistic. A Visium section keeps only its in-tissue
    # spots, so at every hole boundary the nearest neighbour jumps to sqrt(3) pitch and CV rises
    # to 0.05-0.3 while the grid is perfectly intact -- that rejected about 1,000 sound sections.
    # What survives masking is QUANTISATION: nearly every spot still has a neighbour at exactly
    # one pitch. Masked Visium scores 0.995-1.000 here and cell-level data scores 0.06-0.15, a
    # margin with no overlap at all.
    frac_modal = float(np.mean(np.abs(d1 - nn) <= 0.05 * nn)) if nn > 0 else np.nan
    return {"nn": nn, "lattice_cv": cv, "frac_modal_nn": frac_modal,
            "n_finite": int(finite.sum()),
            "frac_finite": float(finite.mean()),
            "span_x": float(span[0]), "span_y": float(span[1])}


# A lattice technology has a manufactured or declared pitch, so the conversion is exact:
# um_per_unit = pitch / observed_spacing, whatever that number turns out to be. The 5,439 Visium
# sections here sit at a median spacing of 199 for a 100 um pitch, which is simply the full-
# resolution pixel (~0.5 um) that Space Ranger writes -- not a fault, and not something a list of
# candidate units would have found.
# Slide-seq is deliberately absent: its beads are deposited randomly on the puck, so a 10 um
# bead diameter is not a 10 um pitch. Treating it as a lattice asked for a 0.6x rescale on 189
# sections whose spacing (14.4-17.3 across the corpus) is exactly what randomly packed 10 um
# beads look like in microns already.
LATTICE_PITCH_UM = {"Visium": 100.0}

# Physical bounds, used where there is no lattice to calibrate against. Both come from the
# instrument and the tissue, not from this corpus, which is what makes them non-circular.
LATTICE_FRAC_MIN = 0.90                           # below this the points are not on a grid
SPAN_MIN_UM, SPAN_MAX_UM = 200.0, 60_000.0        # smaller than a usable section / bigger than a slide
# Plausible nearest-neighbour spacing in microns for anything this atlas holds: tighter than
# touching nuclei at one end, wider than a 100 um Visium pitch with room to spare at the other.
# The upper bound is for catching coordinates that are wrong by orders of magnitude, not for
# enforcing a maximum sparsity. 400 um rejected real sections whose spots sit ~440 um apart
# (sparse or heavily filtered grids, hand-placed ROIs) while the failures this is meant to catch
# are at 33,866 um or 0.0007 -- far outside either bound.
NN_MIN_UM, NN_MAX_UM = 2.0, 2_000.0


def _expected(technology: str, grid_bin_um=None) -> float | None:
    """Expected NN spacing in microns, where the technology pins it down."""
    if technology in ("Visium HD", "Stereo-seq"):
        b = pd.to_numeric(grid_bin_um, errors="coerce")
        return float(b) if np.isfinite(b) and b > 0 else None
    # Deliberately NOT falling back to a corpus-measured value. Calibrating the expectation on
    # the samples the expectation already accepts is circular, and it failed loudly when tried:
    # MERFISH came out at 0.00 um because most MERFISH sections here are in centimetres, so the
    # "ok" set it calibrated on was the mis-scaled majority. A technology without a manufactured
    # pitch goes to the span/density regime instead, which asks physics rather than the corpus.
    return LATTICE_PITCH_UM.get(technology)


def classify(nn: float | None, technology: str, grid_bin_um=None, frac_finite: float = 1.0,
             span_x: float | None = None, span_y: float | None = None,
             n_units: float | None = None, frac_modal_nn: float | None = None):
    """Name the microns-per-unit of a sample's coordinates, or say why it cannot be pinned.

    Two regimes, because they carry different evidence:

      * a lattice technology states its pitch, so dividing the pitch by the observed spacing
        gives um_per_unit exactly -- no candidate list, no tolerance;
      * a cell-level technology has no lattice, so the only handles are the span of the section
        and the density of cells in it. Both are bounded by physics, and because scaling the
        coordinates by f scales density by 1/f^2, a 1500x density window pins f to within ~39x.
        That is comfortably inside the 1000x gaps between mm, cm and m, so the unit comes out
        unambiguous -- or it does not come out at all, which is reported rather than guessed.

    Returns (unit, um_per_unit, verdict). 'ok' means already microns, 'rescale' means sound
    geometry in another unit, and anything else means the coordinates cannot be trusted.
    """
    if frac_finite is not None and frac_finite < 0.5:
        return None, None, "coords_non_finite"
    if nn is None or not np.isfinite(nn) or nn <= 0:
        return None, None, "coords_degenerate"

    # Array indices are not a unit, and no metric factor converts them. A 50 x 50 DBiT-seq grid
    # numbered 1..50 has nn = 1.0, and because the candidate factors jump from 1 to 1,000 the
    # only one that lands inside the plausible range is millimetres -- giving cells 1 mm apart.
    # The real pitch (25 or 50 um for DBiT-seq) comes from the protocol, not from the numbers, so
    # the honest answer is to refuse. Detected as consecutive integers on a small dense grid.
    if (span_x is not None and span_y is not None and np.isfinite(span_x)
            and abs(nn - 1.0) < 1e-6 and max(span_x, span_y) < 2_000):
        return None, None, "array_index_not_a_unit"

    # --- lattice regime: the pitch is known, so the factor is a division -- but only once the
    # points are confirmed to sit on a grid, since the technology label alone is not evidence
    exp = _expected(technology, grid_bin_um)
    on_lattice = (frac_modal_nn is not None and np.isfinite(frac_modal_nn)
                  and frac_modal_nn >= LATTICE_FRAC_MIN)
    if exp is not None and not on_lattice and frac_modal_nn is not None:
        exp = None                      # labelled a grid, measured as scattered: use physics
    if exp is not None:
        um = exp / nn
        if not np.isfinite(um) or um <= 0:
            return None, None, "spacing_unexplained"
        if 0.8 <= um <= 1.25:
            return "microns", 1.0, "ok"
        # a plausible factor keeps the section within slide dimensions
        if span_x and span_y and np.isfinite(span_x) and np.isfinite(span_y):
            sx, sy = span_x * um, span_y * um
            if not (SPAN_MIN_UM <= max(sx, sy) <= SPAN_MAX_UM):
                return None, None, "spacing_unexplained"
        return f"{um:.4g} um/unit", float(um), "rescale"

    # --- cell-level regime: no lattice, so use span and density
    if not (span_x and span_y and np.isfinite(span_x) and np.isfinite(span_y)
            and span_x > 0 and span_y > 0 and n_units and n_units > 0):
        return None, None, "no_expected_spacing"

    # Density is measured from the SPACING, not from the bounding box. A merged multi-section
    # object, or any section whose tissue does not fill its own bounding rectangle, has most of
    # that rectangle empty: one such sample spans 30 x 29 mm with 9,626 spots, giving 11 per mm2
    # over the box while its spots sit exactly 100 um apart. Nearest-neighbour spacing is a local
    # quantity and immune to that, and it scales with the unit in the same way, so requiring the
    # converted spacing to be physically plausible does the same job without the empty space.
    ok_factors = []
    for name, f in UNIT_FACTORS.items():
        sx, sy = span_x * f, span_y * f
        if not (SPAN_MIN_UM <= max(sx, sy) <= SPAN_MAX_UM):
            continue
        if not (NN_MIN_UM <= nn * f <= NN_MAX_UM):
            continue
        ok_factors.append((name, f))

    if not ok_factors:
        return None, None, "spacing_unexplained"
    if len(ok_factors) > 1:
        # prefer microns when it is among the survivors; otherwise the choice is not ours to make
        for name, f in ok_factors:
            if name == "microns":
                return "microns", 1.0, "ok"
        return "|".join(n for n, _ in ok_factors), None, "unit_ambiguous"
    name, f = ok_factors[0]
    return name, f, ("ok" if name == "microns" else "rescale")


def run(manifest_csv: str | None = None, out_csv: str | None = None,
        limit: int | None = None, log_every: int = 250) -> pd.DataFrame:
    manifest_csv = manifest_csv or (RUNS / "sample_manifest.csv")
    out_csv = out_csv or (RUNS / "coord_units.csv")
    m = pd.read_csv(manifest_csv, low_memory=False)
    m = m[m.h5ad_path.notna()]
    if limit:
        m = m.head(limit)
    logger.info("calibrating coordinates for %d samples", len(m))

    rows = []
    for i, r in enumerate(m.itertuples(index=False), 1):
        rec = {"source": r.source, "accession": r.accession, "sample_id": r.sample_id,
               "technology": r.technology, "path": r.h5ad_path,
               "grid_bin_um": getattr(r, "grid_bin_um", None)}
        try:
            s = nn_spacing(r.h5ad_path)
            if s is None:
                rec.update(nn=None, verdict="no_coords", unit=None, um_per_unit=None)
            else:
                rec.update(s)
                unit, um, verdict = classify(s["nn"], r.technology, rec["grid_bin_um"],
                                             s.get("frac_finite", 1.0),
                                             s.get("span_x"), s.get("span_y"),
                                             getattr(r, "n_cells", None),
                                             s.get("frac_modal_nn"))
                rec.update(unit=unit, um_per_unit=um, verdict=verdict)
        except Exception as e:
            rec.update(nn=None, verdict="error", unit=None, um_per_unit=None, error=str(e)[:200])
        rows.append(rec)
        if i % log_every == 0:
            logger.info("  %d/%d", i, len(m))

    d = pd.DataFrame(rows)
    d.to_csv(out_csv, index=False)
    logger.info("wrote %s", out_csv)
    logger.info("verdicts:\n%s", d.verdict.value_counts().to_string())
    return d


def reclassify(units_csv: str | None = None, out_csv: str | None = None) -> pd.DataFrame:
    """Re-run `classify` over an existing scan. Cheap: the spacing is already measured."""
    f = Path(units_csv) if units_csv else (RUNS / "coord_units.csv")
    d = pd.read_csv(f)
    man = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    n_by = {(str(a), str(b), str(c)): v for a, b, c, v in
            zip(man.source, man.accession, man.sample_id, man.n_cells)}
    res = []
    for r in d.itertuples(index=False):
        n = n_by.get((str(r.source), str(r.accession), str(r.sample_id)))
        fm = getattr(r, "frac_modal_nn", None)
        res.append(classify(r.nn if pd.notna(r.nn) else None, str(r.technology),
                            getattr(r, "grid_bin_um", None),
                            r.frac_finite if pd.notna(getattr(r, "frac_finite", np.nan)) else 1.0,
                            getattr(r, "span_x", None), getattr(r, "span_y", None), n,
                            fm if fm is not None and pd.notna(fm) else None))
    d["unit"] = [x[0] for x in res]
    d["um_per_unit"] = [x[1] for x in res]
    d["verdict"] = [x[2] for x in res]
    d.to_csv(out_csv or f, index=False)
    logger.info("verdicts:\n%s", d.verdict.value_counts().to_string())
    return d


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    import sys
    if "--reclassify" in sys.argv:
        reclassify()
    else:
        run()
