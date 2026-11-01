"""Resolve the coordinate units that classify() could not name, by testing a pixel hypothesis.

161 sections came back `unknown`. They are not ambiguous -- they are in PIXELS, and classify()
has no pixel hypothesis, so they fell through. A CosMx section reads a span of 106,000 units,
which as microns is 10.6 cm and is rejected; at 0.1202 um/px it is 12.8 mm, which is a slide.

The test is physical and has to pass on BOTH axes at once, which is what makes it safe: a wrong
scale factor can put the cell spacing in range or the section span in range, but not both,
because they differ by three orders of magnitude and a single factor moves them together. A
section is converted only if, at the manufacturer's pixel size, neighbouring cells sit 2-30 um
apart AND the section spans 0.5-40 mm.

Nearest-neighbour spacing is measured on the full coordinate set here rather than a subsample.
The subsampled estimate is inflated by sqrt(N/n) -- for a 1.8M-cell CosMx section that is a
factor of 9.6, enough to push a correct conversion outside the gate.
"""
from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .config import RUNS

logger = logging.getLogger(__name__)
OUT = RUNS / "coord_pixels.csv"

# manufacturer pixel sizes, um per pixel
PIXEL_UM = {"CosMx": 0.1202, "Xenium": 0.2125, "Visium HD": 0.2125,
            "MERFISH": 0.108, "Visium": 0.2125}
NN_MIN_UM, NN_MAX_UM = 2.0, 30.0
SPAN_MIN_MM, SPAN_MAX_MM = 0.5, 40.0
MAX_PTS = 300_000          # full set where it fits; a cap only for the largest sections


def measure(path: str) -> dict:
    a = ad.read_h5ad(path, backed="r")
    xy = np.asarray(a.obsm["spatial"], dtype=float)[:, :2]
    xy = xy[np.isfinite(xy).all(1)]
    if len(xy) < 50:
        return {"n": len(xy)}
    sub = xy if len(xy) <= MAX_PTS else xy[np.random.default_rng(0).choice(
        len(xy), MAX_PTS, replace=False)]
    nn = float(np.median(cKDTree(sub).query(sub, k=2)[0][:, 1]))
    # a subsample thins the plane, so the spacing it reports is too large by sqrt(N/n)
    nn *= np.sqrt(len(sub) / len(xy))
    return {"n": len(xy), "nn_raw": nn,
            "span_x_raw": float(np.ptp(xy[:, 0])), "span_y_raw": float(np.ptp(xy[:, 1]))}


def plan(limit: int | None = None) -> pd.DataFrame:
    r = pd.read_csv(RUNS / "review_samples.csv", low_memory=False)
    unk = r[r.coord_coord_units.astype(str).isin(["unknown", "nan"])]
    logger.info("%d sections with unresolved coordinate units", len(unk))
    rows = []
    for i, (p, tech) in enumerate(zip(unk.path, unk.technology_resolved_y), 1):
        rec = {"path": p, "technology": tech, "px_um": PIXEL_UM.get(str(tech))}
        try:
            rec |= measure(p)
        except Exception as e:
            rec["verdict"] = f"error: {str(e)[:50]}"
            rows.append(rec)
            continue
        px = rec.get("px_um")
        if not px or "nn_raw" not in rec:
            rec["verdict"] = "no pixel size for this platform" if not px else "too few points"
        else:
            nn_um = rec["nn_raw"] * px
            span_mm = max(rec["span_x_raw"], rec["span_y_raw"]) * px / 1000
            rec |= {"nn_um": round(nn_um, 2), "span_mm": round(span_mm, 2)}
            ok_nn = NN_MIN_UM <= nn_um <= NN_MAX_UM
            ok_sp = SPAN_MIN_MM <= span_mm <= SPAN_MAX_MM
            rec["verdict"] = ("pixels" if (ok_nn and ok_sp) else
                              f"rejected (nn {'ok' if ok_nn else 'bad'}, "
                              f"span {'ok' if ok_sp else 'bad'})")
        rows.append(rec)
        if i % 40 == 0:
            logger.info("  %d/%d", i, len(unk))
        if limit and i >= limit:
            break
    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    logger.info("verdicts: %s", d.verdict.map(lambda v: str(v).split(" (")[0])
                .value_counts().to_dict())
    return d


def apply(dry_run: bool = True) -> pd.DataFrame:
    from .coord_rescale import _apply
    from pathlib import Path

    d = pd.read_csv(OUT)
    go = d[d.verdict == "pixels"]
    logger.info("%s %d sections from pixels to microns",
                "DRY RUN: would convert" if dry_run else "converting", len(go))
    if dry_run:
        return go
    out = []
    for _, r in go.iterrows():
        status = _apply(Path(r.path), float(r.px_um), f"pixels ({r.px_um} um/px)")
        out.append({**r.to_dict(), "status": status})
    res = pd.DataFrame(out)
    logger.info("status: %s", res.status.value_counts().to_dict())
    res.to_csv(RUNS / "coord_pixels_applied.csv", index=False)
    return res
