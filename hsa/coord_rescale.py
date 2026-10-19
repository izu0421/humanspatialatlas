"""Put coordinates recorded in the wrong unit onto microns, in place.

Why this exists: `coord_units` separates samples whose coordinates are simply labelled wrong
(a MERFISH section written in centimetres) from samples whose coordinate column was picked wrong
(nearest-neighbour spacing that matches no unit at all). The first group needs one multiplication
and nothing else -- the geometry is already correct, so re-downloading 4 GB to recover it would
be waste. The second group cannot be repaired here at all and is left for re-materialisation.

Everything spatial downstream assumes microns: a SpatialDM bandwidth, a COMMOT distance cutoff and
a neighbour radius are all written in microns, so a section in centimetres silently puts every
cell outside every neighbourhood and the analysis returns nothing rather than failing.

Two safety rules, both learned the hard way. The spacing is re-measured here instead of trusted
from the CSV, because a file re-materialised since that CSV was written would be rescaled against
a stale measurement. And the factor is recorded in `uns` so a second run is a no-op rather than a
second multiplication.
"""
from __future__ import annotations

import logging
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from .config import RUNS
from .coord_units import classify, nn_spacing

logger = logging.getLogger(__name__)
OUT = RUNS / "coord_rescale.csv"


def _already(path: Path) -> float | None:
    """The factor a previous run applied, if any."""
    try:
        with h5py.File(path, "r") as h:
            g = h.get("uns/hsa")
            if g is None or "coord_rescale_factor" not in g:
                return None
            v = g["coord_rescale_factor"][()]
            return float(v if not isinstance(v, bytes) else v.decode())
    except Exception:
        return None


def _apply(path: Path, factor: float, unit: str) -> str:
    with h5py.File(path, "r+") as h:
        if "obsm/spatial" not in h:
            return "no obsm/spatial"
        d = h["obsm/spatial"]
        xy = np.asarray(d[...], dtype=np.float64) * factor
        if d.dtype != np.float64 or d.shape != xy.shape:
            del h["obsm/spatial"]
            h["obsm"].create_dataset("spatial", data=xy)
        else:
            d[...] = xy
        g = h.require_group("uns/hsa")
        g.attrs["encoding-type"] = "dict"; g.attrs["encoding-version"] = "0.1.0"
        for k, v in (("coord_rescale_factor", float(factor)),
                     ("coord_um_per_unit", 1.0),
                     ("coord_coord_units", "microns"),
                     ("coord_original_units", unit)):
            if k in g:
                del g[k]
            if isinstance(v, str):
                g.create_dataset(k, data=v, dtype=h5py.string_dtype(encoding="utf-8"))
            else:
                g.create_dataset(k, data=v)
    return "ok"


def run(units_csv: str | None = None, dry: bool = False, limit: int | None = None) -> pd.DataFrame:
    f = Path(units_csv or (RUNS / "coord_units.csv"))
    u = pd.read_csv(f)
    todo = u[u.verdict.eq("rescale")].copy()
    if limit:
        todo = todo.head(limit)
    logger.info("%d samples flagged for rescale", len(todo))

    rows = []
    for i, r in enumerate(todo.itertuples(index=False), 1):
        p = Path(str(r.path))
        rec = {"accession": r.accession, "sample_id": r.sample_id, "technology": r.technology,
               "path": str(p), "unit_csv": r.unit, "factor_csv": r.um_per_unit}
        if not p.exists():
            rows.append({**rec, "status": "missing"}); continue
        prev = _already(p)
        if prev is not None:
            rows.append({**rec, "status": "already_rescaled", "factor": prev}); continue

        # re-measure rather than trust the CSV: the file may have been rewritten since
        s = nn_spacing(str(p))
        if s is None or s.get("nn") is None:
            rows.append({**rec, "status": "no_coords_now"}); continue
        n_units = getattr(r, "n_finite", None) or s.get("n_finite")
        unit, um, verdict = classify(s["nn"], r.technology, getattr(r, "grid_bin_um", None),
                                     s.get("frac_finite", 1.0),
                                     s.get("span_x"), s.get("span_y"), n_units)
        if verdict != "rescale" or um is None:
            rows.append({**rec, "status": f"recheck_says_{verdict}", "nn_now": s["nn"]}); continue
        # compare the FACTOR, not its printed name: a lattice factor is a computed number
        # ("0.4516 um/unit") that will not render identically between two runs
        try:
            csv_f = float(r.um_per_unit)
            if not (0.9 <= (um / csv_f) <= 1.1):
                rows.append({**rec, "status": f"recheck_factor_differs_{um:.4g}",
                             "nn_now": s["nn"]})
                continue
        except (TypeError, ValueError):
            pass

        rec.update(nn_before=s["nn"], factor=um, unit=unit)
        if dry:
            rows.append({**rec, "status": "dry"}); continue
        try:
            st = _apply(p, float(um), str(unit))
            if st == "ok":
                after = nn_spacing(str(p))
                rec["nn_after"] = after["nn"] if after else None
                _, _, v2 = classify(rec["nn_after"], r.technology,
                                    getattr(r, "grid_bin_um", None), 1.0,
                                    after.get("span_x") if after else None,
                                    after.get("span_y") if after else None, n_units)
                st = "ok" if v2 == "ok" else f"after_rescale_{v2}"
            rows.append({**rec, "status": st})
        except Exception as e:
            rows.append({**rec, "status": f"{type(e).__name__}: {e}"[:90]})
        if i % 25 == 0:
            logger.info("  %d/%d", i, len(todo))
            pd.DataFrame(rows).to_csv(OUT, index=False)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    logger.info("statuses:\n%s", d.status.value_counts().to_string())
    return d


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    import sys
    run(dry="--dry" in sys.argv)
