"""Is a sample fit for ligand-receptor analysis? Three independent checks.

LR methods need more than "the file opens". They compute a spatially weighted correlation
between two genes, so each of these can silently invalidate a result rather than fail:

  1. counts     a normalised or residual matrix passes every downstream call and returns
                confident nonsense; 584 files in this atlas are not integer counts
  2. geometry   coordinates must be finite, non-degenerate, and on a plausible scale --
                a kernel bandwidth is derived from nearest-neighbour distance, so a sample
                whose spots are all at one point or spread over a nonsense range produces
                a weight matrix that is either empty or fully connected
  3. panel      the pair must be measurable: a gene absent from the panel is NOT a zero,
                and a sample sharing few pairs with the reference contributes mostly NaN.
                Judged against what the panel COULD support: a 300-gene Xenium panel can
                only ever reach ~40 of 2,699 reference pairs, so a flat threshold marked 839
                targeted samples defective for being targeted. Only a whole-transcriptome
                sample with few measurable pairs indicates a real fault (an identifier dialect).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent.parent
RUNS = HERE / "runs"
OUT = RUNS / "lr_ready.csv"

MIN_UNITS = 200            # a Moran statistic on fewer is not worth computing
MIN_PAIRS = 50             # pairs a whole-transcriptome sample should reach
MIN_PAIRS_TARGETED = 5     # below this even a designed panel cannot support an analysis
MIN_MEDIAN_COUNTS = 10     # below this the per-unit signal is mostly zeros
WTA_GENES = 15_000         # above this the sample is whole-transcriptome, not a panel


# Gene identifiers arrive in several dialects. A barnyard deposit writes "GRCh38_ACTB", an
# Ensembl deposit writes "ENSG00000075624.13" with a version suffix, and a probe panel writes
# "ACTB-v1-prime". Matching the raw string scored 276 samples -- several with 60,286 genes --
# as having ZERO measurable pairs, which is the implausible answer that exposed it.
_PREFIX = re.compile(r"^(GRCH38|GRCH37|HG38|HG19|HUMAN|MM10|GRCM38|GRCM39|MOUSE)[_-]")
_SUFFIX = re.compile(r"(\.\d+|-V\d+(-PRIME|-ALT)?)$")


def _norm(x) -> str:
    s_ = str(x).upper()
    s_ = _PREFIX.sub("", s_)
    return _SUFFIX.sub("", s_)


def _col(g, key):
    o = g[key]
    if isinstance(o, h5py.Group):
        if "categories" in o:                      # categorical encoding
            cats = o["categories"][:]; codes = o["codes"][:]
            cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in cats])
            return cats[codes]
        o = o["values"]
    v = o[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def _names(h, axis):
    g = h[axis]
    o = g[str(g.attrs.get("_index", "_index"))]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def _scale_check(xy, technology, grid_bin_um=None, max_pts: int = 20_000):
    """Nearest-neighbour spacing, and whether it matches a physical unit."""
    from scipy.spatial import cKDTree

    from .coord_units import classify

    rng = np.random.default_rng(0)
    p = xy[rng.choice(len(xy), max_pts, replace=False)] if len(xy) > max_pts else xy
    if len(p) < 50:
        return None, None, "too_few_points"
    d, _ = cKDTree(p).query(p, k=2)
    nn = float(np.median(d[:, 1]))
    unit, _um, verdict = classify(nn, str(technology), grid_bin_um, 1.0)
    # a technology with no lattice (GeoMx ROIs, unknown bin size) cannot be judged this way
    if verdict in ("no_expected_spacing",):
        return nn, None, "ok"
    return nn, unit, verdict


def check_one(path: str, pairs: pd.DataFrame, technology: str = "",
              grid_bin_um=None) -> dict:
    r = {"path": str(path), "n_units": None, "n_genes": None, "n_pairs": None,
         "is_counts": None, "geometry": None, "nn_spacing": None, "coord_unit": None,
         "panel_class": None, "verdict": "", "reasons": ""}
    bad = []
    with h5py.File(path, "r") as h:
        var = _names(h, "var")
        r["n_genes"] = len(var)
        r["n_units"] = len(_names(h, "obs"))

        # 1. counts
        X = h.get("X")
        data = X["data"] if isinstance(X, h5py.Group) else X
        s = np.asarray(data[: min(200_000, data.shape[0])], dtype=float)
        r["is_counts"] = bool(np.all(s >= 0) and np.allclose(s, np.round(s)))
        if not r["is_counts"]:
            bad.append("not integer counts")

        # 2. geometry
        if "obsm/spatial" not in h:
            bad.append("no coordinates"); r["geometry"] = "missing"
        else:
            xy = np.asarray(h["obsm/spatial"][:], dtype=float)[:, :2]
            fin = np.isfinite(xy).all(1)
            uniq = len(np.unique(xy[fin], axis=0)) if fin.any() else 0
            if fin.mean() < 0.5:
                bad.append("coordinates mostly non-finite"); r["geometry"] = "non-finite"
            elif uniq < max(10, 0.1 * fin.sum()):
                bad.append("coordinates degenerate"); r["geometry"] = "degenerate"
            elif not np.all((np.nanmax(xy[fin], 0) - np.nanmin(xy[fin], 0)) > 0):
                bad.append("zero spatial span"); r["geometry"] = "zero span"
            else:
                # Scale, not just finiteness. Four GSE282714 sections spanned 2e9 um -- 2 km --
                # and passed every test above, because nothing here looked at how far apart the
                # units actually sit. Nearest-neighbour spacing is the one quantity the
                # technology pins down, so it is what separates a sound section recorded in
                # centimetres (one multiplication away) from a wrong coordinate column.
                nn, unit, scale_v = _scale_check(xy[fin], technology, grid_bin_um)
                r["nn_spacing"] = nn
                r["coord_unit"] = unit
                if scale_v == "rescale":
                    bad.append(f"coordinates in {unit}, not microns")
                    r["geometry"] = f"wrong unit ({unit})"
                elif scale_v == "spacing_unexplained":
                    bad.append(f"spacing {nn:.3g} matches no unit; coordinate column suspect")
                    r["geometry"] = "implausible scale"
                else:
                    r["geometry"] = "ok"

        # 3. panel coverage against the reference pair list.
        # 595 samples index var by Ensembl ID and keep the symbol in a side column
        # (feature_name / gene_ids / gene_symbol). Matching var_names alone scored every one of
        # them as "0 reference pairs measurable", which is an implausible answer for a
        # 36,601-gene whole-transcriptome sample and the reason this check is worth running.
        present = set(_norm(x) for x in var)
        for col in ("feature_name", "gene_symbols", "gene_symbol", "gene_ids", "symbol"):
            if col in h["var"]:
                try:
                    present |= set(_norm(x) for x in _col(h["var"], col))
                except Exception:
                    pass
        hit = (pairs.gene_a.map(_norm).isin(present)
               & pairs.gene_b.map(_norm).isin(present))
        r["n_pairs"] = int(hit.sum())
        r["panel_class"] = ("whole_transcriptome" if r["n_genes"] > WTA_GENES
                            else "targeted")
        if r["panel_class"] == "whole_transcriptome":
            if r["n_pairs"] < MIN_PAIRS:
                bad.append(f"only {r['n_pairs']} reference pairs measurable from "
                           f"{r['n_genes']} genes; identifiers may not be matching")
        elif r["n_pairs"] < MIN_PAIRS_TARGETED:
            bad.append(f"only {r['n_pairs']} reference pairs measurable")

    if r["n_units"] and r["n_units"] < MIN_UNITS:
        bad.append(f"only {r['n_units']} units")
    r["reasons"] = " | ".join(bad)
    if not bad:
        r["verdict"] = "ready"
    elif any("not microns" in b for b in bad):
        # the geometry is correct, only its unit is mislabelled -- coord_rescale fixes it in place
        r["verdict"] = "fixable_rescale"
    elif all("not integer" in b or "only" in b for b in bad):
        r["verdict"] = "usable_with_care"
    else:
        r["verdict"] = "not_ready"
    return r


def run(pairs_csv: str | Path | None = None, limit: int | None = None) -> pd.DataFrame:
    man = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    led = pd.read_csv(RUNS / "materialise_ledger.csv")
    ok = led[led.status.astype(str).str.startswith("ok")
             | led.status.astype(str).eq("already_done")]
    key = ["source", "accession", "sample_id"]
    j = man.merge(ok[key + ["path"]], on=key).query("usable")
    pairs_csv = pairs_csv or (HERE.parent / "lr_atlas" / "pair_universe.csv")
    pairs = pd.read_csv(pairs_csv)
    pairs = pairs[pairs.label.eq(1)] if "label" in pairs.columns else pairs
    logger.info("checking %d samples against %d reference pairs", len(j), len(pairs))
    if limit:
        j = j.head(limit)
    rows = []
    for i, (_, r) in enumerate(j.iterrows(), 1):
        try:
            rec = check_one(r.path, pairs, technology=r.technology_resolved,
                            grid_bin_um=r.get("grid_bin_um"))
        except Exception as e:
            rec = {"path": r.path, "verdict": "not_ready",
                   "reasons": f"unreadable: {type(e).__name__}"}
        rows.append({**{k: r[k] for k in key},
                     "technology": r.technology_resolved, **rec})
        if i % 1000 == 0:
            logger.info("  %d/%d", i, len(j))
    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    logger.info("%s", d.verdict.value_counts().to_dict())
    return d
