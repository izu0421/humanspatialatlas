"""Panel registry: which samples actually share a gene space.

Reported panel size is not panel identity. Among 62 surveyed Xenium samples there were 38 distinct
gene sets, and the 164 samples reporting "541 genes" span four different panels that share as
little as 27% of their genes. Grouping by size would pool different assays, so cohorts have to be
keyed on the gene set itself.

Reading gene names does not require the matrix. In every format they sit in one small place, so
this is ~2 orders of magnitude cheaper than load_sample() and can be run over the whole corpus.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .config import RUNS
from .standardise import _find, sample_files

logger = logging.getLogger(__name__)
REGISTRY = RUNS / "panel_registry.csv"
PANELS = RUNS / "panels.json"


# Xenium and CellRanger h5 files carry control codewords alongside real genes: one 541-feature
# file held 300 Gene Expression plus 241 Unassigned/Negative Control entries. Reported n_genes in
# sample_qc.csv therefore includes controls for some samples and not others, which makes panel size
# useless as a grouping key until they are stripped.
CONTROL = re.compile(r"NEGCONTROL|UNASSIGNED|CODEWORD|BLANK|ANTISENSE|DEPRECATED"
                     r"|GENOMICCONTROL|INTERGENIC", re.I)


def _drop_controls(genes) -> list[str]:
    return [g for g in genes if not CONTROL.search(str(g))]


def _h5_gene_names(p: Path) -> list[str]:
    import h5py
    with h5py.File(p, "r") as f:
        if "matrix/features/name" in f:
            names = [x.decode() if isinstance(x, bytes) else str(x)
                     for x in f["matrix/features/name"][:]]
            if "matrix/features/feature_type" in f:     # authoritative when present
                ft = [x.decode() if isinstance(x, bytes) else str(x)
                      for x in f["matrix/features/feature_type"][:]]
                gex = [n for n, t in zip(names, ft) if t == "Gene Expression"]
                if gex:
                    return gex
            return names
        if "matrix/features/id" in f:
            return [x.decode() if isinstance(x, bytes) else str(x)
                    for x in f["matrix/features/id"][:]]
        for g in f.keys():                       # CellRanger v2 layout: /<genome>/gene_names
            for key in (f"{g}/gene_names", f"{g}/genes"):
                if key in f:
                    return [x.decode() if isinstance(x, bytes) else str(x) for x in f[key][:]]
    raise KeyError("no feature names in h5")


def _h5ad_gene_names(p: Path) -> list[str]:
    import h5py
    with h5py.File(p, "r") as f:
        var = f["var"]
        idx = var.attrs.get("_index", "_index")
        idx = idx.decode() if isinstance(idx, bytes) else str(idx)
        d = var[idx]
        if d.dtype.kind in "OS":                 # plain string array
            return [x.decode() if isinstance(x, bytes) else str(x) for x in d[:]]
        cats = var[idx]["categories"][:]         # categorical index
        codes = var[idx]["codes"][:]
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in cats]
        return [cats[c] for c in codes]


def _tsv_gene_names(p: Path) -> list[str]:
    """Features from an mtx trio's features.tsv.

    Note this returns the gene SET, so duplicated symbols collapse. load_sample() instead calls
    var_names_make_unique(), which renames duplicates to GENE-1 and so reports a slightly larger
    var count (33,538 vs 33,514 on one Visium HD sample). For panel identity the set is correct;
    the renamed duplicates are an artefact of needing unique axis labels.
    """
    op = gzip.open if p.name.endswith(".gz") else open
    with op(p, "rt") as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    if not rows:
        raise ValueError("empty features file")
    col = 1 if len(rows[0]) > 1 else 0           # 10x features.tsv: id, name, type
    return [r[col] for r in rows]


def _csv_gene_names(p: Path, max_header: int = 4_000_000) -> list[str]:
    """Only the header line, so a multi-GB expression CSV costs one read."""
    op = gzip.open if p.name.endswith(".gz") else open
    with op(p, "rt", errors="replace") as fh:
        head = fh.readline(max_header)
    sep = "\t" if head.count("\t") > head.count(",") else ","
    cells = [c.strip().strip('"') for c in head.rstrip("\n").split(sep)]
    return [c for c in cells[1:] if c]           # drop the index column label


def gene_names(source: str, accession: str, sample_id: str) -> tuple[list[str], str]:
    """-> (gene names, which reader was used). Raises on failure; callers record the error."""
    paths = sample_files(source, accession, sample_id)
    if not paths:
        raise FileNotFoundError("no files")
    if (p := _find(paths, r"\.h5ad$")):
        return _h5ad_gene_names(p), "h5ad"
    if (p := _find(paths, r"(cell_feature_matrix|filtered_feature_bc_matrix|feature_bc_matrix)\.h5$")):
        return _h5_gene_names(p), "10x_h5"
    if (p := _find(paths, r"matrix\.mtx(\.gz)?$")):
        q = _find(paths, r"(features|genes)\.tsv(\.gz)?$")
        if q is None:
            raise FileNotFoundError("mtx without features file")
        return _tsv_gene_names(q), "mtx_trio"
    if (p := _find(paths, r"(exprmat|cell_by_gene|stdata|_counts?\.|expression)")):
        return _csv_gene_names(p), "csv_matrix"
    raise FileNotFoundError("no recognised matrix")


def panel_id(genes) -> str:
    """Stable id for a gene SET, case- and order-insensitive, controls excluded."""
    norm = sorted({str(g).strip().upper() for g in _drop_controls(genes) if str(g).strip()})
    return hashlib.sha1("\n".join(norm).encode()).hexdigest()[:12], norm


def build(technologies=("Xenium", "Visium HD"), resume: bool = True) -> pd.DataFrame:
    qc = pd.read_csv(RUNS / "sample_qc.csv")
    t = qc[qc.status.astype(str).str.startswith("ok") & qc.technology.isin(technologies)]
    rows, panels = [], {}
    if resume and REGISTRY.exists():
        rows = pd.read_csv(REGISTRY).to_dict("records")
        panels = json.loads(PANELS.read_text()) if PANELS.exists() else {}
    seen = {(r["source"], r["accession"], str(r["sample_id"])) for r in rows}

    for i, (_, r) in enumerate(t.iterrows(), 1):
        key = (r.source, r.accession, str(r.sample_id))
        if key in seen:
            continue
        rec = {"source": r.source, "accession": r.accession, "sample_id": str(r.sample_id),
               "technology": r.technology, "reported_n_genes": r.n_genes}
        try:
            g, reader = gene_names(*key)
            pid, norm = panel_id(g)
            rec.update(reader=reader, n_genes=len(norm), panel_id=pid, status="ok")
            panels.setdefault(pid, norm)
        except Exception as e:
            rec.update(status=f"{type(e).__name__}: {str(e)[:80]}")
        rows.append(rec)
        if i % 100 == 0:
            pd.DataFrame(rows).to_csv(REGISTRY, index=False)
            PANELS.write_text(json.dumps(panels))
            logger.info("[%d/%d] %d panels so far", i, len(t), len(panels))
    df = pd.DataFrame(rows)
    df.to_csv(REGISTRY, index=False)
    PANELS.write_text(json.dumps(panels))
    logger.info("registry: %d samples, %d distinct panels", len(df), len(panels))
    return df


# Visium HD bin size, done properly. Geometry alone cannot settle it: coordinates are in pixels
# for some depositors and microns for others, and the per-sample px/um ratio ranges 0.25-6.17, so
# every pitch-based estimator plateaued near 57% accuracy against the samples that declare a size.
# The grid itself is exact, though -- Space Ranger lays a fixed number of bins over the capture
# area, so array_row's extent is an unambiguous key. Verified 10/10 against declared sizes.
# Visium HD ships in two capture areas, 6.5 mm and 11 mm, and Space Ranger lays a fixed bin count
# over each. Observed grids of 5524/1381/691 rows are the 11 mm slide at 2/8/16 um, not odd
# variants of the 6.5 mm one. Classic Visium's 78 x 128 spot grid also turns up here, in samples
# whose depositor labelled them "Visium HD" -- those are a technology mislabel, not a bin size.
HD_GRID_UM = {3350: 2.0, 1675: 4.0, 838: 8.0, 419: 16.0,          # 6.5 mm capture area
              5500: 2.0, 2750: 4.0, 1375: 8.0, 688: 16.0}         # 11 mm capture area
NOT_HD_GRID = {78, 128, 120, 64}                                  # classic Visium spot grids


def _grid_rows(source: str, accession: str, sample_id: str) -> int | None:
    """Number of bin rows from the positions table, without reading the matrix."""
    from .standardise import _read_coord_table
    p = _find(sample_files(source, accession, sample_id), r"tissue_positions")
    if p is None:
        return None
    df = _read_coord_table(p)
    col = next((c for c in df.columns if "array_row" in str(c).lower()), None)
    if col is None:
        return None
    return int(pd.to_numeric(df[col], errors="coerce").max()) + 1


def _grid_from_coords(source: str, accession: str, sample_id: str):
    """Bins-across from coordinates alone, for h5ad deposits with no positions table.

    Unit-free by construction: on a regular grid the modal step between unique x values is the
    pitch, so extent/pitch is the number of bins across in whatever units the file uses, and that
    count is matched against the known Space Ranger grids -- pixels-vs-microns never arises.

    It must be able to refuse. Continuous coordinates (segmented cells, not bins) have no modal
    step, and an earlier version that took the median of the smallest quartile of steps collapsed
    the pitch to near-zero on those, reporting grids of 3.4 million and 14 million rows against
    cell counts of 700k and 1.2M. The guards below reject anything that is not plausibly a grid.
    """
    import h5py
    import numpy as np
    p = _find(sample_files(source, accession, sample_id), r"\.h5ad$")
    if p is None:
        return None
    with h5py.File(p, "r") as f:
        if "obsm/spatial" not in f:
            return None
        xy = f["obsm/spatial"][:, :2].astype(float)
    n_obs = xy.shape[0]
    best = None
    for j in (0, 1):
        u = np.unique(xy[:, j])
        # a bin grid has few distinct coordinates relative to its bins; continuous data does not
        if u.size < 10 or u.size > max(50, 0.5 * n_obs):
            continue
        d = np.diff(u)
        d = d[d > 0]
        if d.size < 5:
            continue
        # modal step: round to 3 significant figures of the smallest step and take the commonest
        q = np.round(d / d.min()).astype(int)
        step = np.bincount(q[q <= 20]).argmax() if (q <= 20).any() else 1
        pitch = np.median(d[q == step]) if (q == step).any() else np.median(d)
        if not np.isfinite(pitch) or pitch <= 0:
            continue
        n = int(round((u[-1] - u[0]) / pitch)) + 1
        if not (50 <= n <= 6000):            # no Space Ranger grid lies outside this
            continue
        best = n if best is None else max(best, n)
    return best


def grid_geometry(technology: str = "Visium HD") -> pd.DataFrame:
    """Per-sample bin size, from the bin grid rather than from the filename or the coordinates.

    Bin size is the dominant uniformity problem for Visium HD: only 96 of 320 samples declare one
    in the sample id, and the rest span 2 um to >100 um, which includes samples that are not on an
    HD grid at all. Pooling them because they share a gene panel would mix resolutions.
    """
    qc = pd.read_csv(RUNS / "sample_qc.csv")
    t = qc[qc.status.astype(str).str.startswith("ok") & (qc.technology == technology)]
    rows = []
    for _, r in t.iterrows():
        rec = {"source": r.source, "accession": r.accession, "sample_id": str(r.sample_id),
               "technology": r.technology, "n_cells": r.n_cells}
        m = re.search(r"square_(\d{3})um", str(r.sample_id))
        rec["declared_um"] = float(m.group(1)) if m else None
        try:
            g = _grid_rows(r.source, r.accession, str(r.sample_id))
        except Exception as e:
            g, rec["note"] = None, f"{type(e).__name__}"
        if g is None or g <= 1:
            try:
                g = _grid_from_coords(r.source, r.accession, str(r.sample_id))
                rec["grid_from"] = "coords"
            except Exception as e:
                rec["note"] = f"coords:{type(e).__name__}"
        rec["grid_rows"] = g
        # nearest known grid, accepting a few percent of slack for cropped capture areas
        rec["bin_um"] = next((v for k, v in HD_GRID_UM.items() if g and abs(g - k) / k < 0.04), None)
        rec["capture_mm"] = None if rec["bin_um"] is None else (
            11.0 if any(abs(g - k) / k < 0.04 for k in (5500, 2750, 1375, 688)) else 6.5)
        if g in NOT_HD_GRID:
            rec["verdict"] = "not_visium_hd"            # classic Visium spot grid
        elif rec["bin_um"] is not None:
            rec["verdict"] = "hd_grid"
        elif g is None:
            rec["verdict"] = "no_positions_table"
        else:
            rec["verdict"] = "unrecognised_grid"
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(RUNS / "grid_geometry.csv", index=False)
    return df


def coordinate_scale(technology: str = "Visium HD") -> pd.DataFrame:
    """Microns per coordinate unit for every sample, derived from the bin grid.

    225 of 320 Visium HD samples store coordinates in full-resolution image pixels and 92 in
    microns, at image scales whose pixels-per-micron ranges 0.25-6.17 -- so no single conversion
    factor works and the depositor's convention cannot be assumed. The grid settles it per sample:
    a known bin count at a known bin size spans a known physical distance, so

        um_per_unit = (grid_rows * bin_um) / coordinate_extent

    Samples already in microns come out at ~1.0, which is the check that this is right rather than
    a rescaling of nonsense.
    """
    geo = pd.read_csv(RUNS / "grid_geometry.csv")
    qc = pd.read_csv(RUNS / "sample_qc.csv")
    key = ["source", "accession", "sample_id"]
    for d in (geo, qc):
        d[key] = d[key].astype(str)
    m = geo.merge(qc[key + ["x_range", "y_range"]], on=key, how="left")
    x = pd.to_numeric(m.x_range, errors="coerce")
    y = pd.to_numeric(m.y_range, errors="coerce")
    extent = pd.concat([x, y], axis=1).max(axis=1)
    span_um = m.grid_rows * m.bin_um                  # physical width the grid covers
    m["um_per_unit"] = span_um / extent
    m["coord_units"] = np.where(m.um_per_unit.between(0.8, 1.25), "microns",
                        np.where(m.um_per_unit.notna(), "pixels", "unknown"))
    out = m[key + ["technology", "bin_um", "grid_rows", "um_per_unit", "coord_units"]]
    out.to_csv(RUNS / "coordinate_scale.csv", index=False)
    return out


def reclassify(technology: str = "Visium HD") -> pd.DataFrame:
    """Decide each sample's technology from file evidence, not the depositor's label.

    Of 320 samples labelled "Visium HD", only 205 sit on an HD bin grid: 64 are on classic
    Visium's 78 x 128 spot grid, and 49 have one distinct coordinate per cell -- not bins at all,
    with 33 of those reading their coordinates from a Xenium-style cells.parquet. Counting all 320
    as Visium HD overstates the platform and would pool three different assays into one cohort.

    Evidence used, in order of authority: the bin grid, the coordinate convention the loader
    matched, and the panel size. Samples that satisfy none are left unknown rather than guessed.
    """
    geo = pd.read_csv(RUNS / "grid_geometry.csv")
    qc = pd.read_csv(RUNS / "sample_qc.csv")
    key = ["source", "accession", "sample_id"]
    for d in (geo, qc):
        d[key] = d[key].astype(str)
    reg = pd.read_csv(REGISTRY) if REGISTRY.exists() else pd.DataFrame(columns=key + ["n_genes"])
    reg[key] = reg[key].astype(str)
    m = (geo.merge(qc[key + ["n_unique_coords", "coord_source", "loader"]], on=key, how="left")
            .merge(reg[key + ["n_genes"]], on=key, how="left"))
    uniq = pd.to_numeric(m.n_unique_coords, errors="coerce") / m.n_cells

    def call(r, u):
        if r.verdict == "hd_grid":
            return "Visium HD", "bin grid %d um" % r.bin_um
        if r.verdict == "not_visium_hd":
            return "Visium", "classic spot grid (%d rows)" % r.grid_rows
        if u is not None and u > 0.95:                    # one position per cell
            if str(r.coord_source) == "xenium":
                return "Xenium", "cells.parquet, continuous coordinates"
            if pd.notna(r.n_genes) and r.n_genes > 10000:
                return "cell-level WTA", "continuous coordinates, %d genes" % r.n_genes
            return "cell-level (unknown platform)", "continuous coordinates"
        return "unknown", r.verdict

    out = []
    for (_, r), u in zip(m.iterrows(), uniq):
        tech, why = call(r, None if pd.isna(u) else float(u))
        out.append({**{k: r[k] for k in key}, "labelled": technology, "resolved": tech,
                    "evidence": why, "bin_um": r.bin_um, "n_cells": r.n_cells})
    df = pd.DataFrame(out)
    df.to_csv(RUNS / "technology_resolved.csv", index=False)
    return df


def hd_cohorts(min_frac: float = 0.90) -> pd.DataFrame:
    """Visium HD cohorts, stratified by bin size as well as gene panel.

    Gene panel alone is not enough: the 2, 8 and 16 um strata differ 64-fold in bin area, so a
    cohort defined only on shared genes would pool objects that are not comparable. Seeded on the
    most COMMON panel in each stratum -- seeding on the largest instead anchors the cohort to a
    rare full-transcriptome deposit that nothing else covers, and the cohort collapses to 2
    samples. Requires grid_geometry() and reclassify() to have run.
    """
    reg = pd.read_csv(REGISTRY)
    reg = reg[reg.status == "ok"]
    res = pd.read_csv(RUNS / "technology_resolved.csv")
    sets = {k: set(v) for k, v in json.loads(PANELS.read_text()).items()}
    key = ["source", "accession", "sample_id"]
    for d in (reg, res):
        d[key] = d[key].astype(str)
    hd = res[res.resolved == "Visium HD"].merge(reg[key + ["panel_id"]], on=key, how="left")

    rows = []
    for b, g in hd.groupby("bin_um"):
        counts = g.groupby("panel_id").size().sort_values(ascending=False)
        seed = counts.index[0]
        mem = [p for p in counts.index
               if len(sets.get(p, set()) & sets[seed]) / len(sets[seed]) >= min_frac]
        sub = g[g.panel_id.isin(mem)]
        rows.append({"technology": "Visium HD", "bin_um": b, "n_panels": len(mem),
                     "n_samples": len(sub), "n_datasets": sub.accession.nunique(),
                     "n_shared_genes": len(set.intersection(*[sets[p] for p in mem])),
                     "n_bins": int(sub.n_cells.sum()),
                     "n_excluded": int(len(g) - len(sub))})
    df = pd.DataFrame(rows)
    df.to_csv(RUNS / "hd_cohorts.csv", index=False)
    return df


def cohorts(min_frac: float = 0.90, min_samples: int = 20) -> pd.DataFrame:
    """Group samples into cohorts that genuinely share a gene space.

    A panel joins a cohort when it covers at least `min_frac` of the SEED panel's genes. The
    denominator must be the seed, not min(|panel|, |seed|): a 343-gene panel is 90% contained in
    the 18,082-gene Visium HD probe set, and admitting it collapses the cohort's shared gene space
    from 16,789 genes to 1.

    Seeds are the most-used panel per technology above `min_samples`, taken largest-first so a
    cohort is anchored on a full panel rather than a subset of one.
    """
    reg = pd.read_csv(REGISTRY)
    reg = reg[reg.status == "ok"]
    sets = {k: set(v) for k, v in json.loads(PANELS.read_text()).items()}
    rows, assigned = [], set()
    for tech, t in reg.groupby("technology"):
        counts = t.groupby("panel_id").size().sort_values(ascending=False)
        seeds = sorted([p for p in counts.index if counts[p] >= min_samples],
                       key=lambda p: -len(sets[p]))
        for seed in seeds:
            if seed in assigned:
                continue
            mem = [p for p in counts.index
                   if p not in assigned and len(sets[p] & sets[seed]) / len(sets[seed]) >= min_frac]
            if not mem:
                continue
            assigned |= set(mem)
            sub = t[t.panel_id.isin(mem)]
            shared = set.intersection(*[sets[p] for p in mem])
            rows.append({"technology": tech, "cohort_seed": seed, "n_panels": len(mem),
                         "n_samples": len(sub), "n_datasets": sub.accession.nunique(),
                         "n_shared_genes": len(shared), "seed_panel_genes": len(sets[seed])})
    df = pd.DataFrame(rows).sort_values(["technology", "n_samples"], ascending=[True, False])
    df.to_csv(RUNS / "panel_cohorts.csv", index=False)
    return df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    build()
    grid_geometry()
    reclassify()
    print(hd_cohorts().to_string(index=False))
    print(cohorts().to_string(index=False))
