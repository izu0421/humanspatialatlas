"""Find samples deposited more than once, and pick one representative of each.

Why this exists: the same section reaches us through more than one route. Zenodo mints a new
record id per version, so a re-uploaded dataset arrives as two accessions with consecutive ids
(14941390 and 14941391 are the same 39 MERFISH sections). Authors also mirror a GEO series to
Mendeley or ArrayExpress. Harvesting treats each route as a dataset, and the units then get
counted twice -- 29.2M of them before this ran.

Matching on shape alone would be reckless: two 2203 x 36601 Visium sections from one study are
not the same section. So a candidate pair is only called a duplicate once the expression values
themselves agree, which is why this hashes X for the candidates (a few hundred rows) rather
than for the corpus.

Which copy survives matters for metadata, not for counts: the primary archive usually carries the
richer sample annotation, so GEO/ArrayExpress beats a Zenodo or Mendeley mirror, and among Zenodo
versions the highest id (the newest) wins.
"""
from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import RUNS

logger = logging.getLogger(__name__)

# lower number = preferred as the surviving copy
ARCHIVE_RANK = {"GEO": 0, "ArrayExpress": 1, "E-MTAB": 1, "HCA": 2, "CELLxGENE": 2,
                "Zenodo": 3, "Mendeley": 4, "figshare": 4, "10x Datasets": 5}

HASH_ROWS = 2_000        # first N cells is plenty to separate real sections
HASH_VALS = 200_000


def _data_hash(path: str) -> str | None:
    """MD5 of the leading expression values, for deciding whether two rows are the same data."""
    import anndata as ad
    try:
        a = ad.read_h5ad(path, backed="r")
        X = a.X[: min(HASH_ROWS, a.shape[0])]
        d = X.data[:HASH_VALS] if sp.issparse(X) else np.asarray(X).ravel()[:HASH_VALS]
        return hashlib.md5(np.ascontiguousarray(d, dtype=np.float64).tobytes()).hexdigest()
    except Exception as e:
        logger.warning("hash failed for %s: %s", path, e)
        return None


def _rank(row) -> tuple:
    """Sort key; the smallest wins. Richer metadata first, then primary archive, then newest id."""
    src = str(row.get("source", ""))
    acc = str(row.get("accession", ""))
    meta_gaps = sum(
        1 for c in ("tissue_ontology_term_id", "donor_id", "disease_state", "age_group", "sex")
        if pd.isna(row.get(c)) or str(row.get(c)) in ("", "unknown", "nan")
    )
    # among Zenodo versions the later record is the corrected one, so negate the numeric id
    digits = "".join(ch for ch in acc if ch.isdigit())
    newest = -int(digits) if digits else 0
    return (meta_gaps, ARCHIVE_RANK.get(src, 9), newest, acc)


def find(manifest_csv: str | None = None, out_csv: str | None = None) -> pd.DataFrame:
    """One row per materialised sample that is a redundant copy of another."""
    manifest_csv = manifest_csv or (RUNS / "sample_manifest.csv")
    out_csv = out_csv or (RUNS / "duplicates.csv")
    m = pd.read_csv(manifest_csv, low_memory=False)
    m = m[m.h5ad_path.notna() & m.n_cells.notna()].copy()

    n_vars = pd.to_numeric(m.get("h5ad_n_vars"), errors="coerce")
    m["_key"] = (m.sample_id.astype(str) + "|" + m.n_cells.astype("Int64").astype(str)
                 + "|" + n_vars.astype("Int64").astype(str))

    counts = m._key.value_counts()
    cand = m[m._key.isin(counts[counts > 1].index)].copy()
    logger.info("%d candidate rows in %d shape-identical groups", len(cand), cand._key.nunique())
    if not len(cand):
        pd.DataFrame().to_csv(out_csv, index=False)
        return pd.DataFrame()

    logger.info("hashing expression values to confirm")
    cand["_hash"] = [_data_hash(p) for p in cand.h5ad_path]

    drop = []
    for (_, h), g in cand.groupby(["_key", "_hash"], dropna=True):
        if len(g) < 2:
            continue                                   # same shape, different data: keep both
        if g.accession.nunique() == 1:
            continue                                   # one deposit listing a sample twice is a
                                                       # different problem; not ours to resolve
        ordered = sorted(g.to_dict("records"), key=_rank)
        keep = ordered[0]
        for r in ordered[1:]:
            drop.append({"source": r["source"], "accession": r["accession"],
                         "sample_id": r["sample_id"], "technology": r["technology"],
                         "n_cells": r["n_cells"], "h5ad_path": r["h5ad_path"],
                         "duplicate_of_accession": keep["accession"],
                         "duplicate_of_source": keep["source"], "data_md5": h[:16]})

    d = pd.DataFrame(drop)
    d.to_csv(out_csv, index=False)
    if len(d):
        logger.info("%d redundant copies, %s redundant units", len(d), f"{d.n_cells.sum():,.0f}")
        logger.info("by duplicated-from accession:\n%s",
                    d.duplicate_of_accession.value_counts().head(10).to_string())
    return d


AGG_RATIO = 1.05   # an aggregate holding no more units than its siblings adds nothing


def series_aggregates(manifest_csv: str | None = None, out_csv: str | None = None) -> pd.DataFrame:
    """Rows that are a whole series in one file, alongside that series' per-sample rows.

    GEO series often carry a merged, QC-filtered object in the series-level supplementary files
    as well as per-GSM matrices. Both get harvested, and the merged copy then contributes units a
    second time -- GSE330224 added 590k spots already present in its 319 per-GSM sections, with
    array indices for coordinates instead of pixels, so it was also the worse copy to keep.

    Only flagged when the aggregate holds no more units than its siblings do: an aggregate that
    adds units is contributing sections that were never deposited per-sample, and keeping it is
    then the right call even though it overlaps.
    """
    manifest_csv = manifest_csv or (RUNS / "sample_manifest.csv")
    out_csv = out_csv or (RUNS / "series_aggregates.csv")
    m = pd.read_csv(manifest_csv, low_memory=False)
    m = m[m.h5ad_path.notna()].copy()
    sib_n = m.groupby("accession").sample_id.nunique()

    rows = []
    agg = m[m.sample_id.astype(str) == m.accession.astype(str)]
    for r in agg.itertuples(index=False):
        if sib_n.get(r.accession, 0) <= 1:
            continue                                    # the only row: it IS the sample
        g = m[m.accession == r.accession]
        sibs = g[g.sample_id != r.sample_id]
        su = pd.to_numeric(sibs.n_cells, errors="coerce").fillna(0).sum()
        au = float(pd.to_numeric(pd.Series([r.n_cells]), errors="coerce").fillna(0).iloc[0])
        if su <= 0 or au <= 0:
            continue
        ratio = au / su
        if ratio <= AGG_RATIO:
            rows.append({"source": r.source, "accession": r.accession, "sample_id": r.sample_id,
                         "technology": r.technology, "n_cells": au, "h5ad_path": r.h5ad_path,
                         "n_siblings": len(sibs), "sibling_units": su, "ratio": round(ratio, 3),
                         "reason": "series-level aggregate of its own per-sample rows"})

    d = pd.DataFrame(rows)
    d.to_csv(out_csv, index=False)
    if len(d):
        logger.info("%d series-level aggregates, %s redundant units",
                    len(d), f"{d.n_cells.sum():,.0f}")
    return d


def redundant_keys(csv: str | None = None, agg_csv: str | None = None) -> set[tuple[str, str, str]]:
    """(source, accession, sample_id) of every copy that should not be counted."""
    from pathlib import Path
    csv = csv or (RUNS / "duplicates.csv")
    agg_csv = agg_csv or (RUNS / "series_aggregates.csv")
    out = set()
    for f in (csv, agg_csv):
        if not Path(f).exists():
            continue
        d = pd.read_csv(f)
        if not len(d):
            continue
        out |= {(str(r.source), str(r.accession), str(r.sample_id))
                for r in d.itertuples(index=False)}
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    find()
    series_aggregates()


# ---------------------------------------------------------------------------------------
# Rule 3: the same archive sample under two accessions.
#
# find() keys on sample_id|n_cells|n_vars, so it only catches a copy that was re-deposited
# untouched. It misses the common case: a mirror that reprocessed the data. GSE240429's
# GSM7697868 is 4,992 spots in GEO and 2,378 in the Mendeley BIG-600K mirror -- the same
# tissue section, measured once, counted twice, and invisible to a shape key.
#
# The identity used here is the archive accession embedded in the sample id (GSM, E-MTAB,
# SRS), which survives reprocessing because it names the biological sample rather than the
# matrix. Where the shapes also match, the expression values are hashed as find() does, so an
# exact mirror is confirmed rather than assumed; where they do not, the pair is reported as
# `reprocessed` and the ratio of unit counts is given, because that is a judgement the
# headline depends on: two copies of one section are one section's worth of tissue.
ARCHIVE_ID_PATTERNS = [r"(GSM\d{5,})", r"(E-MTAB-\d+)", r"(E-GEOD-\d+)", r"(SRS\d{6,})"]


def archive_id(sample_id) -> str | None:
    """The archive accession named inside a sample id, if there is one."""
    import re

    s = str(sample_id)
    for pat in ARCHIVE_ID_PATTERNS:
        m = re.search(pat, s, re.I)
        if m:
            return m.group(1).upper()
    return None


def cross_archive(manifest_csv: str | None = None,
                  out_csv: str | None = None) -> pd.DataFrame:
    """One row per section that re-deposits an archive sample already held under another
    accession. The surviving copy is chosen by archive rank, then by unit count."""
    manifest_csv = manifest_csv or (RUNS / "sample_manifest.csv")
    out_csv = out_csv or (RUNS / "duplicates_cross_archive.csv")
    m = pd.read_csv(manifest_csv, low_memory=False)
    m = m[m.h5ad_path.notna() & m.verdict.ne("quarantined") & m.n_cells.notna()].copy()
    m["_aid"] = m.sample_id.map(archive_id)
    g = m.dropna(subset=["_aid"])
    multi = g.groupby("_aid").accession.nunique()
    ids = multi[multi > 1].index
    logger.info("%d sections carry an archive id; %d ids appear under more than one accession",
                len(g), len(ids))
    if not len(ids):
        pd.DataFrame().to_csv(out_csv, index=False)
        return pd.DataFrame()

    drop = []
    for aid, sub in g[g._aid.isin(ids)].groupby("_aid"):
        shapes_match = sub.n_cells.nunique() == 1 and sub.h5ad_n_vars.nunique() == 1
        kind, md5 = "reprocessed", ""
        if shapes_match:
            hashes = {p: _data_hash(p) for p in sub.h5ad_path}
            if len(set(h for h in hashes.values() if h)) == 1:
                kind, md5 = "identical", next(iter(hashes.values()))[:16]
        # _rank() is not reused here. Its third term prefers the HIGHER numeric accession,
        # which is right for Zenodo (a new record id is a new version of the same deposit) and
        # wrong for GEO (GSE320041 and GSE320042 are two series, not two versions). Applied to
        # this rule it kept the smaller copy in 16 of 28 cases, discarding 11,013 units of real
        # data. Here the order is: richer metadata, then primary archive, then MORE units.
        def _pick(r):
            gaps = sum(1 for c in ("tissue_ontology_term_id", "donor_id", "disease_state",
                                   "age_group", "sex")
                       if pd.isna(r.get(c)) or str(r.get(c)) in ("", "unknown", "nan"))
            return (gaps, ARCHIVE_RANK.get(str(r.get("source", "")), 9),
                    -float(r["n_cells"]), str(r["accession"]))

        ordered = sorted(sub.to_dict("records"), key=_pick)
        keep = ordered[0]
        for r in ordered[1:]:
            drop.append({"archive_id": aid, "kind": kind,
                         "source": r["source"], "accession": r["accession"],
                         "sample_id": r["sample_id"], "n_cells": r["n_cells"],
                         "h5ad_path": r["h5ad_path"],
                         "duplicate_of_accession": keep["accession"],
                         "duplicate_of_source": keep["source"],
                         "duplicate_of_path": keep["h5ad_path"],
                         "keep_n_cells": keep["n_cells"],
                         "unit_ratio": round(float(r["n_cells"]) / max(1.0,
                                                                      float(keep["n_cells"])), 3),
                         "data_md5": md5})
    d = pd.DataFrame(drop)
    d.to_csv(out_csv, index=False)
    if len(d):
        logger.info("%d redundant copies (%s), %s redundant units",
                    len(d), d.kind.value_counts().to_dict(), f"{d.n_cells.sum():,.0f}")
    return d


def apply_cross_archive(dry_run: bool = True, csv: str | None = None,
                        manifest_csv: str | None = None) -> pd.DataFrame:
    """Mark the redundant copies in the manifest and remove their .h5ad files.

    Deleting is irreversible, so each row is gated on the SURVIVOR being readable first: a
    duplicate is only removed once the copy that replaces it has been opened and found to hold
    the cells the manifest claims. A survivor that fails that check leaves its duplicate in
    place and is reported, because losing both copies is far worse than carrying one twice.
    """
    import anndata as ad

    csv = csv or (RUNS / "duplicates_cross_archive.csv")
    manifest_csv = manifest_csv or (RUNS / "sample_manifest.csv")
    d = pd.read_csv(csv)
    if not len(d):
        logger.info("nothing to apply")
        return d

    checked = []
    for _, r in d.iterrows():
        ok, why = False, ""
        try:
            a = ad.read_h5ad(r.duplicate_of_path, backed="r")
            ok = a.n_obs > 0 and a.n_vars > 0
            why = "" if ok else f"survivor empty ({a.n_obs}x{a.n_vars})"
        except Exception as e:
            why = f"survivor unreadable: {str(e)[:60]}"
        checked.append({**r.to_dict(), "survivor_ok": ok, "skip_reason": why})
    d = pd.DataFrame(checked)
    bad = d[~d.survivor_ok]
    if len(bad):
        logger.warning("%d duplicates kept because their survivor failed the check", len(bad))
    go = d[d.survivor_ok]

    if dry_run:
        logger.info("DRY RUN: would mark and delete %d sections (%s units); %d held back",
                    len(go), f"{go.n_cells.sum():,.0f}", len(bad))
        return d

    m = pd.read_csv(manifest_csv, low_memory=False)
    drop_paths = set(go.h5ad_path)
    hit = m.h5ad_path.isin(drop_paths)
    m.loc[hit, "verdict"] = "duplicate"
    m.loc[hit, "duplicate_of"] = m.loc[hit, "h5ad_path"].map(
        go.set_index("h5ad_path").duplicate_of_accession)
    m.to_csv(manifest_csv, index=False)
    logger.info("marked %d manifest rows verdict=duplicate", int(hit.sum()))

    freed = 0
    for p in go.h5ad_path:
        try:
            fp = Path(p)
            if fp.exists():
                freed += fp.stat().st_size
                fp.unlink()
        except Exception as e:
            logger.warning("could not delete %s: %s", str(p)[-40:], str(e)[:50])
    logger.info("deleted %d files, freed %.1f GB", len(go), freed / 1e9)
    d.to_csv(csv, index=False)
    return d
