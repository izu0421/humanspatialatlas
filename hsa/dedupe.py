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

import numpy as np
import pandas as pd
import scipy.sparse as sp

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


def find(manifest_csv: str = "runs/sample_manifest.csv",
         out_csv: str = "runs/duplicates.csv") -> pd.DataFrame:
    """One row per materialised sample that is a redundant copy of another."""
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


def series_aggregates(manifest_csv: str = "runs/sample_manifest.csv",
                      out_csv: str = "runs/series_aggregates.csv") -> pd.DataFrame:
    """Rows that are a whole series in one file, alongside that series' per-sample rows.

    GEO series often carry a merged, QC-filtered object in the series-level supplementary files
    as well as per-GSM matrices. Both get harvested, and the merged copy then contributes units a
    second time -- GSE330224 added 590k spots already present in its 319 per-GSM sections, with
    array indices for coordinates instead of pixels, so it was also the worse copy to keep.

    Only flagged when the aggregate holds no more units than its siblings do: an aggregate that
    adds units is contributing sections that were never deposited per-sample, and keeping it is
    then the right call even though it overlaps.
    """
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


def redundant_keys(csv: str = "runs/duplicates.csv",
                   agg_csv: str = "runs/series_aggregates.csv") -> set[tuple[str, str, str]]:
    """(source, accession, sample_id) of every copy that should not be counted."""
    from pathlib import Path
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
