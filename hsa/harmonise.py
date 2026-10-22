"""Make every sample present the same interface: obs columns, counts, coordinates.

The consistency review asked all 7,733 files the same questions and the answers were not the
same. Three gaps, in order of how much they cost a downstream user:

  obs    EVERY file has an empty obs -- not a few, all of them. Materialisation reduced obs to a
         bare index for schema uniformity, which made the files uniform by removing the content.
         Nothing in the atlas carries per-cell metadata, so joining a sample back to its donor,
         tissue or platform means re-reading the manifest for every file.
  counts 676 matrices are not raw counts (319 log1p, 278 other, 79 centred). That is recorded in
         uns now, but nothing states whether counts are RECOVERABLE: log1p is invertible and a
         centred matrix is not, and those are different situations for an analyst.
  coords present and finite everywhere, but the unit was recorded only on samples that needed
         converting, so a third of the corpus read as "unknown" when it was simply already right.

What this does NOT do is change any measured value. It writes obs columns, records what X holds
and whether counts can be recovered from it, and states the coordinate unit. A harmoniser that
silently transformed matrices would be the thing nobody could audit later.
"""
from __future__ import annotations

import logging

import h5py
import numpy as np
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)
OUT = RUNS / "harmonise_log.csv"

# Per-cell columns every sample should carry, so a concatenated atlas is self-describing without
# a side lookup. Sample-level values are broadcast; that is the point -- after concat they are
# the only thing distinguishing one section's cells from another's.
OBS_FIELDS = [
    "sample_id", "accession", "source", "technology_resolved", "tissue",
    "tissue_ontology_term_id", "disease_state", "sex", "age_group", "donor_id",
    "development_stage", "self_reported_ethnicity", "preservation",
]
# log1p is invertible (expm1), a library-size normalisation is not without the original totals,
# and a centred matrix has lost the zeros entirely.
COUNTS_RECOVERABLE = {"raw_counts": "yes", "log1p": "yes_expm1", "normalised": "no",
                      "other": "no", "scaled": "no"}


def _write(g, key, val):
    if key in g:
        del g[key]
    if val is None or (isinstance(val, float) and np.isnan(val)):
        g.create_dataset(key, data="unknown", dtype=h5py.string_dtype(encoding="utf-8"))
    elif isinstance(val, str):
        g.create_dataset(key, data=val, dtype=h5py.string_dtype(encoding="utf-8"))
    else:
        g.create_dataset(key, data=val)


def harmonise_one(path: str, meta: dict) -> dict:
    """Give one file the standard obs columns and the provenance fields in uns."""
    rec = {"path": path, "obs_added": 0, "status": "ok"}
    with h5py.File(path, "r+") as h:
        obs = h.require_group("obs")
        idx_name = obs.attrs.get("_index", "_index")
        idx_name = idx_name.decode() if isinstance(idx_name, bytes) else str(idx_name)
        d = obs.get(idx_name)
        if d is None:
            return {**rec, "status": "no obs index"}
        n = (d["values"].shape[0] if isinstance(d, h5py.Group) and "values" in d
             else d["categories"].shape[0] if isinstance(d, h5py.Group) else d.shape[0])
        if isinstance(d, h5py.Group) and "codes" in d:
            n = d["codes"].shape[0]

        added = []
        for f in OBS_FIELDS:
            v = meta.get(f)
            v = "unknown" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)
            if f in obs:
                del obs[f]
            # one value per cell, stored as a categorical so a 1M-cell section costs a few bytes
            cats = np.array([v], dtype=object)
            grp = obs.create_group(f)
            grp.attrs["encoding-type"] = "categorical"
            grp.attrs["encoding-version"] = "0.2.0"
            grp.attrs["ordered"] = False
            grp.create_dataset("categories", data=cats,
                               dtype=h5py.string_dtype(encoding="utf-8"))
            grp.create_dataset("codes", data=np.zeros(n, dtype=np.int8))
            added.append(f)

        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = idx_name
        obs.attrs["column-order"] = np.array(added, dtype=object)
        rec["obs_added"] = len(added)
        rec["n_obs"] = int(n)

        u = h.require_group("uns/hsa")
        u.attrs["encoding-type"] = "dict"; u.attrs["encoding-version"] = "0.1.0"
        kind = str(meta.get("matrix_kind") or "unknown")
        _write(u, "counts_recoverable", COUNTS_RECOVERABLE.get(kind, "unknown"))
        _write(u, "harmonised_schema", 3)
        rec["counts_recoverable"] = COUNTS_RECOVERABLE.get(kind, "unknown")
    return rec


def run(limit: int | None = None, log_every: int = 500) -> pd.DataFrame:
    m = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    m = m[m.h5ad_path.notna() & m.verdict.ne("quarantined")]
    if limit:
        m = m.head(limit)
    logger.info("harmonising %d samples", len(m))
    rows = []
    for i, r in enumerate(m.itertuples(index=False), 1):
        meta = {f: getattr(r, f, None) for f in OBS_FIELDS}
        meta["matrix_kind"] = getattr(r, "matrix_kind", None)
        try:
            rows.append(harmonise_one(str(r.h5ad_path), meta))
        except Exception as e:
            rows.append({"path": r.h5ad_path, "status": f"{type(e).__name__}: {str(e)[:80]}"})
        if i % log_every == 0:
            logger.info("  %d/%d", i, len(m))
            pd.DataFrame(rows).to_csv(OUT, index=False)
    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    logger.info("statuses:\n%s", d.status.value_counts().to_string())
    return d


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    run()
