"""Decide a sample's species from the genes it measured, for the ones metadata never stated.

Why this exists: 253 loaded samples carry species "unknown", and because the atlas counts only
human material they were being dropped from every headline figure without appearing anywhere as
an exclusion -- 25M measured units, silently. The depositor's metadata cannot settle it (that is
why they are unknown), but the matrix can: a gene vocabulary is species-specific whether it is
written as Ensembl ids or as symbols.

Two signals, strongest first. An Ensembl id names its species outright (ENSG / ENSMUSG), and is
definitive. Failing that, symbol casing separates the two conventions -- SOX17 against Sox17 --
which is the same test hsa/audit.py uses to keep mouse out of the atlas, and it is reused here
rather than reimplemented so the two cannot drift apart and disagree about one sample.
"""
from __future__ import annotations

import logging
import re

import h5py
import numpy as np
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)

UPPER_SYM = re.compile(r"^[A-Z][A-Z0-9\-\.]{1,14}$")
TITLE_SYM = re.compile(r"^[A-Z][a-z0-9]{1,13}[0-9a-z]$")
MIN_SYMBOLS = 10


def _var_names(h) -> list[str]:
    g = h.get("var")
    if g is None:
        return []
    idx = g.attrs.get("_index", "_index")
    idx = idx.decode() if isinstance(idx, bytes) else str(idx)
    d = g.get(idx)
    if d is None:
        return []
    # An h5ad index arrives in three encodings and only one is a plain array: a categorical is a
    # group of categories+codes, and a nullable string array is a group of values+mask. Reading
    # only the plain case returned nothing for all 216 samples here, which looked like "no gene
    # names" rather than like a reader that did not understand the file.
    if isinstance(d, h5py.Group):
        v = d.get("categories") if "categories" in d else d.get("values")
        if v is None:
            return []
    else:
        v = d
    out = v[: min(4000, v.shape[0])]
    return [x.decode() if isinstance(x, bytes) else str(x) for x in out]


def from_h5ad(path: str) -> tuple[str | None, str]:
    """(species, evidence). None when the genes do not settle it."""
    try:
        with h5py.File(path, "r") as h:
            names = _var_names(h)
            if not names:
                return None, "no var names"
            ensg = sum(1 for x in names if x.upper().startswith("ENSG"))
            ensm = sum(1 for x in names if x.upper().startswith("ENSMUS"))
            if ensg + ensm >= MIN_SYMBOLS:
                frac = ensg / (ensg + ensm)
                if frac >= 0.95:
                    return "homo_sapiens", f"{ensg} ENSG ids ({frac:.0%} of {ensg + ensm})"
                if frac <= 0.05:
                    return "mus_musculus", f"{ensm} ENSMUSG ids"
                return None, f"mixed Ensembl: {frac:.0%} human"

            sym = [x for x in names if not x.upper().startswith(("ENSG", "ENSMUS"))]
            up = sum(1 for x in sym if UPPER_SYM.match(x))
            ti = sum(1 for x in sym if TITLE_SYM.match(x))
            if up + ti < MIN_SYMBOLS:
                return None, f"only {up + ti} usable symbols"
            # same thresholds as hsa/audit.py: a small panel must be purer before it is called
            need = 0.8 if up + ti >= 30 else 0.95
            frac_mouse = ti / (up + ti)
            if frac_mouse > need:
                return "mus_musculus", f"{frac_mouse:.0%} mouse-style of {up + ti} symbols"
            if frac_mouse < (1 - need):
                return "homo_sapiens", f"{1 - frac_mouse:.0%} human-style of {up + ti} symbols"
            return None, f"ambiguous casing: {frac_mouse:.0%} mouse-style"
    except Exception as e:
        return None, f"{type(e).__name__}: {str(e)[:60]}"


def resolve_unknown(write: bool = False, limit: int | None = None) -> pd.DataFrame:
    """Fill species for every loaded sample whose metadata never stated one."""
    from . import db

    sm = pd.DataFrame(db.query(
        "SELECT source, accession, sample_id, species FROM sample_meta"))
    man = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    key = ["source", "accession", "sample_id"]
    for d in (sm, man):
        d[key] = d[key].astype(str)
    j = sm.merge(man[key + ["h5ad_path"]], on=key, how="inner")
    todo = j[(j.species.isna() | j.species.isin(["unknown", ""])) & j.h5ad_path.notna()]
    if limit:
        todo = todo.head(limit)
    logger.info("resolving species for %d samples with a materialised matrix", len(todo))

    rows = []
    for i, r in enumerate(todo.itertuples(index=False), 1):
        sp, ev = from_h5ad(str(r.h5ad_path))
        rows.append({**{k: getattr(r, k) for k in key}, "species": sp, "evidence": ev})
        if i % 50 == 0:
            logger.info("  %d/%d", i, len(todo))

    d = pd.DataFrame(rows)
    if len(d):
        d.to_csv(RUNS / "species_resolved.csv", index=False)
        logger.info("resolved:\n%s", d.species.value_counts(dropna=False).to_string())
        if write:
            n = 0
            for r in d[d.species.notna()].itertuples(index=False):
                db.execute("UPDATE sample_meta SET species=? WHERE source=? AND accession=?"
                           " AND sample_id=?",
                           (r.species, r.source, r.accession, r.sample_id))
                n += 1
            logger.info("wrote species for %d samples", n)
    return d


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    import sys
    resolve_unknown(write="--write" in sys.argv)
