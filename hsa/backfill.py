"""Bring every written .h5ad's uns/hsa block up to the current manifest, in place.

The files were produced over two weeks by code that gained fields as it went, so 16 of the 93
uns/hsa fields are missing from some files -- `preservation` and `multi_piece` exist in 283 of
7,613 because they were added after most files were written. A downstream user reading one file
cannot tell whether a missing field means "not applicable" or "written before we recorded it".

Only the uns group is touched. Re-writing whole AnnData objects to add two metadata fields
would move 330 GB of matrix through memory for no reason; h5py can replace the group in place,
which costs kilobytes per file.
"""
from __future__ import annotations

import logging
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)
OUT = RUNS / "backfill.csv"


def _write_scalar(g, key, val):
    """Store one metadata value the way anndata expects to read it back."""
    if key in g:
        del g[key]
    if val is None or (isinstance(val, float) and np.isnan(val)):
        # h5py.Empty round-trips as an Empty object, not as a null a reader can test. The
        # pipeline already writes "unknown" for an absent label everywhere else, so a missing
        # field reads the same whether it was never curated or curated as unknown -- which is
        # the point: the file should not imply a distinction the metadata cannot support.
        d = g.create_dataset(key, data="unknown", dtype=h5py.string_dtype(encoding="utf-8"))
        d.attrs["encoding-type"] = "string"; d.attrs["encoding-version"] = "0.2.0"
        return
    if isinstance(val, (bool, np.bool_)):
        d = g.create_dataset(key, data=bool(val))
        d.attrs["encoding-type"] = "numeric-scalar"
    elif isinstance(val, (int, np.integer)):
        d = g.create_dataset(key, data=int(val))
        d.attrs["encoding-type"] = "numeric-scalar"
    elif isinstance(val, (float, np.floating)):
        d = g.create_dataset(key, data=float(val))
        d.attrs["encoding-type"] = "numeric-scalar"
    else:
        d = g.create_dataset(key, data=str(val),
                             dtype=h5py.string_dtype(encoding="utf-8"))
        d.attrs["encoding-type"] = "string"
    d.attrs["encoding-version"] = "0.2.0"


def run(limit: int | None = None, dry: bool = False) -> pd.DataFrame:
    man = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    led = pd.read_csv(RUNS / "materialise_ledger.csv")
    ok = led[led.status.astype(str).str.startswith("ok")
             | led.status.astype(str).eq("already_done")]
    key = ["source", "accession", "sample_id"]
    j = man.merge(ok[key + ["path"]], on=key, how="inner")
    if limit:
        j = j.head(limit)
    # the union of fields any sample carries is the schema every file should present
    fields = [c for c in man.columns if c not in ("path",)]
    rows = []
    for i, (_, r) in enumerate(j.iterrows(), 1):
        p = Path(str(r.path))
        if not p.exists():
            rows.append({"path": str(p), "status": "missing"}); continue
        try:
            if dry:
                with h5py.File(p, "r") as h:
                    have = set(h["uns/hsa"].keys()) if "uns/hsa" in h else set()
                rows.append({"path": str(p), "status": "dry",
                             "missing": len(set(fields) - have)})
                continue
            with h5py.File(p, "r+") as h:
                uns = h.require_group("uns")
                if "hsa" in uns:
                    before = set(uns["hsa"].keys())
                else:
                    before = set()
                g = uns.require_group("hsa")
                g.attrs["encoding-type"] = "dict"; g.attrs["encoding-version"] = "0.1.0"
                added = 0
                for c in fields:
                    v = r[c]
                    if isinstance(v, float) and np.isnan(v):
                        v = None
                    _write_scalar(g, c, v)
                    if c not in before:
                        added += 1
                _write_scalar(g, "hsa_schema", 2)
            rows.append({"path": str(p), "status": "ok", "added": added})
        except Exception as e:
            rows.append({"path": str(p), "status": f"{type(e).__name__}: {e}"[:90]})
        if i % 500 == 0:
            logger.info("[%d/%d] backfilled", i, len(j))
            pd.DataFrame(rows).to_csv(OUT, index=False)
    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    logger.info("%d files, %d ok", len(d), int((d.status == "ok").sum()))
    return d
