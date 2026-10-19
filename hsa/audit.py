"""Standing integrity checks over the materialised atlas.

These began as one-off investigations after a stuck conversion batch, and each one caught a
defect that no existing check would have noticed:

  present       a ledger row saying "ok" is not evidence of a file; one Xenium sample was
                recorded already_done against a truncated h5ad.
  readable      an existing file is not evidence of a complete one -- the truncated file opened
                at the HDF5 level and failed only when a group was read.
  gene_axis     var_names can silently become positional integers when a GEO series shares one
                features.tsv across sibling GSM directories.
  id_as_gene    a dense CSV read with index_col=0 leaves CosMx's cell_ID looking like a gene; it
                carried up to 85% of the library in affected samples.
  empty_axis    spatial-ATAC deposits arrive shaped (N, 0) and are not transcriptomes.

Checks run off the files themselves rather than off titles or recorded metadata, because a title
cannot tell you what a file contains: a keyword scan for "chromatin|methylat" flagged 45 perfectly
good Visium samples belonging to multi-omic studies.
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
LEDGER = HERE / "runs" / "materialise_ledger.csv"
OUT = HERE / "runs" / "atlas_audit.csv"

ID_LIKE = {"fov", "cell_id", "cell", "cellid", "cell_index", "barcode", "cell_barcode",
           "x", "y", "x_global_px", "y_global_px", "area", "width", "height"}
SYMBOL = re.compile(r"^[A-Za-z][A-Za-z0-9\-\.@_]{0,29}$")
ENSEMBL = re.compile(r"^ENSG\d{11}")
UPPER_SYM = re.compile(r"^[A-Z][A-Z0-9\-]{1,}$")        # SOX17, human
TITLE_SYM = re.compile(r"^[A-Z][a-z]{2,}[0-9A-Za-z\-]*$")  # Sox17, mouse


def _var_names(h: h5py.File) -> np.ndarray:
    g = h["var"]
    o = g[str(g.attrs.get("_index", "_index"))]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]       # nullable-string-array
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def check_file(path: str | Path) -> dict:
    """Open one h5ad and return the problems found, or an empty 'problem' field."""
    p = Path(path)
    rec = {"path": str(p), "n_var": -1, "frac_identifier": np.nan,
           "frac_mouse_like": np.nan, "problem": ""}
    if not p.exists():
        return {**rec, "problem": "missing"}
    try:
        with h5py.File(p, "r") as h:
            v = _var_names(h)
            rec["n_var"] = len(v)
            if len(v) == 0:
                return {**rec, "problem": "empty_gene_axis"}
            low = [str(x).strip().lower() for x in v]
            ids = [x for x in low if x in ID_LIKE]
            rec["frac_identifier"] = len(ids) / len(v)
            if ids:
                return {**rec, "problem": f"identifier_as_gene: {sorted(set(ids))[:4]}"}
            s = v[:500]
            gene_like = np.mean([bool(SYMBOL.match(str(x)) or ENSEMBL.match(str(x))) for x in s])
            if gene_like < 0.5:
                kind = ("positional_integers"
                        if np.mean([str(x).isdigit() for x in s]) > 0.5 else "unrecognised")
                return {**rec, "problem": f"var_axis_not_genes: {kind}"}
            # Species, from the symbols rather than from the recorded metadata. 50 wholly mouse
            # samples reached this human atlas; 48 carried no species field at all and 2 were
            # labelled homo_sapiens, so a check that trusts the metadata would have passed all 50.
            sym = [str(x) for x in v if not str(x).startswith(("ENSG", "ENSMUS"))][:600]
            up = sum(1 for x in sym if UPPER_SYM.match(x))
            ti = sum(1 for x in sym if TITLE_SYM.match(x))
            # A 30-symbol minimum excluded every small panel from the species check, which is
            # where two mouse datasets hid: a 33-gene osmFISH mouse cortex panel got in twice,
            # and a 35-plex mouse CODEX panel with it. Ten symbols are plenty to tell SOX17 from
            # Sox17, and a small panel is required to be purer before it is called.
            if up + ti >= 10:
                need = 0.8 if up + ti >= 30 else 0.95
                frac_mouse = ti / (up + ti)
                rec["frac_mouse_like"] = round(frac_mouse, 3)
                if frac_mouse > need:
                    return {**rec, "problem": f"not_human: {frac_mouse:.0%} mouse-style symbols "
                                              f"({up + ti} symbols examined)"}
                if frac_mouse > 0.2 and up + ti >= 30:
                    return {**rec, "problem": f"mixed_species: {frac_mouse:.0%} mouse-style symbols"}
    except Exception as e:                                          # truncated or corrupt
        return {**rec, "problem": f"unreadable: {type(e).__name__}: {str(e)[:70]}"}
    return rec


def run(write: bool = True) -> pd.DataFrame:
    """Audit every sample the ledger claims to have converted."""
    led = pd.read_csv(LEDGER)
    ok = led[led.status.astype(str).str.startswith("ok")
             | led.status.astype(str).eq("already_done")]
    rows = []
    for i, (_, r) in enumerate(ok.iterrows(), 1):
        rec = check_file(r.path)
        rows.append({**rec, "source": r.source, "accession": r.accession,
                     "sample_id": r.sample_id})
        if i % 2000 == 0:
            logger.info("audited %d/%d", i, len(ok))
    d = pd.DataFrame(rows)
    bad = d[d.problem.astype(str).ne("")]
    logger.info("audited %d samples, %d problems", len(d), len(bad))
    if len(bad):
        for kind, n in bad.problem.str.split(":").str[0].value_counts().items():
            logger.warning("  %s: %d", kind, n)
    if write:
        d.to_csv(OUT, index=False)
    return d
