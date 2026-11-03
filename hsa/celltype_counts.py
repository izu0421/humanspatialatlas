"""Aggregate the per-section cell-type calls into one table the dashboard can read.

`annotate_apply.py` writes `obs/cell_type` into each section and a one-line summary into
`runs/annotate_log.csv`, but that log records only the DOMINANT type per section. The
composition of the atlas -- how many epithelial cells, how many T cells, across how many
organs -- is not in it, and the dashboard cannot open 1,400 h5ad files to find out.

So the counts are read once here, per (section, cell type), and written to
`runs/celltype_counts.csv`. Reading `obs/cell_type` with h5py takes about 10 ms a section,
so the whole corpus is a few seconds; it is the opening of the files, not the counting, that
would be slow if the page did it.

Two things are carried alongside the counts because a count on its own would be read as more
certain than it is:

  frac_agrees   the fraction of cells whose own call survived its neighbourhood consensus.
                CellTypist's per-cell probability saturates on a targeted panel (97% of cells
                above 0.999 on the first section tested), so it is useless as a confidence;
                whether majority voting overruled the cell is the quantity that moves.
  n_genes_used  the panel the call was made on. A 50-gene panel and a 6,106-gene panel both
                produce a "cell type" column, and the difference belongs next to the number.
"""
from __future__ import annotations

import argparse
import collections
import logging
from pathlib import Path

import h5py
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)
OUT = RUNS / "celltype_counts.csv"
LOG = RUNS / "annotate_log.csv"


def _has_calls(path: str) -> bool:
    """Cheap test: does this section already carry cell-type calls?"""
    try:
        with h5py.File(path, "r") as h:
            return "cell_type" in h.get("obs", {})
    except Exception:
        return False


def _read(path: str) -> tuple[dict, dict] | None:
    """(type -> n cells, section-level facts) for one annotated section, or None."""
    with h5py.File(path, "r") as h:
        if "cell_type" not in h.get("obs", {}):
            return None
        d = h["obs"]["cell_type"]
        # anndata writes a plain string-array here, but a nullable-string-array is a GROUP of
        # values+mask. Assuming the dataset form is what made a corpus-wide scan of var_names
        # return a confident zero for every file in it.
        vals = d["values"][()] if isinstance(d, h5py.Group) else d[()]
        counts = collections.Counter(x.decode() if isinstance(x, bytes) else str(x) for x in vals)
        u = h.get("uns/hsa") or {}
        fact = {k: (u[k][()].decode() if isinstance(u[k][()], bytes) else u[k][()])
                for k in ("celltype_model", "celltype_frac_agrees", "celltype_n_genes_used",
                          "celltype_n_pcs", "celltype_n_types") if k in u}
        fact["n_cells"] = int(len(vals))
    return dict(counts), fact


def build(force: bool = False) -> pd.DataFrame:
    """Every section carrying a cell_type column, taken from DISK rather than from the log.

    The log is incomplete by design: annotate_apply skips a section that already carries
    `uns/hsa/celltype_model` and writes no row for it, so a re-run to retry failures silently
    drops the rows of everything it skipped. After the 2026-11-02 retry that was 37 sections --
    annotated, on disk, and absent from the page -- which is the log under-reporting the atlas,
    not the atlas missing them. The annotation itself is in the files, so the files are the
    source of truth and the log only supplies the organ.
    """
    if not LOG.exists():
        logger.warning("no %s: nothing annotated yet", LOG)
        return pd.DataFrame()
    log = pd.read_csv(LOG)
    organ_of = dict(zip(log.path, log.organ))
    paths = set(log.loc[log.status == "ok", "path"])
    try:
        from .annotate_run import targets
        t = targets()
        for pth, org in zip(t.h5ad_path, t.organ):
            if pth not in paths and _has_calls(pth):
                paths.add(pth)
                organ_of.setdefault(pth, org)
    except Exception as e:                      # the page must not die on a target-list failure
        logger.warning("could not scan targets for unlogged annotations: %s", str(e)[:80])
    paths = sorted(paths)
    rows, skipped = [], []
    for p in paths:
        try:
            got = _read(p)
        except Exception as e:                  # a section being written right now, or moved
            skipped.append((p, str(e)[:50]))
            continue
        if got is None:
            skipped.append((p, "no obs/cell_type"))
            continue
        counts, fact = got
        organ = organ_of.get(p, "unknown")
        for t, n in counts.items():
            rows.append({"path": p, "organ": organ, "cell_type": t, "n": n, **fact})
    d = pd.DataFrame(rows)
    # A per-file try/except turns a total failure into an exit-0 run with an empty table, so the
    # read rate is asserted rather than reported.
    rate = len(set(d.path)) / max(1, len(paths))
    logger.info("read %d of %d annotated sections (%.1f%%), %d (section, type) rows",
                len(set(d.path)) if len(d) else 0, len(paths), 100 * rate, len(d))
    if skipped:
        logger.warning("%d sections skipped, e.g. %s", len(skipped), skipped[:2])
    if rate < 0.9:
        raise SystemExit(f"only {100 * rate:.0f}% of annotated sections could be read -- "
                         "the aggregation is broken, not the corpus")
    d.to_csv(OUT, index=False)
    logger.info("wrote %s", OUT)
    return d


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    d = build(force=ap.parse_args().force)
    if len(d):
        print(d.groupby("cell_type").n.sum().sort_values(ascending=False).head(15).to_string())
