"""Phase 2 of annotation: train a model per (organ, panel) and label every section.

`annotate_run.py` only builds the Census references. Nothing has ever called `train()` and
`annotate_one()` across the corpus, which is why celltype_models/ was empty and no sample
carried a label.

Work is grouped by (organ, panel_key) rather than by section, because a model is a function of
the organ's reference and the section's gene panel, and hundreds of sections share one panel.
Training happens once per group; annotation then reuses it.

The run is resumable and idempotent: every section already carrying `uns/hsa/celltype_model`
is skipped unless --force, and the log is appended after each group, so a run that dies
partway loses at most one group.
"""
from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import anndata as ad
import h5py
import pandas as pd

from .annotate import ORGAN_MAP, annotate_one, panel_key, train
from .annotate_run import targets
from .config import RUNS

logger = logging.getLogger(__name__)
LOG = RUNS / "annotate_log.csv"


def _panel_of(path: str) -> str | None:
    try:
        return panel_key([str(g) for g in ad.read_h5ad(path, backed="r").var_names])
    except Exception as e:
        logger.warning("panel read failed %s: %s", Path(path).name[:40], str(e)[:60])
        return None


def _already(path: str) -> bool:
    try:
        with h5py.File(path, "r") as h:
            return "celltype_model" in h.get("uns/hsa", {})
    except Exception:
        return False


def plan(limit_organ: list[str] | None = None) -> pd.DataFrame:
    c = targets()
    if limit_organ:
        c = c[c.organ.isin(limit_organ)]
    logger.info("resolving panels for %d sections", len(c))
    c = c.assign(panel=[_panel_of(p) for p in c.h5ad_path])
    c = c.dropna(subset=["panel"])
    logger.info("%d sections, %d organs, %d (organ, panel) groups",
                len(c), c.organ.nunique(), c.groupby(["organ", "panel"]).ngroups)
    return c


def run(limit_organ: list[str] | None = None, force: bool = False,
        max_sections: int | None = None) -> pd.DataFrame:
    c = plan(limit_organ)
    done = set()
    if LOG.exists() and not force:
        done = set(pd.read_csv(LOG).path)
    rows, n = [], 0
    groups = sorted(c.groupby(["organ", "panel"]), key=lambda kv: -len(kv[1]))
    for gi, ((organ, pk), sub) in enumerate(groups, 1):
        logger.info("[%d/%d] %s panel %s: %d sections", gi, len(groups), organ, pk, len(sub))
        try:
            genes = [str(g) for g in ad.read_h5ad(sub.h5ad_path.iloc[0],
                                                  backed="r").var_names]
            t0 = time.time()
            train(organ, genes)                      # cached after the first call
            logger.info("    model ready in %.0fs", time.time() - t0)
        except Exception as e:
            logger.warning("    train failed: %s", str(e)[:120])
            for p in sub.h5ad_path:
                rows.append({"path": p, "organ": organ, "status": f"train: {str(e)[:80]}"})
            continue
        for p in sub.h5ad_path:
            if p in done or (not force and _already(p)):
                continue
            try:
                r = annotate_one(p, organ)
                r["status"] = "ok"
            except Exception as e:
                r = {"path": p, "organ": organ, "status": str(e)[:120]}
                logger.warning("    %s: %s", Path(p).name[:38], str(e)[:80])
            rows.append(r)
            n += 1
            if max_sections and n >= max_sections:
                break
        pd.DataFrame(rows).to_csv(LOG, index=False)
        ok = sum(1 for r in rows if r.get("status") == "ok")
        logger.info("    cumulative: %d annotated, %d failed", ok, len(rows) - ok)
        if max_sections and n >= max_sections:
            break
    d = pd.DataFrame(rows)
    d.to_csv(LOG, index=False)
    return d


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--organs", nargs="*", default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--max-sections", type=int, default=None)
    a = ap.parse_args()
    d = run(limit_organ=a.organs, force=a.force, max_sections=a.max_sections)
    if len(d):
        print(d.status.value_counts().head(10).to_string())
