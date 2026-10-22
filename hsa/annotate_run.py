"""Build every organ's reference, then train and apply a model per (organ, panel).

Split into two phases on purpose. The Census fetch is network-bound and highly variable (blocks
of 5,000 cells took anywhere from 24 seconds to 13 minutes on one organ), while training and
annotation are CPU-bound and local. Running the fetches concurrently and the rest afterwards
keeps the slow phase overlapped instead of serialised behind each panel.
"""
from __future__ import annotations

import argparse
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

import anndata as ad
import pandas as pd

from .annotate import CELL_TECH, ORGAN_MAP, cached_reference, panel_key
from .config import RUNS

logger = logging.getLogger(__name__)


def targets() -> pd.DataFrame:
    m = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    c = m[m.h5ad_path.notna() & m.verdict.ne("quarantined")
          & m.technology_resolved.isin(CELL_TECH)].copy()
    c["organ"] = c.tissue.map(ORGAN_MAP)
    return c.dropna(subset=["organ"])


def _panel_genes(path: str) -> list[str]:
    return [str(g) for g in ad.read_h5ad(path, backed="r").var_names]


def organ_gene_union(c: pd.DataFrame, organ: str, max_sections: int = 40) -> list[str]:
    """Union of the panels in one organ, sampled across sections rather than read exhaustively."""
    sub = c[c.organ == organ]
    # one section per (technology, n_vars) is enough to see every distinct panel in practice
    probe = sub.drop_duplicates(["technology_resolved", "h5ad_n_vars"]).head(max_sections)
    genes: set[str] = set()
    for p in probe.h5ad_path:
        try:
            genes |= set(_panel_genes(p))
        except Exception as e:
            logger.warning("panel read failed %s: %s", str(p)[-40:], str(e)[:60])
    return sorted(genes)


def _fetch(organ: str, genes: list[str]) -> dict:
    import time
    t0 = time.time()
    try:
        a = cached_reference(organ, genes)
        return {"organ": organ, "status": "ok", "n_cells": int(a.n_obs),
                "n_genes": int(a.n_vars), "secs": round(time.time() - t0, 1)}
    except Exception as e:
        return {"organ": organ, "status": f"{type(e).__name__}: {str(e)[:110]}",
                "secs": round(time.time() - t0, 1)}


def prefetch(workers: int = 5, organs: list[str] | None = None) -> pd.DataFrame:
    c = targets()
    todo = organs or sorted(c.organ.unique())
    logger.info("building references for %d organs", len(todo))
    jobs = {}
    for o in todo:
        g = organ_gene_union(c, o)
        logger.info("  %-22s %4d sections, panel union %d genes",
                    o, int((c.organ == o).sum()), len(g))
        jobs[o] = g
    rows = []
    with ProcessPoolExecutor(workers) as ex:
        futs = {ex.submit(_fetch, o, g): o for o, g in jobs.items()}
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            logger.info("[%d/%d] %s %s", len(rows), len(futs), r["organ"], r["status"])
            pd.DataFrame(rows).to_csv(RUNS / "reference_build.csv", index=False)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--organs", nargs="*", default=None)
    a = ap.parse_args()
    prefetch(workers=a.workers, organs=a.organs)
