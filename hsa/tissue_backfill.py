"""Fill in tissue for sections the harvest left as `unknown`, by rules that leave evidence.

All 270 remaining unknowns sit in 37 accessions, so this is a per-dataset metadata gap rather
than a per-sample one, and it is fixed per dataset. Every rule records what it matched on, so a
wrong call can be found and reverted rather than becoming an unattributable fact in the atlas.

Two rules, in order of how much they assume:

  token_join  -- the section's id contains a cohort id (ACH0010, 6332STDY10289520, UKF296) that
                 another section in the atlas already carries WITH a resolved tissue. This is an
                 internal join: nothing is fetched and nothing is guessed. It is also a
                 duplicate signal, flagged for review, since the same cohort id in two
                 accessions usually means the same material deposited twice -- but the ids are
                 lab-internal rather than archive accessions, so they are not acted on by
                 dedupe.cross_archive() without a human looking.

  title_rule  -- the dataset title names the tissue unambiguously. Only titles where the tissue
                 is the subject are listed; "colorectal cancer metastases" is NOT one of them,
                 because the metastases are in liver or lung and the deposit does not say which
                 per sample. A guess there would be a coin flip recorded as a fact.
"""
from __future__ import annotations

import logging
import re

import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)
LOG = RUNS / "tissue_backfill.csv"

UNKNOWN = {"unknown", "nan", "none", ""}
TOKEN_RE = r"(UKF\d+|STDY\d+|ACH\d+|10xGenomics\d+|10X\d{3,}|GSM\d{5,})"

# accession -> (tissue, why). Only where the title settles it.
TITLE_RULES = {
    "zenodo:21782884": ("brain", "title: pan-cancer metastatic BRAIN tumors -- the "
                                 "metastases are in brain, whatever the primary was"),
}


def _is_unknown(s) -> bool:
    return str(s).strip().lower() in UNKNOWN


def _tokens(s) -> set[str]:
    return {t.upper() for t in re.findall(TOKEN_RE, str(s), re.I)}


def plan(manifest_csv: str | None = None) -> pd.DataFrame:
    manifest_csv = manifest_csv or (RUNS / "sample_manifest.csv")
    m = pd.read_csv(manifest_csv, low_memory=False)
    live = m[m.h5ad_path.notna() & ~m.verdict.isin(["quarantined", "duplicate"])]
    known = live[~live.tissue.map(_is_unknown)]
    unk = live[live.tissue.map(_is_unknown)]
    logger.info("%d sections without a tissue, in %d accessions",
                len(unk), unk.accession.nunique())

    lut: dict[str, tuple[str, str]] = {}
    for sid, acc, tis in zip(known.sample_id, known.accession, known.tissue):
        for t in _tokens(sid):
            lut.setdefault(t, (str(tis), str(acc)))

    rows = []
    for idx, r in unk.iterrows():
        for t in _tokens(r.sample_id):
            if t in lut:
                tis, src_acc = lut[t]
                rows.append({"index": idx, "accession": r.accession,
                             "sample_id": r.sample_id, "tissue": tis,
                             "rule": "token_join", "evidence": f"{t} also in {src_acc}",
                             "duplicate_suspect": src_acc})
                break
        else:
            if r.accession in TITLE_RULES:
                tis, why = TITLE_RULES[r.accession]
                rows.append({"index": idx, "accession": r.accession,
                             "sample_id": r.sample_id, "tissue": tis,
                             "rule": "title_rule", "evidence": why,
                             "duplicate_suspect": ""})
    d = pd.DataFrame(rows)
    if len(d):
        logger.info("resolvable: %d (%s)", len(d), d.rule.value_counts().to_dict())
        logger.info("tissues assigned: %s", d.tissue.value_counts().to_dict())
    return d


def apply(dry_run: bool = True, manifest_csv: str | None = None) -> pd.DataFrame:
    manifest_csv = manifest_csv or (RUNS / "sample_manifest.csv")
    d = plan(manifest_csv)
    if not len(d):
        return d
    if dry_run:
        logger.info("DRY RUN: would set tissue on %d sections", len(d))
        return d
    m = pd.read_csv(manifest_csv, low_memory=False)
    before = int(m.loc[d["index"], "tissue"].map(_is_unknown).sum())
    assert before == len(d), "a row to be filled already carries a tissue; refusing to overwrite"
    m.loc[d["index"], "tissue"] = d.set_index("index").tissue
    m.loc[d["index"], "tissue_source"] = d.set_index("index").rule
    m.to_csv(manifest_csv, index=False)
    d.to_csv(LOG, index=False)
    logger.info("filled %d tissues; %d sections still unknown", len(d),
                int(m[m.h5ad_path.notna()
                      & ~m.verdict.isin(["quarantined", "duplicate"])].tissue
                    .map(_is_unknown).sum()))
    return d
