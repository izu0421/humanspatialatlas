"""Queue every HUMAN entry of pretraining_data_overview.xlsx (HST-Corpus-112M) as curation candidates."""
import re

import pandas as pd

from . import db
from .config import ROOT

XLSX = ROOT.parent / "pretraining_data_overview.xlsx"

HOSTS = [("10xgenomics", "10x Datasets"), ("vizgen", "Vizgen"), ("nanostring", "NanoString"), ("cellxgene", "CELLxGENE"),
         ("zenodo", "Zenodo"), ("ncbi.nlm.nih.gov", "GEO"), ("drive.google", "Google Drive"),
         ("broadinstitute", "Broad SCP"), ("datadryad", "Dryad"), ("brainimagelibrary", "Brain Image Library")]


def _source(link):
    return next((name for k, name in HOSTS if k in link), "other")


def _accession(src, link):
    if src == "GEO" and (m := re.search(r"(GSE\d+)", link)):
        return m.group(1)
    if src == "Zenodo" and (m := re.search(r"(?:records?|zenodo\.)/?(\d+)", link)):
        return f"zenodo:{m.group(1)}"
    return re.sub(r"^https?://(www\.)?", "", link).rstrip("/")[:200]


def run():
    o = pd.read_excel(XLSX, sheet_name=0)
    h = o[o.species == "homo_sapiens"].copy()
    h["link"] = h.data_access_link.fillna(h.download_url).fillna(h.raw_on_farm.map(lambda x: f"offline:{x}"))
    h["link"] = h["link"].fillna(h.download_url).astype(str).str.split(";").str[0].str.strip()
    db.execute("""CREATE TABLE IF NOT EXISTS xlsx_rows (row INTEGER PRIMARY KEY, folder_name TEXT, sample TEXT, study TEXT, assay TEXT,
                  tissue TEXT, link TEXT, download_url TEXT, source TEXT, accession TEXT)""")
    n = 0
    for link, g in h.groupby("link"):
        src = _source(link)
        acc = _accession(src, link)
        for i, r in g.iterrows():
            db.execute("INSERT OR REPLACE INTO xlsx_rows VALUES (?,?,?,?,?,?,?,?,?,?)",
                       (int(i), None if pd.isna(r.folder_name) else r.folder_name, str(r["sample"]), str(r.study), r.assay, r.tissue, link,
                        None if pd.isna(r.download_url) else r.download_url, src, acc))
        if src == "CELLxGENE" or link.startswith("offline:"):
            continue   # covered by the deterministic CELLxGENE harvest; checked in coverage report
        dls = [u for u in g.download_url.dropna().unique()][:25]
        hint = (f"HST-Corpus-112M (must cover): {len(g)} human samples, assay={','.join(g.assay.unique())}, "
                f"tissue={','.join(map(str, g.tissue.unique()))}; sample folders={list(g.folder_name)[:25]}; "
                f"known bundle URLs={dls}")
        cur = db.execute("INSERT OR IGNORE INTO candidates (source, accession, title, url, reason, added) VALUES (?,?,?,?,?,?)",
                         (src, acc, str(g.study.iloc[0]), link, hint, 0.0))   # added=0 -> curated first
        if cur.rowcount == 0:   # already queued (e.g. GEO harvest): attach the hint
            db.execute("UPDATE candidates SET reason = reason || ' | ' || ? WHERE source=? AND accession=?", (hint, src, acc))
        n += 1
    print(f"xlsx: {len(h)} human rows, {h.link.nunique()} links, {n} candidates queued/updated")
    print(pd.Series([_source(l) for l in h.link.astype(str)]).value_counts().to_string())
