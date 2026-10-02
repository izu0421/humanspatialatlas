"""Cross-check HSA against other spatial databases' sample lists (HEST-1k, DeepSpaceDB, SODB).

Each list is reduced to source accessions (GEO series, Zenodo, 10x dataset pages, ...). Accessions HSA has never
seen are queued for curation at top priority; the overlap is reported per reference database.

  python -m hsa.harvest_refs
"""
import io
import re
from collections import defaultdict

import pandas as pd

from . import db, tools_api as T

HEST_CSV = "https://raw.githubusercontent.com/mahmoodlab/HEST/main/assets/HEST_v1_1_0.csv"
DSDB = "https://genomics.virus.kyoto-u.ac.jp/deepspacedb/get_data"
SODB = "https://gene.ai.tencent.com/SpatialOmics/api/pysodb/info"


def norm_link(url):
    """URL -> (source, accession) in HSA's conventions."""
    u = str(url or "").strip()
    if m := re.search(r"(GSE\d+)", u):
        return "GEO", m.group(1)
    if m := re.search(r"zenodo\.org/(?:records?|record)/(\d+)|zenodo\.(\d+)", u):
        return "Zenodo", f"zenodo:{m.group(1) or m.group(2)}"
    if m := re.search(r"figshare\.com/.*?/(\d{6,})", u):
        return "figshare", f"figshare:{m.group(1)}"
    if "10xgenomics.com/datasets/" in u or "10xgenomics.com/resources/datasets/" in u:
        slug = re.sub(r"[?#].*", "", u).rstrip("/").rsplit("/", 1)[-1]
        return "10x Datasets", f"10xgenomics.com/datasets/{slug}"
    if m := re.search(r"(E-MTAB-\d+)", u):
        return "ArrayExpress", m.group(1)
    if m := re.search(r"(HTA\d+|HT\d{3}[A-Z0-9]+)", u):
        return "HTAN", m.group(1)
    if u.startswith("http"):
        return "other", re.sub(r"^https?://(www\.)?", "", u).rstrip("/")[:200]
    return None, None


def hest():
    h = pd.read_csv(io.StringIO(T.get(HEST_CSV).text))
    h = h[h.species == "Homo sapiens"]
    out = defaultdict(set)
    for _, r in h.iterrows():
        for col in ("download_page_link1", "study_link"):
            src, acc = norm_link(r.get(col))
            if src and src != "other":
                out[(src, acc)].add(r["id"])
                break
        else:
            src, acc = norm_link(r.get("download_page_link1"))
            if src:
                out[(src, acc)].add(r["id"])
    return out, len(h)


def deepspacedb():
    d = T.get(DSDB, timeout=180).json()["data"]
    out = defaultdict(set)
    n = 0
    for r in d:
        if r.get("organism") not in ("human", "Homo sapiens"):
            continue
        n += 1
        src, acc = ("GEO", r["series_id"]) if str(r.get("series_id") or "").startswith("GSE") else norm_link(r.get("URL"))
        if src:
            out[(src, acc)].add(r["sample_id"] or r["ID"])
    return out, n


def _gsm_to_gse(gsms):
    """Map GSM -> GSE via Entrez (batched)."""
    res = {}
    for i in range(0, len(gsms), 150):
        chunk = gsms[i:i + 150]
        r = T.get(f"{T.EU}/esearch.fcgi", params={"db": "gds", "retmode": "json", "retmax": 500,
                                                   "term": " OR ".join(f"{g}[ACCN]" for g in chunk)})
        ids = r.json()["esearchresult"].get("idlist", [])
        for j in range(0, len(ids), 150):
            js = T.get(f"{T.EU}/esummary.fcgi", params={"db": "gds", "retmode": "json",
                                                        "id": ",".join(ids[j:j + 150])}).json()["result"]
            for uid in js.get("uids", []):
                s = js[uid]
                if s.get("accession", "").startswith("GSM") and s.get("gse"):
                    res[s["accession"]] = "GSE" + str(s["gse"]).split(";")[0]
    return res


def sodb():
    rows = T.get(SODB, timeout=180).json()["data"]
    st = [r for r in rows if r[0] in ("Spatial Transcriptomics", "Spatial MultiOmics")]
    gsm_of = {}
    for cat, ds, exp in st:
        if m := re.search(r"(GSM\d+)", exp):
            gsm_of[(ds, exp)] = m.group(1)
    m = _gsm_to_gse(sorted(set(gsm_of.values())))
    out = defaultdict(set)
    for (ds, exp), g in gsm_of.items():
        if g in m:
            out[("GEO", m[g])].add(exp)
    no_acc = sorted({ds for cat, ds, exp in st if (ds, exp) not in gsm_of})
    return out, len(st), no_acc


def run():
    known = {(r["source"], r["accession"]): r["status"] for r in db.query("SELECT source, accession, status FROM candidates")}
    known.update({(r["source"], r["accession"]): "curated" for r in db.query("SELECT source, accession FROM datasets")})
    verdict = {(r["source"], r["accession"]): r["verdict"] for r in db.query("SELECT source, accession, verdict FROM datasets")}
    report = []
    refs = {"HEST-1k": hest, "DeepSpaceDB": deepspacedb, "SODB": sodb}
    for name, fn in refs.items():
        res = fn()
        accs, n_items = res[0], res[1]
        new = [k for k in accs if k not in known]
        ready = sum(verdict.get(k) == "ready" for k in accs)
        for src, acc in new:
            url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}" if src == "GEO" else \
                  ("https://" + acc if not acc.startswith(("zenodo:", "figshare:")) else "")
            db.execute("INSERT OR IGNORE INTO candidates (source, accession, title, url, reason, added) VALUES (?,?,?,?,?,?)",
                       (src, acc, "", url, f"ref:{name} ({len(accs[(src, acc)])} samples there)", 0.0))
        report.append({"reference": name, "items": n_items, "accessions": len(accs),
                       "already_in_HSA": len(accs) - len(new), "ready_in_HSA": ready, "new_queued": len(new),
                       "new_by_source": dict(pd.Series([s for s, _ in new]).value_counts()) if new else {}})
        if name == "SODB":
            report[-1]["sodb_datasets_without_GEO_accession"] = len(res[2])
        for k in accs:
            known.setdefault(k, "pending")
        print(report[-1], flush=True)
    return report


if __name__ == "__main__":
    run()
