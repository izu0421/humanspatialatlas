"""Export the catalogue to CSV + a self-contained HTML dashboard.

  python -m hsa.export          -> exports/*.csv, index.html (GitHub Pages)
"""
import base64
import json
import re
import shutil
import time
from pathlib import Path

import pandas as pd

from . import db
from .config import DATA, ROOT

OUT = ROOT / "exports"
PRIVATE = ROOT / "exports_private"     # gitignored: reference-list coverage, internal paths
DASH = ROOT                   # index.html at repo root is what GitHub Pages serves
TEMPLATE = Path(__file__).with_name("dashboard_template.html")

TECH = [("visium hd", "Visium HD"), ("visium", "Visium"), ("xenium", "Xenium"), ("cosmx", "CosMx"),
        ("merscope", "MERFISH"), ("merfish", "MERFISH"), ("stereo", "Stereo-seq"), ("slide", "Slide-seq"),
        ("curio", "Slide-seq"), ("geomx", "GeoMx"), ("dbit", "DBiT-seq"), ("starmap", "STARmap"),
        ("seqfish", "seqFISH"), ("cartana", "ISS/Cartana"), ("in situ seq", "ISS/Cartana"), ("hdst", "HDST"),
        ("spatial transcriptomics", "ST (legacy)"), ("legacy st", "ST (legacy)")]


def norm_tech(t):
    t = (t or "").lower()
    return next((v for k, v in TECH if k in t), "Other" if t else "Unknown")


def norm_tissue(t):
    t = (t or "").split(";")[0].split(",")[0].strip().lower()
    return t[:40] or "unspecified"


def tables():
    ds = pd.DataFrame(db.query("SELECT * FROM datasets"))
    fs = pd.DataFrame(db.query("SELECT * FROM files"))
    src = pd.DataFrame(db.query("SELECT * FROM sources"))
    cand = pd.DataFrame(db.query("SELECT source, accession, status FROM candidates"))
    keep = fs[fs.role.isin(["matrix", "coords", "matrix+coords"])]
    agg = keep.groupby(["source", "accession"]).agg(
        n_files=("url", "size"), n_files_downloaded=("dl_status", lambda s: int((s == "done").sum())),
        gb=("size_bytes", lambda s: round(s.fillna(0).sum() / 1e9, 3)),
        n_samples_with_files=("sample_id", "nunique")).reset_index()
        # bundles are counted separately
    bun = fs[fs.role == "bundle"].groupby(["source", "accession"]).size().rename("n_bundles").reset_index()
    ds = ds.merge(agg, how="left", on=["source", "accession"]).merge(bun, how="left", on=["source", "accession"])
    for c in ["n_files", "n_files_downloaded", "n_samples_with_files", "n_bundles"]:
        ds[c] = ds[c].fillna(0).astype(int)
    ds["gb"] = ds["gb"].fillna(0.0)
    # an archive-only dataset whose matrix/coords members have been extracted is effectively ready
    ds["verdict_curator"] = ds.verdict
    ds.loc[(ds.verdict == "bundle_only") & (ds.n_files_downloaded > 0), "verdict"] = "ready"
    ds["technology_norm"] = ds.technology.map(norm_tech)
    ds["tissue_norm"] = ds.tissue.map(norm_tissue)
    ds["updated"] = pd.to_datetime(ds.updated, unit="s").dt.strftime("%Y-%m-%d %H:%M")
    ds = ds.sort_values(["verdict", "source", "accession"])
    return ds, fs, src, cand


def xlsx_coverage(ds):
    try:
        x = pd.DataFrame(db.query("SELECT * FROM xlsx_rows"))
    except Exception:
        return pd.DataFrame()
    v = ds.set_index(["source", "accession"]).verdict.to_dict()
    cxg_links = " ".join(ds[ds.source == "CELLxGENE"].links.astype(str))

    def status(r):
        if str(r.link).startswith("offline:"):
            return "offline (internal only)"
        if r.source == "CELLxGENE":
            m = re.search(r"collections/([0-9a-f\-]{36})", r.link)
            return "ready" if m and m.group(1) in cxg_links else "not matched"
        return v.get((r.source, r.accession), "pending curation")
    x["hsa_status"] = x.apply(status, axis=1)
    return x


def run():
    OUT.mkdir(exist_ok=True)
    DASH.mkdir(exist_ok=True)
    ds, fs, src, cand = tables()
    cols = ["source", "accession", "title", "technology_norm", "technology", "tissue", "disease", "verdict", "verdict_curator",
            "n_samples", "n_samples_with_files", "n_files", "n_files_downloaded", "n_bundles", "gb", "url",
            "publication", "notes", "links", "curated_by", "updated"]
    ds[cols].to_csv(OUT / "hsa_datasets.csv", index=False)
    fs.drop(columns=["local_path"]).to_csv(OUT / "hsa_files.csv", index=False)
    src.assign(updated=pd.to_datetime(src.updated, unit="s").dt.strftime("%Y-%m-%d")).to_csv(OUT / "hsa_sources.csv", index=False)
    xc = xlsx_coverage(ds)
    if len(xc):
        PRIVATE.mkdir(exist_ok=True)
        xc.to_csv(PRIVATE / "hst_corpus_coverage.csv", index=False)

    # per-sample harmonised metadata
    try:
        sm = pd.DataFrame(db.query("SELECT * FROM sample_meta"))
        conf = pd.DataFrame(db.query("SELECT * FROM xlsx_conflicts"))
    except Exception:
        sm, conf = pd.DataFrame(), pd.DataFrame()
    samples_payload = {}
    if len(sm):
        sm = sm.drop(columns=["llm_json"])
        sm["updated"] = pd.to_datetime(sm.updated, unit="s").dt.strftime("%Y-%m-%d %H:%M")
        sm.to_csv(OUT / "hsa_samples.csv", index=False)
        hum = sm[sm.species.isin(["homo_sapiens", "human_in_mouse_xenograft"])]
        known = lambda c: int((~hum[c].isin(["unknown", "", None]) & hum[c].notna()).sum())
        top = lambda c, n: [[k, int(v)] for k, v in hum.loc[~hum[c].isin(["unknown", ""]), c].value_counts().head(n).items()]
        samples_payload = {
            "n": len(hum), "n_nonhuman": int(len(sm) - len(hum)),
            "completeness": [[lab, known(c)] for lab, c in [("Tissue", "tissue_ontology_term_id"), ("Disease", "disease_state"),
                             ("Sex", "sex"), ("Age group", "age_group"), ("Development stage", "development_stage_ontology_term_id"),
                             ("Ethnicity", "self_reported_ethnicity_ontology_term_id"), ("Donor ID", "donor_id")]],
            "sex": [[k, int(v)] for k, v in hum.sex.value_counts().items()],
            "age_group": [[k, int(v)] for k, v in hum.age_group.value_counts().items()],
            "disease_state": [[k, int(v)] for k, v in hum.disease_state.value_counts().items()],
            "disease_level1": top("disease_level1", 12), "disease_level2": top("disease_level2", 14),
            "tissue": top("tissue", 16),
            "source": [[k, int(v)] for k, v in hum.metadata_source.value_counts().items()],
            "conflicts": len(conf)}
    if len(conf):
        PRIVATE.mkdir(exist_ok=True)
        conf.to_csv(PRIVATE / "xlsx_conflicts.csv", index=False)

    usage = db.query("SELECT ROUND(SUM(usd),2) usd FROM usage")[0]["usd"] or 0
    disk_gb = sum(p.stat().st_size for p in DATA.rglob("*") if p.is_file()) / 1e9 if DATA.exists() else 0
    payload = {
        "generated": time.strftime("%Y-%m-%d %H:%M"),
        "kpi": {"datasets": len(ds), "ready": int((ds.verdict == "ready").sum()),
                "samples": int(ds.loc[ds.verdict == "ready", "n_samples_with_files"].sum()),
                "files_done": int((fs.dl_status == "done").sum()),
                "files_total": int(fs.role.isin(["matrix", "coords", "matrix+coords"]).sum()),
                "disk_gb": round(disk_gb, 1), "usd": usage,
                "cand_total": len(cand), "cand_done": int((cand.status != "pending").sum()) if len(cand) else 0},
        "datasets": ds[["source", "accession", "title", "technology_norm", "tissue_norm", "disease", "verdict",
                        "n_samples_with_files", "n_files", "n_files_downloaded", "n_bundles", "gb", "url"]]
                    .fillna("").to_dict("records"),
        "sources": src.fillna("").to_dict("records"),
        "samples": samples_payload,
    }
    from .logo import svg
    body = (TEMPLATE.read_text().replace("__HSA_LOGO__", svg(inline=True))
            .replace("__HSA_DATA__", json.dumps(payload, default=str).replace("</", "<\\/")))
    favicon = "data:image/svg+xml;base64," + base64.b64encode(svg().encode()).decode()
    (DASH / "index.html").write_text(        # standalone page (GitHub Pages / local)
        "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
        f"<link rel=\"icon\" href=\"{favicon}\">\n"
        "<style>*,*::before,*::after{box-sizing:border-box}:root{color-scheme:light;"
        "padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}"
        "body{margin:0;-webkit-font-smoothing:antialiased}img{max-width:100%}[hidden]{display:none!important}</style>\n"
        "</head>\n<body>\n" + body + "\n</body></html>\n")
    (DASH / ".nojekyll").touch()
    print(f"exported {len(ds)} datasets, {len(fs)} files -> {OUT}, {DASH/'index.html'}")


if __name__ == "__main__":
    run()
