"""Export the catalogue to CSV + a self-contained HTML dashboard.

  python -m hsa.export          -> exports/*.csv, index.html (GitHub Pages)
"""
import base64
import json
import re
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import db
from .config import DATA, ROOT, RUNS

OUT = ROOT / "exports"
PRIVATE = ROOT / "exports_private"     # gitignored: reference-list coverage, internal paths
DASH = ROOT                   # index.html at repo root is what GitHub Pages serves
TEMPLATE = Path(__file__).with_name("dashboard_template.html")

TECH = [("atera", "Atera"), ("visium hd", "Visium HD"), ("visium", "Visium"), ("xenium", "Xenium"), ("cosmx", "CosMx"),
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


# technologies that measure a spot or a bin, not a cell. Visium HD alone contributes ~260M 2um
# bins, so folding them into a "cells" total would overstate the atlas by a factor of three.
# Exact names, not substrings. "ST" as a substring matched STARmap, a single-cell imaging
# platform, and filed its cells as spots. Matching is exact against the resolved technology.
SPOT_TECH = ("Visium", "Visium HD", "Slide-seq", "Stereo-seq", "DBiT-seq", "GeoMx",
             "ST (legacy)", "ST")
# the resolver emits these for samples whose coordinates are one-per-cell rather than a grid
CELL_TECH = ("cell-level WTA", "cell-level (unknown platform)")


def headline():
    """One sentence: how many cells, across how many donors, in how many tissues.

    Counted only over human samples whose matrix actually loads, since a sample we cannot open
    contributes no cells. Donors are (accession, donor_id) pairs because donor ids are local to a
    study, and the count is a lower bound: not every sample states one.
    """
    qc = RUNS / "sample_qc.csv"
    if not qc.exists():
        return {}
    q = pd.read_csv(qc)
    try:
        sm = pd.DataFrame(db.query("SELECT source, accession, sample_id, species, donor_id,"
                                   " tissue_ontology_term_id FROM sample_meta"))
    except Exception:
        return {}
    if not len(sm):
        return {}
    m = q.merge(sm, on=["source", "accession", "sample_id"], how="left")
    # Prefer the technology resolved from file evidence over the depositor's label. Of 320 samples
    # labelled "Visium HD", only 205 sit on an HD bin grid: 64 are classic Visium and 49 are
    # cell-level with one coordinate per cell. Counting all 320 as HD both overstates that platform
    # and misfiles ~21M measured units on the wrong side of the cell/spot split.
    res_f = RUNS / "technology_resolved.csv"
    if res_f.exists():
        res = pd.read_csv(res_f)
        key = ["source", "accession", "sample_id"]
        for d in (m, res):
            d[key] = d[key].astype(str)
        m = m.merge(res[key + ["resolved"]], on=key, how="left")
        m["technology"] = m.resolved.where(m.resolved.notna() & ~m.resolved.eq("unknown"),
                                           m.technology)
    hum = m[m.species.isin(["homo_sapiens", "human_in_mouse_xenograft"])]
    # Drop sections that reached us twice. Zenodo mints a record id per version and authors mirror
    # GEO series elsewhere, so 92 rows are byte-identical copies of another row -- 29.2M units that
    # would otherwise be counted twice. See hsa/dedupe.py; copies are confirmed by hashing X, not
    # by matching shapes.
    from . import dedupe
    red = dedupe.redundant_keys()
    if red:
        k = list(zip(hum.source.astype(str), hum.accession.astype(str), hum.sample_id.astype(str)))
        hum = hum[[x not in red for x in k]]
    # A sample is in the atlas only if it is ON DISK. 103 rows pass QC with no .h5ad behind them
    # -- 59 quarantined, 37 that were never materialised, 7 bundle-only -- and counting them put
    # 47 samples in the headline that no one can open, which is also why the headline and the
    # load panel disagreed by exactly that number.
    mfp = RUNS / "sample_manifest.csv"
    if mfp.exists():
        mm = pd.read_csv(mfp, low_memory=False)
        kk = ["source", "accession", "sample_id"]
        mm[kk] = mm[kk].astype(str)
        on_disk = set(map(tuple, mm.loc[mm.h5ad_path.notna(), kk].values))
        # Quarantine has to be applied here too, not only via "has no file". Every quarantined
        # sample so far happened to lack an .h5ad, so the on-disk filter hid the gap; the first
        # quarantined sample that DID have a file (GSE300380/GSM9059356, coordinates assigned by
        # row order and demonstrably random) stayed in the headline while the load panel excluded
        # it.
        if "verdict" in mm.columns:
            on_disk -= set(map(tuple, mm.loc[mm.verdict.eq("quarantined"), kk].values))
        k = list(zip(hum.source.astype(str), hum.accession.astype(str), hum.sample_id.astype(str)))
        hum = hum[[x in on_disk for x in k]]
    ok = hum[hum.status.eq("ok") & hum.n_cells.notna()].copy()
    is_spot = ok.technology.astype(str).apply(
        lambda t: t not in CELL_TECH and t in SPOT_TECH)
    d = hum[hum.donor_id.notna() & ~hum.donor_id.isin(["unknown", ""])]
    t = hum[hum.tissue_ontology_term_id.str.startswith("UBERON", na=False)]
    per_tech = (ok.assign(unit=np.where(is_spot, "spot/bin", "cell"))
                  .groupby(["technology", "unit"]).n_cells.agg(["size", "sum"])
                  .reset_index().sort_values("sum", ascending=False))
    # Two different questions, and conflating them overstates the atlas. "cells"/"spots" below
    # count every unit materialised. But a Visium HD section that bin2cell segmented appears
    # BOTH as millions of bins and as its segmented cells, and the tissue is the same tissue --
    # so the honest answer to "how much of this atlas is single cells" subtracts those bins and
    # adds the segmented cells instead. 552M of the 598M spot total is Visium HD, and 67 of those
    # sections are segmented, so the distinction moves hundreds of millions of units.
    seg_cells, seg_n = 0, 0
    seg_dir = RUNS.parent / "visium_hd_cells"
    if seg_dir.exists():
        import h5py
        for f in seg_dir.glob("*.h5ad"):
            try:
                with h5py.File(f, "r") as h:
                    g = h["obs"]; idx = g[str(g.attrs.get("_index", "_index"))]
                    v = idx["values"] if isinstance(idx, h5py.Group) else idx
                    seg_cells += int(v.shape[0]); seg_n += 1
            except Exception:
                pass
    # which HD sections have a segmentation: the manifest knows, this table does not
    seg_keys = set()
    mf = RUNS / "sample_manifest.csv"
    if mf.exists():
        mm = pd.read_csv(mf, low_memory=False)
        if "hd_cells_path" in mm.columns:
            seg_keys = {(str(a), str(sid)) for a, sid in
                        zip(mm.loc[mm.hd_cells_path.notna(), "accession"],
                            mm.loc[mm.hd_cells_path.notna(), "sample_id"])}
    hd_seg_mask = pd.Series(
        [(str(a), str(sid)) in seg_keys for a, sid in zip(ok.accession, ok.sample_id)],
        index=ok.index)
    bins_unsegmented = int(ok.loc[is_spot & ~hd_seg_mask, "n_cells"].sum())

    # The headline pair is the segmentation-aware one: a Visium HD section that bin2cell
    # segmented contributes its CELLS, and its bins are not also counted, because they are the
    # same tissue measured twice. The unadjusted totals stay beside them so the adjustment is
    # visible rather than implied -- it moves 423M units between the two.
    cells_raw = int(ok.loc[~is_spot, "n_cells"].sum())
    spots_raw = int(ok.loc[is_spot, "n_cells"].sum())
    return {
        "cells": cells_raw + seg_cells,
        "spots": bins_unsegmented,
        "cells_raw": cells_raw,
        "spots_raw": spots_raw,
        "hd_segmented_samples": seg_n,
        "hd_segmented_cells": seg_cells,
        "donors": int(d.set_index(["accession", "donor_id"]).index.nunique()),
        "donor_frac": round(100 * len(d) / max(1, len(hum)), 1),
        "tissues": int(t.tissue_ontology_term_id.nunique()),
        "samples_loaded": int(len(ok)), "datasets": int(ok.accession.nunique()),
        "per_tech": [[r.technology, r.unit, int(r["size"]), int(r["sum"])]
                     for _, r in per_tech.iterrows()],
    }


def load_status():
    """Per technology: how many samples reached the atlas, and why the rest did not.

    The "in atlas" bucket is deliberately the SAME population the headline counts, so the two
    agree: human, deduplicated, and with an .h5ad actually on disk. Before that it was simply
    every row whose status began with "ok", which came to 7,820 against a headline of 7,662 --
    two numbers for one quantity, on one page, with nothing to explain the gap.

    The difference is itself worth showing rather than quietly subtracting, so samples that load
    but are kept out get an "excluded" bucket: quarantined, non-human, or a second copy of a
    section already counted. A reader can then see that the atlas is smaller than the number of
    files that opened, and why.
    """
    qc = RUNS / "sample_qc.csv"
    if not qc.exists():
        return []
    q = pd.read_csv(qc)
    key = ["source", "accession", "sample_id"]
    q[key] = q[key].astype(str)

    # the three things that keep a loadable sample out of the atlas
    from . import dedupe
    red = dedupe.redundant_keys()
    have_h5ad, quarantined = set(), set()
    mf = RUNS / "sample_manifest.csv"
    if mf.exists():
        mm = pd.read_csv(mf, low_memory=False)
        mm[key] = mm[key].astype(str)
        have_h5ad = set(map(tuple, mm.loc[mm.h5ad_path.notna(), key].values))
        if "verdict" in mm.columns:
            quarantined = set(map(tuple, mm.loc[mm.verdict.eq("quarantined"), key].values))
    try:
        sm = pd.DataFrame(db.query("SELECT source, accession, sample_id, species FROM sample_meta"))
        sm[key] = sm[key].astype(str)
        human = set(map(tuple, sm.loc[sm.species.isin(
            ["homo_sapiens", "human_in_mouse_xenograft"]), key].values))
    except Exception:
        human = None

    def bucket(r):
        x = str(r.status)
        if not x.startswith("ok"):
            if x.startswith("no_coords"):
                return "no coordinates"
            if x.startswith("no_recognised_matrix"):
                return "no matrix found"
            return "read error"
        k = (r.source, r.accession, r.sample_id)
        if k in quarantined or k in red or (human is not None and k not in human):
            return "excluded"
        if k not in have_h5ad:
            return "excluded"            # passes QC but was never written; not in the atlas
        return "in atlas"

    q["bucket"] = [bucket(r) for r in q.itertuples(index=False)]
    tab = q.pivot_table(index="technology", columns="bucket", aggfunc="size", fill_value=0)
    tab["n"] = tab.sum(1)
    # Every technology, not the top 12. There are 15, and the three smallest hold 9 samples --
    # truncating hid exactly those, so the bars summed to 2 fewer than the headline and the page
    # had a discrepancy with no cause a reader could find.
    tab = tab.sort_values("n", ascending=False)
    keys = ["in atlas", "no coordinates", "no matrix found", "read error", "excluded"]
    return [[t] + [int(tab.loc[t, k]) if k in tab.columns else 0 for k in keys]
            + [int(tab.loc[t, "n"])] for t in tab.index]


def hd_cells():
    """Cell-level Visium HD produced by bin2cell, which the QC table does not know about.

    These are a derived product: 2 um bins segmented into cells, so they are not rows in
    sample_qc.csv and would otherwise be invisible on the page. The label source matters and is
    reported alongside -- H&E segmentation sees nuclei directly, whereas the gene-expression route
    infers them from transcript density and recovers roughly half as many on a matched sample.
    """
    import h5py
    out = {"files": 0, "cells": 0, "by_method": {}}
    d = ROOT / "visium_hd_cells"
    if not d.exists():
        return out
    for f in sorted(d.glob("*.h5ad")):
        try:
            with h5py.File(f, "r") as h:
                cols = list(h["obs"].attrs.get("column-order", []))
                n = h["obs"][cols[0]].shape[0] if cols else 0
                u = h.get("uns/hsa")
                cc = u["cell_calling"][()].decode() if u and "cell_calling" in u else ""
        except Exception:
            continue
        k = "H&E + GEX salvage" if "versatile_he" in cc else (
            "gene expression only" if "fluo" in cc else "unknown")
        out["files"] += 1
        out["cells"] += int(n)
        m = out["by_method"].setdefault(k, {"files": 0, "cells": 0})
        m["files"] += 1
        m["cells"] += int(n)
    out["by_method"] = [[k, v["files"], v["cells"]] for k, v in
                        sorted(out["by_method"].items(), key=lambda kv: -kv[1]["cells"])]
    return out


def run():
    OUT.mkdir(exist_ok=True)
    # Rebuild the master sheet first: it joins the QC table, the biological labels, the panel and
    # grid derivations and the materialisation ledger, so it must not lag behind the page.
    try:
        from . import manifest
        manifest.build()
    except Exception as e:                       # a dashboard refresh should not die on this
        print(f"warning: master sheet not rebuilt ({type(e).__name__}: {e})")
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
                # numerator and denominator must span the SAME roles: counting every
                # downloaded file against only the matrix/coords ones reported 106%.
                "files_done": int((fs.role.isin(["matrix", "coords", "matrix+coords"])
                                   & fs.dl_status.eq("done")).sum()),
                "files_total": int(fs.role.isin(["matrix", "coords", "matrix+coords"]).sum()),
                "disk_gb": round(disk_gb, 1), "usd": usage,
                "cand_total": len(cand), "cand_done": int((cand.status != "pending").sum()) if len(cand) else 0},
        "datasets": ds[["source", "accession", "title", "technology_norm", "tissue_norm", "disease", "verdict",
                        "n_samples_with_files", "n_files", "n_files_downloaded", "n_bundles", "gb", "url"]]
                    .fillna("").to_dict("records"),
        "sources": src.fillna("").to_dict("records"),
        "samples": samples_payload,
        "headline": headline(),
        "hd_cells": hd_cells(),
        "load_status": load_status(),
    }
    from .logo import svg
    body = (TEMPLATE.read_text().replace("__HSA_LOGO__", svg(inline=True))
            .replace("__HSA_DATA__", json.dumps(payload, default=str).replace("</", "<\\/")))
    favicon = "data:image/svg+xml;base64," + base64.b64encode(svg().encode()).decode()
    full = ("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
            f"<link rel=\"icon\" href=\"{favicon}\">\n"
            "<style>*,*::before,*::after{box-sizing:border-box}:root{color-scheme:light;"
            "padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}"
            "body{margin:0;-webkit-font-smoothing:antialiased}img{max-width:100%}[hidden]{display:none!important}</style>\n"
            "</head>\n<body>\n" + body + "\n</body></html>\n")
    from .protect import protect
    locked = protect(full, favicon, svg(inline=False).replace('xmlns=', 'class="logo" role="img" aria-label="HSA logo" xmlns='))
    (DASH / "index.html").write_text(locked or full)   # encrypted page when .dashboard_password exists
    (DASH / ".nojekyll").touch()
    print(f"exported {len(ds)} datasets, {len(fs)} files -> {OUT}, {DASH/'index.html'}")


if __name__ == "__main__":
    run()
