"""Per-sample metadata harmonisation (schema follows HST-Corpus / CELLxGENE conventions).

  raw text (GEO characteristics, dataset pages)  --Sonnet-->  labels  --OLS (deterministic)-->  ontology IDs
  CELLxGENE: ontology-coded already (no LLM).   xlsx (HST-Corpus): override + validation of the LLM.

  python -m hsa.metadata pilot N   |  run  |  validate
"""
import json
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import quote

import anthropic

from . import db, tools_api as T
from .config import MODEL, PRICE, load_key

load_key()
CLIENT = anthropic.Anthropic(max_retries=6, timeout=900)
OLS = "https://www.ebi.ac.uk/ols4/api"

db.CON.executescript("""
CREATE TABLE IF NOT EXISTS sample_meta (
  source TEXT, accession TEXT, sample_id TEXT,
  species TEXT, organism_ontology_term_id TEXT, assay TEXT, assay_ontology_term_id TEXT,
  tissue TEXT, tissue_ontology_term_id TEXT, tissue_match TEXT,
  sex TEXT, sex_ontology_term_id TEXT,
  age_num REAL, age_unit TEXT, age_group TEXT, development_stage TEXT, development_stage_ontology_term_id TEXT,
  disease_state TEXT, disease_level1 TEXT, disease_level2 TEXT, disease_level3 TEXT,
  disease_ontology_term_id_level1 TEXT, disease_ontology_term_id_level2 TEXT, disease_ontology_term_id_level3 TEXT,
  disease_match TEXT, condition TEXT,
  self_reported_ethnicity TEXT, self_reported_ethnicity_ontology_term_id TEXT, donor_id TEXT,
  metadata_source TEXT, confidence TEXT, evidence TEXT, llm_json TEXT, updated REAL,
  PRIMARY KEY (source, accession, sample_id));
CREATE TABLE IF NOT EXISTS ontology_cache (ontology TEXT, query TEXT, obo_id TEXT, label TEXT, match TEXT,
  PRIMARY KEY (ontology, query));
CREATE TABLE IF NOT EXISTS ancestor_cache (obo_id TEXT PRIMARY KEY, ancestors TEXT);
""")

# ---------------------------------------------------------------- fixed vocabularies
ORGANISM = {"homo_sapiens": "NCBITaxon:9606", "mus_musculus": "NCBITaxon:10090"}
SEX = {"female": "PATO:0000383", "male": "PATO:0000384"}
ASSAY = {"Visium": ("visium", "EFO:0010961"), "Visium HD": ("visium_hd", "EFO:0920058"),
         "Xenium": ("xenium", "EFO:0022615"), "CosMx": ("cosmx", "EFO:0022994"), "MERFISH": ("merfish", "EFO:0008992"),
         "Stereo-seq": ("stereo_seq", "EFO:0920125"), "Slide-seq": ("slideseqv2", "EFO:0030062"),
         "GeoMx": ("geomx", "EFO:0920109"), "STARmap": ("starmap", "unknown"), "seqFISH": ("seqfish", "unknown"),
         "DBiT-seq": ("dbit_seq", "unknown"), "ISS/Cartana": ("cartana", "unknown"), "HDST": ("hdst", "unknown"),
         "ST (legacy)": ("spatial_transcriptomics", "EFO:0030005")}
HEALTHY = ("healthy", "PATO:0000461")
# level-1 categories (priority order) and level-2 categories; resolved to MONDO IDs via OLS at first use
LEVEL1 = ["cancer", "benign neoplasm", "nervous system disorder", "cardiovascular disorder", "respiratory system disorder",
          "digestive system disorder", "immune system disorder", "urinary system disorder", "reproductive system disorder",
          "integumentary system disorder", "musculoskeletal system disorder", "endocrine system disorder",
          "hematologic disease", "infectious disease", "metabolic disease", "injury"]
LEVEL2 = ["breast cancer", "lung cancer", "colorectal cancer", "liver cancer", "pancreatic neoplasm", "prostate cancer",
          "kidney cancer", "skin cancer", "central nervous system cancer", "malignant glioma", "female reproductive organ cancer",
          "male reproductive organ cancer", "stomach cancer", "esophageal cancer", "head and neck cancer", "bladder cancer",
          "thyroid cancer", "lymphoma", "leukemia", "sarcoma", "bone cancer", "neuroblastoma", "biliary tract cancer",
          "neurodegenerative disease", "autoimmune disorder of the nervous system", "epilepsy", "psychiatric disorder",
          "autoimmune disease", "inflammatory bowel disease", "chronic kidney disease", "diabetes mellitus",
          "cardiomyopathy", "coronary artery disorder", "interstitial lung disease", "chronic obstructive pulmonary disease",
          "asthma", "viral infectious disease", "bacterial infectious disease", "parasitic infection", "psoriasis",
          "dermatitis", "arthritis", "fibrosis", "liver disease", "lymphoid system disorder", "gastroenteritis",
          "pregnancy disorder", "congenital abnormality", "obesity disorder"]


def snake(s):
    return re.sub(r"[^a-z0-9]+", "_", (s or "").lower()).strip("_")


# ---------------------------------------------------------------- OLS resolver
def ols(ontology, query):
    """label/synonym -> (obo_id, label, match) with match in exact|fuzzy|none. Cached."""
    q = re.sub(r"[_\s]+", " ", (query or "").strip()).lower()
    if not q or q in ("unknown", "na", "n/a", "none", "not reported"):
        return ("unknown", "unknown", "none")
    hit = db.query("SELECT obo_id, label, match FROM ontology_cache WHERE ontology=? AND query=?", (ontology, q))
    if hit:
        return tuple(hit[0].values())
    res = ("unknown", q, "none")
    for exact in ("true", "false"):
        try:
            r = T.get(f"{OLS}/search", params={"q": q, "ontology": ontology, "exact": exact, "rows": 5,
                                               "queryFields": "label,synonym", "type": "class"})
            docs = [d for d in r.json()["response"]["docs"]
                    if d.get("obo_id", "").upper().startswith(ontology.upper()) and not d.get("is_obsolete")]
        except Exception:
            docs = []
        if docs:
            d = next((d for d in docs if d.get("label", "").lower() == q), docs[0])
            res = (d["obo_id"], d["label"], "exact" if exact == "true" else "fuzzy")
            break
    db.execute("INSERT OR REPLACE INTO ontology_cache VALUES (?,?,?,?,?)", (ontology, q, *res))
    return res


def ancestors(obo_id):
    hit = db.query("SELECT ancestors FROM ancestor_cache WHERE obo_id=?", (obo_id,))
    if hit:
        return set(json.loads(hit[0]["ancestors"]))
    pre = obo_id.split(":")[0].lower()
    iri = quote(quote(f"http://purl.obolibrary.org/obo/{obo_id.replace(':', '_')}", safe=""), safe="")
    out = set()
    try:
        js = T.get(f"{OLS}/ontologies/{pre}/terms/{iri}/hierarchicalAncestors?size=500").json()
        out = {t["obo_id"] for t in js.get("_embedded", {}).get("terms", []) if t.get("obo_id")}
    except Exception:
        pass
    db.execute("INSERT OR REPLACE INTO ancestor_cache VALUES (?,?)", (obo_id, json.dumps(sorted(out))))
    return out


_levels = {}
_lv_lock = threading.Lock()


def _level_ids(names):
    with _lv_lock:
        key = tuple(names)
        if key not in _levels:
            _levels[key] = [(n, ols("mondo", n)[0]) for n in names]
            _levels[key] = [(n, i) for n, i in _levels[key] if i != "unknown"]
        return _levels[key]


def disease_levels(mondo_id, label):
    """level3 = the specific term; level1/level2 = the highest-priority / most specific category it falls under."""
    anc = ancestors(mondo_id) | {mondo_id}
    l1 = next(((n, i) for n, i in _level_ids(LEVEL1) if i in anc), ("other_disease", "unknown"))
    cands = [(n, i) for n, i in _level_ids(LEVEL2) if i in anc]
    # keep the most specific level-2 category (drop any that is an ancestor of another candidate)
    cands = [c for c in cands if not any(c[1] in ancestors(o[1]) for o in cands if o != c)]
    l2 = cands[0] if cands else (snake(label), mondo_id)
    return (snake(l1[0]), l1[1]), (snake(l2[0]), l2[1]), (snake(label), mondo_id)


def age_fields(num, unit, prenatal):
    """-> age_group, development_stage label, HsapDv id"""
    if num is None or unit in (None, "unknown"):
        return "unknown", "unknown", "unknown"
    if prenatal:
        weeks = num if unit == "week" else num * 4.345 if unit == "month" else num / 7 if unit == "day" else None
        if weeks is None:
            return "fetal", "unknown", "unknown"
        w = int(round(weeks))
        lab = f"{w}{'st' if w % 10 == 1 and w != 11 else 'nd' if w % 10 == 2 and w != 12 else 'rd' if w % 10 == 3 and w != 13 else 'th'} week post-fertilization stage"
        i, l, _ = ols("hsapdv", lab)
        return ("embryonic" if w < 9 else "fetal"), l, i
    years = num if unit == "year" else num / 12 if unit == "month" else num / 52.18 if unit == "week" else num / 365.25
    grp = ("infant" if years < 2 else "child" if years < 13 else "adolescent" if years < 20 else
           "young_adult" if years < 40 else "middle_aged" if years < 60 else "late_adult")
    if unit == "month" and num < 24:
        i, l, _ = ols("hsapdv", f"{int(num)}-month-old stage")
    elif years >= 1:
        i, l, _ = ols("hsapdv", f"{int(years)}-year-old stage")
    else:
        i, l = "unknown", "unknown"
    return grp, l, i


# ---------------------------------------------------------------- raw metadata collection
def geo_soft(gse):
    r = T.get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&form=text&view=brief", timeout=120)
    out, cur = {}, None
    for line in r.text.splitlines():
        if line.startswith("^SAMPLE"):
            cur = line.split("=", 1)[1].strip()
            out[cur] = {"characteristics": []}
        elif cur and line.startswith("!Sample_"):
            k, _, v = line[8:].partition(" = ")
            if k.startswith("characteristics"):
                out[cur]["characteristics"].append(v)
            elif k in ("title", "source_name_ch1", "organism_ch1", "description", "treatment_protocol_ch1"):
                out[cur][k] = (out[cur].get(k, "") + " " + v).strip()[:400]
    return out


def samples_for(ds):
    rows = db.query("SELECT DISTINCT sample_id FROM files WHERE source=? AND accession=?", (ds["source"], ds["accession"]))
    return [r["sample_id"] for r in rows] or [ds["accession"]]


# ---------------------------------------------------------------- LLM extraction
SAMPLE_SCHEMA = {
    "type": "object", "additionalProperties": False,
    "required": ["sample_id", "species", "tissue", "sex", "age_num", "age_unit", "prenatal", "disease_state", "disease",
                 "condition", "ethnicity", "donor_id", "confidence", "evidence"],
    "properties": {
        "sample_id": {"type": "string"},
        "species": {"type": "string", "enum": ["homo_sapiens", "mus_musculus", "human_in_mouse_xenograft", "other", "unknown"]},
        "tissue": {"type": "string", "description": "most specific anatomical site stated, as an UBERON-style label (e.g. 'dorsolateral prefrontal cortex', 'colonic mucosa'); 'unknown' if not stated"},
        "sex": {"type": "string", "enum": ["female", "male", "unknown"]},
        "age_num": {"type": ["number", "null"], "description": "numeric age if stated (midpoint for ranges)"},
        "age_unit": {"type": "string", "enum": ["year", "month", "week", "day", "unknown"]},
        "prenatal": {"type": "boolean", "description": "true if embryonic/fetal (age is post-conception)"},
        "disease_state": {"type": "string", "enum": ["healthy", "disease", "unknown"],
                          "description": "healthy = normal / control / non-lesional / adjacent normal tissue"},
        "disease": {"type": "string", "description": "most specific diagnosis as a MONDO-style label (e.g. 'invasive ductal breast carcinoma', 'Alzheimer disease'); 'normal' if healthy; 'unknown'"},
        "condition": {"type": "string", "description": "other sample-level context: treatment, lesion vs non-lesional, region, timepoint, FFPE/fresh frozen; '' if none"},
        "ethnicity": {"type": "string", "description": "self-reported ethnicity / ancestry as stated (e.g. 'European', 'African American', 'Han Chinese'); 'unknown'"},
        "donor_id": {"type": "string", "description": "patient/donor identifier if stated, else ''"},
        "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
        "evidence": {"type": "string", "description": "the raw fields you used, briefly"}}}
SCHEMA = {"type": "object", "additionalProperties": False, "required": ["samples"],
          "properties": {"samples": {"type": "array", "items": SAMPLE_SCHEMA}}}

SYSTEM = """You harmonise sample metadata for HSA, a catalogue of human spatial transcriptomics data.
For EVERY listed sample return one record. Use only what the provided metadata states or clearly implies at the \
dataset level (e.g. a series on 'healthy adult human kidney' implies tissue=kidney, disease_state=healthy for all \
samples); never guess sex, age or ethnicity. Prefer the most specific tissue and disease that are stated. \
Tumour-adjacent normal tissue: disease_state=healthy, mention 'tumour-adjacent' in condition. \
Cell lines / organoids: put that in condition. Write labels the way UBERON/MONDO name them."""


def llm_extract(ctx, sample_ids, run):
    out = []
    for i in range(0, len(sample_ids), 80):            # chunk very large series
        chunk = sample_ids[i:i + 80]
        user = json.dumps({"dataset": ctx["dataset"], "samples": {s: ctx["samples"].get(s, {}) for s in chunk}},
                          default=str)[:150000]
        with CLIENT.messages.stream(
                model=MODEL, max_tokens=64000, system=SYSTEM,
                messages=[{"role": "user", "content": user}],
                output_config={"effort": "low", "format": {"type": "json_schema", "schema": SCHEMA}}) as st:
            msg = st.get_final_message()
        u = msg.usage
        step = {"in": u.input_tokens or 0, "out": u.output_tokens or 0, "cache_read": u.cache_read_input_tokens or 0,
                "cache_write": u.cache_creation_input_tokens or 0}
        db.log_usage(run, step, (step["in"] * PRICE["in"] + step["out"] * PRICE["out"]
                                 + step["cache_read"] * PRICE["cache_read"] + step["cache_write"] * PRICE["cache_write"]) / 1e6)
        if msg.stop_reason == "refusal":
            raise RuntimeError(f"refusal {msg.stop_details}")
        out += json.loads(next(b.text for b in msg.content if b.type == "text"))["samples"]
    return out


# ---------------------------------------------------------------- resolve + write
def resolve(ds, rec, source_tag):
    """LLM/xlsx label record -> full harmonised row (deterministic ontology mapping)."""
    sp = rec.get("species", "unknown")
    tech = ds.get("technology_norm") or ""
    assay, assay_id = ASSAY.get(tech, (snake(tech) or "unknown", "unknown"))
    t_id, t_lab, t_match = ols("uberon", rec.get("tissue"))
    sex = rec.get("sex", "unknown") if rec.get("sex") in SEX else "unknown"
    grp, dev, dev_id = age_fields(rec.get("age_num"), rec.get("age_unit"), rec.get("prenatal"))
    state = rec.get("disease_state", "unknown")
    if state == "healthy":
        l1 = l2 = l3 = HEALTHY
        dmatch = "fixed"
    elif rec.get("disease") and rec["disease"].lower() not in ("unknown", "normal"):
        d_id, d_lab, dmatch = ols("mondo", rec["disease"])
        if d_id == "unknown":
            l1 = l2 = l3 = ("unknown", "unknown")
        else:
            l1, l2, l3 = disease_levels(d_id, d_lab)
    else:
        l1 = l2 = l3 = ("unknown", "unknown")
        dmatch = "none"
    e_id, e_lab, _ = ols("hancestro", rec.get("ethnicity"))
    return {"source": ds["source"], "accession": ds["accession"], "sample_id": rec["sample_id"],
            "species": sp, "organism_ontology_term_id": ORGANISM.get(sp, "NCBITaxon:9606" if "xenograft" in sp else "unknown"),
            "assay": assay, "assay_ontology_term_id": assay_id,
            "tissue": snake(t_lab) if t_id != "unknown" else snake(rec.get("tissue")) or "unknown",
            "tissue_ontology_term_id": t_id, "tissue_match": t_match,
            "sex": sex, "sex_ontology_term_id": SEX.get(sex, "unknown"),
            "age_num": rec.get("age_num"), "age_unit": rec.get("age_unit") or "unknown", "age_group": grp,
            "development_stage": dev, "development_stage_ontology_term_id": dev_id,
            "disease_state": state, "disease_level1": l1[0], "disease_level2": l2[0], "disease_level3": l3[0],
            "disease_ontology_term_id_level1": l1[1], "disease_ontology_term_id_level2": l2[1],
            "disease_ontology_term_id_level3": l3[1], "disease_match": dmatch, "condition": rec.get("condition", ""),
            "self_reported_ethnicity": snake(e_lab) if e_id != "unknown" else "unknown",
            "self_reported_ethnicity_ontology_term_id": e_id, "donor_id": rec.get("donor_id", ""),
            "metadata_source": source_tag, "confidence": rec.get("confidence", ""), "evidence": rec.get("evidence", "")[:500],
            "llm_json": json.dumps(rec), "updated": time.time()}


def write(row):
    cols = list(row)
    db.execute(f"INSERT OR REPLACE INTO sample_meta ({','.join(cols)}) VALUES ({','.join('?' * len(cols))})",
               [row[c] for c in cols])


def _ds_rows(where="", args=()):
    from .export import norm_tech
    ds = db.query("SELECT * FROM datasets WHERE verdict IN ('ready','bundle_only') " + where, args)
    for d in ds:
        d["technology_norm"] = norm_tech(d["technology"])
    return ds


def harmonise_dataset(ds):
    sids = samples_for(ds)
    if ds["source"] == "CELLxGENE":
        return harmonise_cellxgene(ds, sids)
    ctx = {"dataset": {k: ds[k] for k in ("source", "accession", "title", "technology", "tissue", "disease", "notes")},
           "samples": {}}
    if ds["source"] == "GEO" and ds["accession"].startswith("GSE"):
        soft = geo_soft(ds["accession"])
        ctx["samples"] = {s: soft.get(s, {}) for s in sids}
        g = T._gds_summaries(T.get(f"{T.EU}/esearch.fcgi", params={"db": "gds", "retmode": "json",
                             "term": f"{ds['accession']}[ACCN] AND gse[Entry Type]"}).json()["esearchresult"]["idlist"][:1])
        if g:
            ctx["dataset"]["summary"] = g[0].get("summary", "")
    recs = llm_extract(ctx, sids, f"meta_{ds['source']}_{ds['accession']}"[:120])
    got = {r["sample_id"]: r for r in recs}
    n = 0
    for s in sids:
        r = got.get(s)
        if r is None:   # model skipped a sample: record as unknown rather than invent
            r = {"sample_id": s, "species": "unknown", "confidence": "low", "evidence": "not returned by model"}
        r["sample_id"] = s
        write(resolve(ds, r, "llm"))
        n += 1
    return n


_CXG = None


def harmonise_cellxgene(ds, sids):
    global _CXG
    if _CXG is None:
        _CXG = {d["dataset_id"]: d for d in T.get("https://api.cellxgene.cziscience.com/curation/v1/datasets",
                                                   timeout=180).json()}
    d = _CXG.get(ds["accession"], {})

    def one(field):
        v = d.get(field) or []
        return ("; ".join(x["label"] for x in v) or "unknown", "; ".join(x["ontology_term_id"] for x in v) or "unknown")
    tech = ds["technology_norm"]
    assay, assay_id = ASSAY.get(tech, (snake(tech), "unknown"))
    t_lab, t_id = one("tissue")
    s_lab, s_id = one("sex")
    dv_lab, dv_id = one("development_stage")
    e_lab, e_id = one("self_reported_ethnicity")
    dz = d.get("disease") or []
    if len(dz) == 1 and dz[0]["ontology_term_id"] == "PATO:0000461":
        l1 = l2 = l3 = HEALTHY
        state = "healthy"
    elif len(dz) == 1:
        l1, l2, l3 = disease_levels(dz[0]["ontology_term_id"], dz[0]["label"])
        state = "disease"
    else:
        lab, ids = one("disease")
        l1 = l2 = l3 = (snake(lab), ids)
        state = "mixed" if dz else "unknown"
    m = re.match(r"(\d+)-year-old", dv_lab)
    age = float(m.group(1)) if m else None
    grp = age_fields(age, "year", False)[0] if age is not None else ("fetal" if "post-fertilization" in dv_lab or "Carnegie" in dv_lab else "unknown")
    for s in sids:
        write({"source": ds["source"], "accession": ds["accession"], "sample_id": s, "species": "homo_sapiens",
               "organism_ontology_term_id": "NCBITaxon:9606", "assay": assay, "assay_ontology_term_id": assay_id,
               "tissue": snake(t_lab), "tissue_ontology_term_id": t_id, "tissue_match": "cellxgene",
               "sex": snake(s_lab), "sex_ontology_term_id": s_id, "age_num": age, "age_unit": "year" if age else "unknown",
               "age_group": grp, "development_stage": dv_lab, "development_stage_ontology_term_id": dv_id,
               "disease_state": state, "disease_level1": l1[0], "disease_level2": l2[0], "disease_level3": l3[0],
               "disease_ontology_term_id_level1": l1[1], "disease_ontology_term_id_level2": l2[1],
               "disease_ontology_term_id_level3": l3[1], "disease_match": "cellxgene", "condition": "",
               "self_reported_ethnicity": snake(e_lab), "self_reported_ethnicity_ontology_term_id": e_id,
               "donor_id": "", "metadata_source": "cellxgene", "confidence": "high", "evidence": "CELLxGENE schema",
               "llm_json": "", "updated": time.time()})
    return len(sids)


# ---------------------------------------------------------------- xlsx (HST-Corpus) override + validation
def xlsx_records():
    """xlsx rows -> {(source, accession): [label-record,...]} in the LLM record shape."""
    import pandas as pd
    from .harvest_xlsx import XLSX
    o = pd.read_excel(XLSX, sheet_name=0)
    rows = {r["row"]: r for r in db.query("SELECT * FROM xlsx_rows")}
    out = {}
    for i, r in o.iterrows():
        if i not in rows:
            continue
        x = rows[i]

        def val(c):
            v = r.get(c)
            return None if pd.isna(v) or str(v).lower() in ("na", "unknown", "nan") else v
        d3 = val("disease_level3")
        healthy = str(val("disease_level1") or "").startswith("healthy") and "disease" not in str(val("disease_level1"))
        unit = val("age_unit")
        rec = {"species": val("species") or "unknown", "tissue": (val("tissue") or "unknown").replace("_", " "),
               "sex": val("sex") or "unknown", "age_num": float(val("age_num")) if val("age_num") is not None and str(val("age_num")).replace('.', '', 1).isdigit() else None,
               "age_unit": unit if unit in ("year", "month", "week", "day") else "unknown",
               "prenatal": str(val("age_group") or "") in ("fetal", "embryonic") or unit == "carnegie",
               "disease_state": "healthy" if healthy else ("disease" if d3 else "unknown"),
               "disease": str(d3).replace("_", " ") if d3 and not healthy else "normal" if healthy else "unknown",
               "condition": "", "ethnicity": "unknown", "donor_id": "", "confidence": "high",
               "evidence": f"pretraining_data_overview.xlsx row {i} ({x['folder_name'] or x['sample']})",
               "_folder": x["folder_name"] or "", "_sample": x["sample"] or ""}
        out.setdefault((x["source"], x["accession"]), []).append(rec)
    return out


def _match_xlsx(sids, recs):
    """Pair xlsx rows to HSA sample ids: by name containment, else 1:1 when both sides have one entry."""
    pairs = {}
    for r in recs:
        for s in sids:
            if r["_folder"] and (r["_folder"].lower() in s.lower() or s.lower() in r["_folder"].lower()):
                pairs[s] = r
    if not pairs and len(sids) == 1 and len({(r['_folder'], r['_sample']) for r in recs}) == 1:
        pairs[sids[0]] = recs[0]
    return pairs


def apply_xlsx():
    """Override LLM rows with the user's curated xlsx values where a sample can be paired; returns validation pairs."""
    xr = xlsx_records()
    comp = []
    for (src, acc), recs in xr.items():
        ds = _ds_rows("AND source=? AND accession=?", (src, acc))
        if not ds:
            continue
        sids = samples_for(ds[0])
        for s, r in _match_xlsx(sids, recs).items():
            old = db.query("SELECT * FROM sample_meta WHERE source=? AND accession=? AND sample_id=?", (src, acc, s))
            new = resolve(ds[0], {**{k: v for k, v in r.items() if not k.startswith("_")}, "sample_id": s}, "xlsx")
            if old and old[0]["metadata_source"] == "llm":
                comp.append((old[0], new))
            if old:   # keep LLM-only fields the xlsx does not carry
                for k in ("condition", "donor_id", "self_reported_ethnicity", "self_reported_ethnicity_ontology_term_id"):
                    new[k] = old[0][k]
                new["llm_json"] = old[0]["llm_json"]
            write(new)
    return comp


def validate(comp):
    fields = ["species", "tissue_ontology_term_id", "sex", "age_group", "disease_state",
              "disease_level1", "disease_level2", "disease_ontology_term_id_level3"]
    rep = {}
    for f in fields:
        both = [(a[f], b[f]) for a, b in comp if b[f] not in ("unknown", None, "")]
        agree = sum(x == y for x, y in both)
        llm_unknown = sum(x in ("unknown", None, "") for x, _ in both)
        rep[f] = {"n_with_reference": len(both), "agree": agree, "llm_unknown": llm_unknown,
                  "pct_agree": round(100 * agree / len(both), 1) if both else None}
    return rep


# ---------------------------------------------------------------- driver
def run(limit=None, workers=8, budget=40.0, only_missing=True, xlsx_only=False):
    ds = _ds_rows()
    if xlsx_only:
        xs = {(r["source"], r["accession"]) for r in db.query("SELECT DISTINCT source, accession FROM xlsx_rows")}
        ds = [d for d in ds if (d["source"], d["accession"]) in xs and d["source"] != "CELLxGENE"]
    done = {(r["source"], r["accession"]) for r in db.query("SELECT DISTINCT source, accession FROM sample_meta")}
    todo = [d for d in ds if not only_missing or (d["source"], d["accession"]) not in done]
    # HST-Corpus datasets and non-GEO first, then GEO
    todo.sort(key=lambda d: (d["source"] == "GEO", d["accession"]))
    if limit:
        todo = todo[:limit]
    spent0 = db.query("SELECT COALESCE(SUM(usd),0) s FROM usage WHERE run LIKE 'meta_%'")[0]["s"]
    print(f"harmonising {len(todo)} datasets (metadata spend so far ${spent0:.2f}, cap ${budget})", flush=True)
    with ThreadPoolExecutor(workers) as ex:
        futs = {}
        for d in todo:
            futs[ex.submit(harmonise_dataset, d)] = d
        for i, fu in enumerate(as_completed(futs)):
            d = futs[fu]
            try:
                n = fu.result()
                msg = f"{n} samples"
            except Exception as e:
                msg = f"ERROR {type(e).__name__}: {str(e)[:200]}"
            spent = db.query("SELECT COALESCE(SUM(usd),0) s FROM usage WHERE run LIKE 'meta_%'")[0]["s"]
            print(f"[{i+1}/{len(todo)}] {d['source']} {d['accession'][:60]}: {msg} (${spent:.2f})", flush=True)
            if spent > budget:
                print("metadata budget reached; stopping", flush=True)
                ex.shutdown(cancel_futures=True)
                break
    comp = apply_xlsx()
    print("validation vs HST-Corpus xlsx (LLM before override):")
    print(json.dumps(validate(comp), indent=1))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "pilot":
        run(limit=int(sys.argv[2]) if len(sys.argv) > 2 else 20)
    elif cmd == "pilot_xlsx":
        run(limit=int(sys.argv[2]) if len(sys.argv) > 2 else 30, xlsx_only=True)
    elif cmd == "run":
        run()
    elif cmd == "validate":
        print(json.dumps(validate(apply_xlsx()), indent=1))
