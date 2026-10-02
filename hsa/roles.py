"""Scout (per database) and Curator (per dataset) agents."""
import re
import threading

from . import db, tools_api as T
from .agent import WEB_SEARCH, run_agent

KEEP_ROLES = {"matrix", "coords", "matrix+coords", "bundle"}
VERDICTS = ["ready", "bundle_only", "restricted", "not_human", "not_spatial", "no_processed_data"]

MISSION = """You work on HSA, the Human Spatial Atlas: a uniformly collected catalogue of every public HUMAN \
spatially-resolved transcriptomics dataset (Visium/Visium HD, Xenium, CosMx, MERFISH/MERSCOPE, Stereo-seq, \
Slide-seq/Curio, DBiT-seq, seqFISH, STARmap, GeoMx, HDST, legacy ST, etc.), modelled on scBaseCount.

The current phase needs only the processed CELL(or spot)-BY-GENE count matrix and the per-cell/spot \
SPATIAL COORDINATES for each sample. Images, raw reads (FASTQ/BAM), and per-transcript tables are out of scope \
for now, but record any links you find. Mouse and other species are excluded.

File roles you assign:
- matrix: count matrix (10x filtered_feature_bc_matrix.h5 or mtx+barcodes+features, Xenium cell_feature_matrix, \
CosMx exprMat_file, MERSCOPE cell_by_gene, csv/tsv count tables). Prefer filtered over raw 10x matrices.
- coords: coordinates (tissue_positions*.csv/parquet + scalefactors_json for Visium; Xenium cells.csv.gz/parquet; \
CosMx metadata_file; MERSCOPE cell_metadata; Slide-seq bead locations).
- matrix+coords: one file holding both (h5ad with obsm['spatial'], Seurat rds with images/coords, Stereo-seq gef).
- bundle: archive (tar/zip) that contains matrix/coords but is not available unpacked.
Be factual: never invent URLs; only use URLs returned by tools."""

# ---------------------------------------------------------------- tool specs
S_STR = {"type": "string"}


def spec(name, desc, props, req):
    return {"name": name, "description": desc,
            "input_schema": {"type": "object", "properties": props, "required": req}}


FILE_ITEM = {"type": "object", "properties": {
    "url": S_STR, "role": {"type": "string", "enum": sorted(KEEP_ROLES)}, "fmt": S_STR,
    "size_bytes": {"type": "integer"}}, "required": ["url", "role"]}

RECORD_PROPS = {
    "source": {"type": "string", "description": "database name, e.g. GEO, Zenodo, CELLxGENE, 10x Datasets, HTAN"},
    "accession": {"type": "string", "description": "stable id, e.g. GSE123456, zenodo:1234567"},
    "title": S_STR, "url": {"type": "string", "description": "landing page"},
    "publication": {"type": "string", "description": "DOI / PMID if known"},
    "technology": {"type": "string", "description": "e.g. Visium, Visium HD, Xenium, CosMx, MERSCOPE, Stereo-seq"},
    "tissue": S_STR, "disease": S_STR,
    "verdict": {"type": "string", "enum": VERDICTS},
    "notes": {"type": "string", "description": "anything a human should know (multi-species, mixed assays, caveats)"},
    "links": {"type": "array", "items": S_STR, "description": "other useful links (raw data, images, code, paper)"},
    "samples": {"type": "array", "description": "explicit per-sample file lists", "items": {
        "type": "object", "properties": {"sample_id": S_STR, "files": {"type": "array", "items": FILE_ITEM}},
        "required": ["sample_id", "files"]}},
    "selections": {"type": "array", "description": (
        "Alternative to `samples` for large datasets: regexes applied (re.search, case-insensitive) to file names "
        "already listed in THIS session. Each matching file is added with the given role; sample_id is the GSM "
        "accession for GEO files, else group(1) of sample_regex, else the accession."), "items": {
        "type": "object", "properties": {"pattern": S_STR, "role": {"type": "string", "enum": sorted(KEEP_ROLES)},
                                         "sample_regex": S_STR}, "required": ["pattern", "role"]}},
}


class Session:
    """Per-agent state: every file URL a listing tool has shown, so record_dataset can validate."""

    def __init__(self, curator=None):
        self.seen = {}           # url -> file dict
        self.recorded = []
        self.curator = curator

    def _remember(self, files):
        for f in files:
            if f.get("url") and f.get("role") != "dir":
                self.seen[f["url"]] = f
        return files

    # --- listing tools
    def geo_search(self, term, retmax=200, retstart=0):
        return T.geo_search(term, min(int(retmax), 500), int(retstart))

    def geo_series(self, gse):
        """Lists EVERY GSM directory up front (so `selections` cover all samples) and returns a compact summary."""
        m = T.geo_series(gse)
        self._remember(m["series_suppl"])
        allf = T.geo_all_sample_files([s["gsm"] for s in m.get("samples", [])])
        pats = {}
        for g, files in allf.items():
            self._remember(files)
            for f in files:
                if f["role"] == "image":
                    continue
                key = re.sub(r"^GSM\d+_?", "", f["name"])
                key = re.sub(r"\d+", "#", key)
                p = pats.setdefault(key, {"pattern": key, "role": f["role"], "n_samples": 0, "example": f["name"]})
                p["n_samples"] += 1
        m["sample_file_patterns"] = sorted(pats.values(), key=lambda x: -x["n_samples"])[:80]
        m["example_sample_files"] = {g: [f for f in v if f["role"] != "image"] for g, v in list(allf.items())[:2]}
        m["n_gsm_listed"] = len(allf)
        return m

    def zip_list(self, url):
        r = T.zip_list(url)
        self._remember(r["files"])
        r["files"] = [f for f in r["files"] if f["role"] not in ("image", "other")]
        return r

    def geo_sample_files(self, gsms):
        out = T.geo_sample_files(gsms)
        for v in out.values():
            self._remember(v)
        # hide images to save tokens; keep a count
        return {g: {"files": [f for f in v if f["role"] not in ("image",)],
                    "n_images_hidden": sum(f["role"] == "image" for f in v)} for g, v in out.items()}

    def zenodo_search(self, query, page=1):
        return T.zenodo_search(query, int(page))

    def zenodo_files(self, record_id):
        r = T.zenodo_files(record_id)
        self._remember(r["files"])
        return r

    def mendeley_files(self, dataset_id, version=None):
        r = T.mendeley_files(dataset_id, version)
        self._remember(r["files"])
        r["files"] = [f for f in r["files"] if f["role"] != "image"]
        return r

    def figshare_search(self, query, page=1):
        return T.figshare_search(query, int(page))

    def figshare_files(self, article_id):
        r = T.figshare_files(article_id)
        self._remember(r["files"])
        return r

    def list_directory(self, url):
        r = T.list_directory(url)
        self._remember(r.get("files", []))
        return r

    def http_request(self, url, method="GET", json_body=None, max_chars=20000):
        r = T.http_request(url, method, json_body, min(int(max_chars), 40000))
        # remember file-looking links found in pages so the agent may record them
        for u in re.findall(r"https?://\S+?\.(?:h5|h5ad|csv|csv\.gz|tsv\.gz|mtx\.gz|parquet|rds|gef|tar\.gz|tar|zip|json)\b",
                            r.get("text", "")):
            self.seen.setdefault(u, {"url": u, "name": u.rsplit("/", 1)[-1], "role": T.classify(u)})
        if method.upper() == "HEAD" and r.get("status") == 200:
            self.seen.setdefault(url, {"url": url, "name": url.rsplit("/", 1)[-1], "role": T.classify(url),
                                       "size_bytes": int(r["content_length"]) if r.get("content_length") else None})
        return r

    def wayback_fetch(self, url, max_chars=30000):
        """Archived copy of a page (for sites behind bot checkpoints, e.g. 10xgenomics.com)."""
        av = T.get(f"https://archive.org/wayback/available?url={url}").json()
        snap = av.get("archived_snapshots", {}).get("closest", {}).get("url")
        if not snap:
            return {"error": "no Wayback snapshot"}
        r = T.http_request(snap, "GET", None, min(int(max_chars), 60000))
        r["text"] = re.sub(r"https?://web\.archive\.org/web/\d+[a-z_]*/", "", r.get("text", ""))
        r["snapshot"] = snap
        for u in set(re.findall(r"https?://[^\s\"'<>|]+", r["text"])):
            if T.classify(u) in KEEP_ROLES or T.classify(u) == "matrix_raw":
                self.seen.setdefault(u, {"url": u, "name": u.rsplit("/", 1)[-1], "role": T.classify(u)})
        return r

    # --- write tools
    def note_source(self, name, url, access_method, est_human_datasets=0, notes=""):
        db.upsert_source(name, url, access_method, int(est_human_datasets or 0), notes)
        return {"ok": True}

    def queue_candidates(self, source, items, reason=""):
        n = db.add_candidates(source, items, reason)
        return {"queued_new": n, "submitted": len(items)}

    def record_dataset(self, **d):
        if d["verdict"] not in VERDICTS:
            raise ValueError(f"verdict must be one of {VERDICTS}")
        samples = {s["sample_id"]: list(s["files"]) for s in d.pop("samples", None) or []}
        for sel in d.pop("selections", None) or []:
            pat = re.compile(sel["pattern"], re.I)
            hit = 0
            for url, f in self.seen.items():
                name = f.get("name") or url.rsplit("/", 1)[-1]
                if not pat.search(name):
                    continue
                gsm = re.search(r"/(GSM\d+)/suppl/", url) or re.search(r"(GSM\d+)_", name)
                sm = re.search(sel["sample_regex"], name) if sel.get("sample_regex") else None
                sid = gsm.group(1) if gsm else (sm.group(1) if sm else d["accession"])
                samples.setdefault(sid, []).append({"url": url, "role": sel["role"],
                                                    "size_bytes": f.get("size_bytes")})
                hit += 1
            if hit == 0:
                raise ValueError(f"selection pattern {sel['pattern']!r} matched no listed file")
        # deterministic gate: unseen URLs, images, transcript tables are rejected
        problems = []
        for sid, files in samples.items():
            for f in files:
                cls = T.classify(f["url"])
                if f["url"] not in self.seen and not d["source"].lower().startswith("cellxgene"):
                    problems.append(f"{f['url']} was never returned by a listing tool (HEAD it or list its directory)")
                if cls in ("image", "transcripts", "segmentation") and f["role"] != "bundle":
                    problems.append(f"{f['url']} looks like {cls}, not {f['role']}")
                f.setdefault("size_bytes", self.seen.get(f["url"], {}).get("size_bytes"))
                f.setdefault("fmt", f["url"].rsplit("/", 1)[-1].split(".", 1)[-1][:20])
        if d["verdict"] == "ready":
            roles = {f["role"] for fs in samples.values() for f in fs}
            if not ("matrix+coords" in roles or {"matrix", "coords"} <= roles):
                problems.append(f"verdict 'ready' needs matrix AND coords (or matrix+coords); got roles {roles}")
        if problems:
            raise ValueError("record rejected:\n- " + "\n- ".join(problems[:15]))
        d["samples"] = [{"sample_id": k, "files": v} for k, v in samples.items()]
        d["n_samples"] = len(samples)
        n = db.record_dataset(d, self.curator or "scout")
        self.recorded.append(d["accession"])
        return {"ok": True, "n_files": n, "n_samples": len(samples)}


TOOL_SPECS = {
    "geo_search": spec("geo_search", "Search GEO series (Entrez gds syntax). Returns total + accession/title/n_samples/taxon/suppfile types.",
                       {"term": S_STR, "retmax": {"type": "integer"}, "retstart": {"type": "integer"}}, ["term"]),
    "geo_series": spec("geo_series", "GEO series metadata, GSM list, series-level files, AND a summary of per-sample file-name patterns across ALL GSMs (digits -> #) with two full example GSM listings. Usually enough to write `selections` without geo_sample_files.",
                       {"gse": S_STR}, ["gse"]),
    "geo_sample_files": spec("geo_sample_files", "List supplementary files for up to 60 GSMs (images hidden).",
                             {"gsms": {"type": "array", "items": S_STR}}, ["gsms"]),
    "zenodo_search": spec("zenodo_search", "Search Zenodo records (Elasticsearch query string).",
                          {"query": S_STR, "page": {"type": "integer"}}, ["query"]),
    "zenodo_files": spec("zenodo_files", "Zenodo record metadata + file list.", {"record_id": S_STR}, ["record_id"]),
    "figshare_search": spec("figshare_search", "Search figshare articles.", {"query": S_STR, "page": {"type": "integer"}}, ["query"]),
    "figshare_files": spec("figshare_files", "figshare article metadata + file list.", {"article_id": S_STR}, ["article_id"]),
    "zip_list": spec("zip_list", "List the members of a remote .zip WITHOUT downloading it (HTTP range reads). Member URLs "
                     "('<zip>#zipmember=<path>') can be recorded as matrix/coords; the downloader extracts only those members. "
                     "Use this for big bundles such as Xenium/CosMx outs.zip.", {"url": S_STR}, ["url"]),
    "wayback_fetch": spec("wayback_fetch", "Fetch the Internet Archive snapshot of a page that is blocked live (429/403 bot checkpoints, "
                          "e.g. www.10xgenomics.com/datasets/...). File links in it become recordable (HEAD them to confirm they are live).",
                          {"url": S_STR}, ["url"]),
    "mendeley_files": spec("mendeley_files", "Mendeley Data dataset metadata + ALL files (folders recursed) with direct download URLs. "
                           "Use for data.mendeley.com/datasets/<id>/<version> links.",
                           {"dataset_id": S_STR, "version": {"type": "integer"}}, ["dataset_id"]),
    "list_directory": spec("list_directory", "List an HTTP directory index (files, sizes, heuristic role).", {"url": S_STR}, ["url"]),
    "http_request": spec("http_request", "GET/POST/HEAD any URL or REST API. HTML is reduced to text + links; binaries are not downloaded (use HEAD for size).",
                         {"url": S_STR, "method": {"type": "string", "enum": ["GET", "POST", "HEAD"]},
                          "json_body": {"type": "object"}, "max_chars": {"type": "integer"}}, ["url"]),
    "note_source": spec("note_source", "Record/overwrite a database in the HSA source registry (keeps the links).",
                        {"name": S_STR, "url": S_STR, "access_method": {"type": "string", "description": "API / bulk download / manual / controlled-access, with the endpoint"},
                         "est_human_datasets": {"type": "integer"}, "notes": S_STR}, ["name", "url", "access_method"]),
    "queue_candidates": spec("queue_candidates", "Queue dataset accessions for per-dataset curation (used for large sources).",
                             {"source": S_STR, "items": {"type": "array", "items": {"type": "object", "properties": {
                                 "accession": S_STR, "title": S_STR, "url": S_STR}, "required": ["accession"]}},
                              "reason": S_STR}, ["source", "items"]),
    "record_dataset": spec("record_dataset", "Write the curated dataset with its matrix/coords files. Validated; on error fix and retry.",
                           RECORD_PROPS, ["source", "accession", "title", "verdict"]),
}


def _bind(session, names):
    return [TOOL_SPECS[n] for n in names], {n: getattr(session, n) for n in names}


# ---------------------------------------------------------------- scout
LANES = {
    "GEO": "NCBI GEO. Do NOT curate series yourself: run many geo_search queries (technology names, 'spatial transcriptomics', "
           "'10x Visium', 'Xenium', 'CosMx', 'MERFISH', 'Stereo-seq', 'Slide-seq', 'GeoMx', 'spatially resolved' ...; always "
           "AND \"Homo sapiens\"[Organism]), page through ALL results with retstart, and queue_candidates every human-plausible series.",
    "CELLxGENE": "CZ CELLxGENE Discover - already harvested deterministically; just note_source with its API.",
    "10x Genomics Datasets": "https://www.10xgenomics.com/datasets - public human Visium/Visium HD/Xenium demo datasets with direct "
                             "download links on cf.10xgenomics.com. Find each human spatial dataset, its output files, and record_dataset.",
    "Zenodo": "Zenodo - zenodo_search for human spatial transcriptomics deposits (Visium, Xenium, CosMx, MERFISH, Stereo-seq, ...). "
              "Page through results; queue_candidates plausible records (accession 'zenodo:<id>').",
    "figshare": "figshare - search for human spatial transcriptomics deposits; queue_candidates plausible ones ('figshare:<id>').",
    "HTAN": "Human Tumor Atlas Network (humantumoratlas.org, data via Synapse/CDS/ISB-CGC). Determine which spatial level-3/4 "
            "matrices are open-access and how to fetch them; record datasets if openly downloadable, else note_source as controlled/registration.",
    "HuBMAP": "HuBMAP (portal.hubmapconsortium.org; search API https://search.api.hubmapconsortium.org/v3/search). Find human spatial "
              "transcriptomics datasets (Visium, Xenium, CosMx, Slide-seq, MERFISH) with processed outputs and public download URLs (assets.hubmapconsortium.org).",
    "Vendor showcases": "Vendor public datasets: Bruker/NanoString CosMx (nanostring.com/products/cosmx .../ffpe-dataset), Vizgen MERSCOPE "
                        "showcase (vizgen.com/data-release-program), STOmics/BGI, Curio. Record openly downloadable human datasets.",
    "Spatial databases": "Aggregator databases: SODB, STOmicsDB, SpatialDB, SOAR, Aquila, SPASCER, CROST, SpatialOmicsDB, SpatialTME, "
                         "Spatial Omics DataBase, 'Museum of spatial transcriptomics'. For each: note_source with access method. Where they "
                         "redistribute processed human data with direct links, record datasets (or queue the underlying GEO accessions).",
    "EBI BioStudies/ArrayExpress": "EBI BioStudies / ArrayExpress (www.ebi.ac.uk/biostudies/api/v1/search?query=...). Find human spatial "
                                   "transcriptomics studies; queue_candidates their accessions (E-MTAB-...).",
    "HCA / Broad SCP / other": "Human Cell Atlas data portal (Azul API service.azul.data.humancellatlas.org), Broad Single Cell Portal "
                               "(spatial studies), Allen Brain Cell atlas (human MERFISH), Wellcome Sanger spatial portals, "
                               "Human Protein Atlas, CZI. note_source each; record openly downloadable human spatial datasets.",
}


def scout(lane, budget_usd):
    s = Session()
    names = list(TOOL_SPECS)
    tools, funcs = _bind(s, names)
    user = (f"Scout lane: **{lane}**.\n{LANES[lane]}\n\nFirst note_source for every database you touch (URL, access method, "
            f"estimated number of human spatial datasets). Be exhaustive within this lane; finish with a short summary of what you "
            f"found, what you could not access, and why.")
    final, usage = run_agent(f"scout_{re.sub(r'[^A-Za-z0-9]+', '_', lane)}", MISSION, user, tools, funcs,
                             effort="high", max_turns=150, budget_usd=budget_usd, server_tools=[WEB_SEARCH])
    return final, usage, s.recorded


# ---------------------------------------------------------------- curator
CURATE_TOOLS = ["geo_series", "geo_sample_files", "zenodo_files", "figshare_files", "zip_list", "wayback_fetch", "mendeley_files", "list_directory",
                "http_request", "record_dataset"]


def curate(cand, budget_usd):
    s = Session(curator="curator")
    tools, funcs = _bind(s, CURATE_TOOLS)
    user = (f"Curate candidate from {cand['source']}: accession={cand['accession']} title={cand.get('title','')!r} "
            f"url={cand.get('url','')}\n"
            + (f"Hints: {cand['reason']}\n" if cand.get("reason") else "") + "\n"
            "Tips: 10x Genomics datasets host per-file outputs next to the bundle (e.g. .../<name>/<name>_cell_feature_matrix.h5, "
            "_cells.csv.gz, _filtered_feature_bc_matrix.h5, _spatial.tar.gz). The live 10x site blocks bots: use wayback_fetch on the dataset page to get exact file URLs, then HEAD them; if only an outs.zip exists, zip_list it "
            "and record its members. Google Drive / login-walled data: verdict=restricted with the link.\n\n"
            "1. Inspect metadata and files (for GEO: geo_series, then geo_sample_files for the GSMs).\n"
            "2. Decide: human? spatial transcriptomics? which technology/tissue/disease? (mixed human+mouse series: keep only human samples).\n"
            "3. Call record_dataset ONCE with the verdict and, for every human spatial sample, its matrix + coords files "
            "(use `selections` regexes for many samples). If matrix/coords exist only inside an archive, verdict=bundle_only "
            "with the archive as role=bundle. If only raw reads or only non-spatial data: the appropriate verdict and no files.\n"
            "Then reply in English with one line: verdict + n samples.")
    final, usage = run_agent(f"curate_{cand['source']}_{cand['accession']}".replace(":", "_").replace("/", "_"),
                             MISSION, user, tools, funcs, effort="medium", max_turns=25, budget_usd=budget_usd)
    db.execute("UPDATE candidates SET status=? WHERE source=? AND accession=?",
               ("curated" if s.recorded else "failed", cand["source"], cand["accession"]))
    return final, usage
