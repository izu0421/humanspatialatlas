"""Deterministic data-access helpers. Each returns plain Python objects (JSON-serialisable).

The agent sees these as tools; the harvester/downloader call them directly.
"""
import html
import json
import re
import threading
import time
from urllib.parse import urljoin

import requests

from .config import HTTP_UA

S = requests.Session()
S.headers["User-Agent"] = HTTP_UA

# ---------------------------------------------------------------- rate limiting
_rl_lock = threading.Lock()
_last = {}
MIN_GAP = {"ncbi.nlm.nih.gov": 0.4, "zenodo.org": 0.6, "figshare.com": 0.3}


def _throttle(url):
    host = next((h for h in MIN_GAP if h in url), None)
    if not host:
        return
    with _rl_lock:
        wait = _last.get(host, 0) + MIN_GAP[host] - time.time()
        if wait > 0:
            time.sleep(wait)
        _last[host] = time.time()


def get(url, **kw):
    for attempt in range(4):
        _throttle(url)
        try:
            r = S.get(url, timeout=kw.pop("timeout", 60), **kw)
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2 ** attempt * 3)
                continue
            return r
        except requests.RequestException:
            time.sleep(2 ** attempt * 3)
    return S.get(url, timeout=60, **kw)


# ---------------------------------------------------------------- filename classifier
# roles: matrix | coords | matrix+coords | transcripts | image | bundle | other
_RULES = [
    ("other", r"differential_expression|^analysis|/analysis/|cluster|umap|pca|projection"),
    ("segmentation", r"cell_boundar|nucleus_boundar|segmentation|masks?\."),
    ("image", r"\.(tiff?|png|jpe?g|svs|ndpi|btf|czi|qptiff|ome\.tif)(\.gz)?$|morphology|_he_image|hires_image|lowres_image"),
    ("transcripts", r"transcripts?\.(csv|parquet|zarr)|detected_transcripts|tx_file|_tx\.csv|molecules|spots\.csv"),
    ("matrix+coords", r"\.h5ad(\.gz)?$|\.gef$|\.cellbin\.gef$|_seurat.*\.rds$|\.rds$|\.h5seurat$|\.zarr\.zip$|\.loom$"),
    ("coords", r"_coords\.|spot_coord|tissue_positions|scalefactors_json|(^|[_\-.])cells\.(csv|parquet)|cell_metadata|metadata_file|"
               r"beadlocations|bead_locations|coordinates|positions|centroids|spatial\.tar\.gz|spatial\.zip|"
               r"_spatial\.csv|xy\.csv|locations?\.csv"),
    ("matrix", r"stdata|_st_data|filtered_feature_bc_matrix|cell_feature_matrix|feature_bc_matrix|matrix\.mtx|barcodes\.tsv|"
               r"features\.tsv|genes\.tsv|exprmat|cell_by_gene|counts?[._\-]|_count\.|expression|\.mtx(\.gz)?$|"
               r"\.h5$|dge\.|_umi"),
    ("bundle", r"\.(tar|tar\.gz|tgz|zip|7z|rar)$"),
]


def classify(name: str) -> str:
    n = name.lower().rsplit("/", 1)[-1]
    if n.startswith("._"):          # macOS resource-fork junk inside archives
        return "other"
    if "raw_feature_bc_matrix" in n:
        return "matrix_raw"      # unfiltered 10x - recorded but not preferred
    for role, pat in _RULES:
        if re.search(pat, n):
            return role
    return "other"


def _size_to_bytes(s):
    m = re.match(r"([\d.]+)\s*([KMGT]?)", s.strip(), re.I)
    if not m:
        return None
    return int(float(m.group(1)) * {"": 1, "K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}[m.group(2).upper()])


def list_directory(url: str):
    """Parse an Apache/nginx-style HTTP index (GEO FTP mirror etc.)."""
    r = get(url)
    if r.status_code != 200:
        return {"url": url, "error": f"HTTP {r.status_code}"}
    out = []
    for m in re.finditer(r'<a href="([^"?#]+)">[^<]*</a>\s*([\d\-]+\s[\d:]+)?\s*([\d.]+[KMGT]?|-)?', r.text):
        href = m.group(1)
        if href.startswith(("/", "http", "..")) or href in ("./",):
            continue
        full = urljoin(url if url.endswith("/") else url + "/", href)
        size = _size_to_bytes(m.group(3)) if m.group(3) and m.group(3) != "-" else None
        out.append({"name": html.unescape(href), "url": full, "size_bytes": size,
                    "role": "dir" if href.endswith("/") else classify(href)})
    return {"url": url, "files": out}


# ---------------------------------------------------------------- GEO
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _geo_dir(acc):
    kind = "series" if acc.startswith("GSE") else "samples"
    stub = acc[:-3] + "nnn" if len(acc) > 6 else acc[:3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/{kind}/{stub}/{acc}/suppl/"


def geo_search(term: str, retmax: int = 200, retstart: int = 0):
    """Search GEO DataSets (db=gds) restricted to series; returns accession, title, n_samples, taxon."""
    r = get(f"{EU}/esearch.fcgi", params={"db": "gds", "term": f"({term}) AND gse[Entry Type]",
                                           "retmax": retmax, "retstart": retstart, "retmode": "json"})
    js = r.json()["esearchresult"]
    ids = js.get("idlist", [])
    res = {"total": int(js.get("count", 0)), "retstart": retstart, "series": []}
    if ids:
        res["series"] = [{k: v for k, v in s.items() if k != "summary"} for s in _gds_summaries(ids)]
    return res


def _gds_summaries(ids):
    out = []
    for i in range(0, len(ids), 100):
        r = get(f"{EU}/esummary.fcgi", params={"db": "gds", "id": ",".join(ids[i:i + 100]), "retmode": "json"})
        js = r.json()["result"]
        for uid in js.get("uids", []):
            d = js[uid]
            out.append({"accession": d.get("accession"), "title": d.get("title"), "taxon": d.get("taxon"),
                        "n_samples": d.get("n_samples"), "gdstype": d.get("gdstype"),
                        "pdat": d.get("pdat"), "summary": d.get("summary", "")[:1500],
                        "suppfile": d.get("suppfile"), "pubmed": d.get("pubmedids"),
                        "samples": [{"gsm": s.get("accession"), "title": s.get("title")}
                                    for s in d.get("samples", [])]})
    return out


def geo_series(gse: str):
    """Series metadata + sample list + series-level supplementary files."""
    r = get(f"{EU}/esearch.fcgi", params={"db": "gds", "term": f"{gse}[ACCN] AND gse[Entry Type]",
                                           "retmode": "json"})
    ids = r.json()["esearchresult"].get("idlist", [])
    meta = _gds_summaries(ids[:1])[0] if ids else {"accession": gse}
    meta["url"] = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
    meta["series_suppl"] = list_directory(_geo_dir(gse)).get("files", [])
    return meta


def geo_all_sample_files(gsms, cap=400):
    return {g: list_directory(_geo_dir(g)).get("files", []) for g in gsms[:cap]}


def geo_sample_files(gsms: list):
    """Supplementary files for up to 60 GSM accessions."""
    return {g: list_directory(_geo_dir(g)).get("files", []) for g in gsms[:60]}


# ---------------------------------------------------------------- Zenodo / Figshare
def zenodo_search(query: str, page: int = 1, size: int = 50):
    r = get("https://zenodo.org/api/records", params={"q": query, "page": page, "size": size,
                                                      "sort": "mostrecent"})
    js = r.json()
    hits = []
    for h in js.get("hits", {}).get("hits", []):
        md = h.get("metadata", {})
        hits.append({"accession": f"zenodo:{h['id']}", "title": md.get("title"), "doi": h.get("doi"),
                     "url": h.get("links", {}).get("self_html"),
                     "description": re.sub("<[^>]+>", " ", md.get("description", ""))[:600],
                     "n_files": len(h.get("files", [])),
                     "total_gb": round(sum(f.get("size", 0) for f in h.get("files", [])) / 1e9, 2)})
    return {"total": js.get("hits", {}).get("total", 0), "page": page, "records": hits}


def zenodo_files(record_id: str):
    rid = str(record_id).replace("zenodo:", "")
    js = get(f"https://zenodo.org/api/records/{rid}").json()
    md = js.get("metadata", {})
    return {"accession": f"zenodo:{rid}", "title": md.get("title"), "doi": js.get("doi"),
            "description": re.sub("<[^>]+>", " ", md.get("description", ""))[:3000],
            "related": [x.get("identifier") for x in md.get("related_identifiers", [])][:20],
            "files": [{"name": f["key"], "url": f["links"]["self"], "size_bytes": f.get("size"),
                       "role": classify(f["key"])} for f in js.get("files", [])]}


def figshare_search(query: str, page: int = 1, page_size: int = 50):
    r = S.post("https://api.figshare.com/v2/articles/search",
               json={"search_for": query, "page": page, "page_size": page_size}, timeout=60)
    return [{"accession": f"figshare:{a['id']}", "title": a.get("title"), "doi": a.get("doi"),
             "url": a.get("url_public_html")} for a in r.json()]


def figshare_files(article_id: str):
    aid = str(article_id).replace("figshare:", "")
    js = get(f"https://api.figshare.com/v2/articles/{aid}").json()
    return {"accession": f"figshare:{aid}", "title": js.get("title"), "doi": js.get("doi"),
            "description": re.sub("<[^>]+>", " ", js.get("description", ""))[:3000],
            "files": [{"name": f["name"], "url": f["download_url"], "size_bytes": f.get("size"),
                       "role": classify(f["name"])} for f in js.get("files", [])]}


# ---------------------------------------------------------------- CELLxGENE
SPATIAL_ASSAYS = ("visium", "slide-seq", "merfish", "xenium", "stereo", "cosmx", "dbit", "hdst",
                  "seqfish", "starmap", "cartana", "resolve", "geomx", "spatial")


def cellxgene_spatial_human():
    js = get("https://api.cellxgene.cziscience.com/curation/v1/datasets", timeout=180).json()
    out = []
    for d in js:
        if not any(o["label"] == "Homo sapiens" for o in d.get("organism", [])):
            continue
        assays = [a["label"] for a in d.get("assay", [])]
        if not any(k in a.lower() for a in assays for k in SPATIAL_ASSAYS):
            continue
        h5ad = [a for a in d.get("assets", []) if a.get("filetype") == "H5AD"]
        out.append({"accession": d["dataset_id"], "collection_id": d["collection_id"],
                    "title": d.get("title"), "collection": d.get("collection_name"),
                    "doi": d.get("collection_doi"), "assay": "; ".join(assays),
                    "tissue": "; ".join(sorted({t["label"] for t in d.get("tissue", [])})),
                    "disease": "; ".join(sorted({t["label"] for t in d.get("disease", [])})),
                    "cell_count": d.get("cell_count"), "explorer_url": d.get("explorer_url"),
                    "h5ad_url": h5ad[0]["url"] if h5ad else None,
                    "h5ad_bytes": h5ad[0].get("filesize") if h5ad else None})
    return out


# ---------------------------------------------------------------- generic HTTP
def http_request(url: str, method: str = "GET", json_body: dict | None = None, max_chars: int = 20000):
    """Fetch any URL / REST endpoint. HTML is stripped to text+links; JSON returned verbatim (truncated)."""
    _throttle(url)
    try:
        if method.upper() == "POST":
            r = S.post(url, json=json_body, timeout=60)
        elif method.upper() == "HEAD":
            r = S.head(url, timeout=60, allow_redirects=True)
            return {"status": r.status_code, "final_url": r.url,
                    "content_length": r.headers.get("Content-Length"),
                    "content_type": r.headers.get("Content-Type")}
        else:
            r = S.get(url, timeout=60, stream=True)
            ct = r.headers.get("Content-Type", "")
            if not any(t in ct for t in ("text", "json", "xml", "javascript")) and ct:
                r.close()
                return {"status": r.status_code, "content_type": ct,
                        "content_length": r.headers.get("Content-Length"),
                        "note": "binary file - not fetched; use HEAD for size"}
    except requests.RequestException as e:
        return {"error": str(e)[:500]}
    ct = r.headers.get("Content-Type", "")
    text = r.text
    if "html" in ct:
        links = re.findall(r'href="([^"#]+)"[^>]*>([^<]{0,80})', text)
        body = re.sub(r"<script.*?</script>|<style.*?</style>", " ", text, flags=re.S)
        body = re.sub(r"\s+", " ", html.unescape(re.sub("<[^>]+>", " ", body)))
        links_txt = "\n".join(f"{urljoin(r.url, h)} | {t.strip()}" for h, t in links[:400])
        text = body[: max_chars // 2] + "\n\nLINKS:\n" + links_txt
    return {"status": r.status_code, "final_url": r.url, "content_type": ct, "text": text[:max_chars]}


# ---------------------------------------------------------------- remote zip (HTTP range requests)
import io
import zipfile


class HttpRangeFile(io.RawIOBase):
    """Seekable read-only file over HTTP Range requests (lets zipfile read only the central directory + members)."""

    def __init__(self, url):
        # 10x CDN stalls on range reads near the end of large files; its S3 origin does not
        url = url.replace("https://cf.10xgenomics.com/", "https://s3-us-west-2.amazonaws.com/10x.files/")
        self.url, self.pos = url, 0
        _throttle(url)
        h = S.head(url, timeout=60, allow_redirects=True)
        self.url = h.url
        if "Content-Length" not in h.headers:        # some servers omit it on HEAD; ask for one byte
            g = S.get(self.url, headers={"Range": "bytes=0-0"}, timeout=60, stream=True)
            cr = g.headers.get("Content-Range", "")
            g.close()
            if "/" not in cr:
                raise IOError("server does not report a size; cannot range-read")
            self.size = int(cr.rsplit("/", 1)[1])
        else:
            self.size = int(h.headers["Content-Length"])
        if h.headers.get("Accept-Ranges", "bytes") == "none":
            raise IOError("server does not support range requests")

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        for attempt in range(6):
            try:
                _throttle(self.url)
                r = S.get(self.url, headers={"Range": f"bytes={self.pos}-{end}"}, timeout=(15, 60))
                if r.status_code == 206:
                    break
            except requests.RequestException:
                pass
            time.sleep(min(60, 3 * 2 ** attempt))
        else:
            raise IOError(f"range read failed {r.status_code}")
        self.pos += len(r.content)
        return r.content

    def readinto(self, b):
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


ZIP_SEP = "#zipmember="


def zip_list(url: str):
    """List members of a remote .zip without downloading it. Member URLs are '<zip url>#zipmember=<path>'."""
    zf = zipfile.ZipFile(io.BufferedReader(HttpRangeFile(url), buffer_size=1 << 16))
    out = [{"name": i.filename, "url": f"{url}{ZIP_SEP}{i.filename}", "size_bytes": i.file_size,
            "role": classify(i.filename)} for i in zf.infolist() if not i.is_dir()]
    return {"zip_url": url, "n_members": len(out), "files": out}


def zip_extract(member_url: str, dest):
    url, name = member_url.split(ZIP_SEP, 1)
    zf = zipfile.ZipFile(io.BufferedReader(HttpRangeFile(url), buffer_size=1 << 20))
    with zf.open(name) as src, open(dest, "wb") as out:
        while chunk := src.read(1 << 22):
            out.write(chunk)


# ---------------------------------------------------------------- Mendeley Data
def mendeley_files(dataset_id: str, version: int | None = None):
    """Mendeley Data public API: metadata + every file (all folders) with direct download URLs."""
    did = re.sub(r".*datasets/", "", str(dataset_id)).split("/")[0].replace("mendeley:", "")
    meta = get(f"https://data.mendeley.com/public-api/datasets/{did}"
               + (f"?version={version}" if version else "")).json()
    if not meta.get("files") and version:      # some datasets only embed files in the unversioned record
        meta = get(f"https://data.mendeley.com/public-api/datasets/{did}").json()
    files = meta.get("files") or get(f"https://data.mendeley.com/public-api/datasets/{did}/files",
                                     params={"folder_id": "root", "version": meta.get("version", 1)}).json()
    out = []
    for f in files:
        cd = f.get("content_details", {})
        name = (f"{f['folder_id'][:8]}/" if f.get("folder_id") not in (None, "root") else "") + f["filename"]
        out.append({"name": name, "url": cd.get("download_url"), "size_bytes": cd.get("size"), "role": classify(f["filename"])})
    return {"accession": f"mendeley:{did}", "version": meta.get("version"), "title": meta.get("name"),
            "doi": (meta.get("doi") or {}).get("id"), "description": (meta.get("description") or "")[:3000], "files": out}
