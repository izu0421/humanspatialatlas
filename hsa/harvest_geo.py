"""Deterministic GEO candidate harvest (no LLM): union of technology queries, human-only series."""
from . import db, tools_api as T

TERMS = ["visium", "\"spatial transcriptomics\"", "\"spatial transcriptome\"", "\"spatially resolved\"", "xenium",
         "cosmx", "merscope", "merfish", "\"stereo-seq\"", "stomics", "\"slide-seq\"", "slideseq", "curio",
         "dbit-seq", "geomx", "seqfish", "starmap", "hdst", "\"spatial gene expression\"", "\"visium hd\"",
         "\"in situ sequencing\"", "\"spatial omics\"", "\"spatial profiling\""]


def run():
    seen = {}
    for t in TERMS:
        term = f'{t} AND "Homo sapiens"[Organism]'
        start, total = 0, None
        while total is None or start < total:
            r = T.geo_search(term, retmax=500, retstart=start)
            total = r["total"]
            for s in r["series"]:
                seen.setdefault(s["accession"], {"accession": s["accession"], "title": s["title"],
                                                 "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={s['accession']}",
                                                 "terms": set()})["terms"].add(t)
            start += 500
        print(f"{t}: {total}  (union {len(seen)})", flush=True)
    for v in seen.values():
        v["reason"] = "geo terms: " + ",".join(sorted(v.pop("terms")))
    n = sum(db.add_candidates("GEO", [v], v["reason"]) for v in seen.values())
    db.upsert_source("GEO", "https://www.ncbi.nlm.nih.gov/geo/",
                     "Entrez E-utilities (db=gds) + HTTPS FTP mirror ftp.ncbi.nlm.nih.gov/geo/{series,samples}/.../suppl/",
                     len(seen), "candidates from deterministic keyword union; LLM curator decides spatial/human/files")
    print(f"queued {n} new GEO candidates ({len(seen)} total)")
