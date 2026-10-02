"""HSA - Human Spatial Atlas pipeline.

  python run_hsa.py harvest              # deterministic: CELLxGENE spatial human datasets
  python run_hsa.py scout [LANE ...]     # Sonnet scouts per database (all lanes if none given)
  python run_hsa.py curate [-n N]        # Sonnet curates queued candidates (parallel)
  python run_hsa.py download             # fetch matrix + coords files
  python run_hsa.py bundles              # extract matrix + coords members from tar/zip bundles
  python run_hsa.py status
"""
import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from hsa import db, tools_api as T


def harvest():
    rows = T.cellxgene_spatial_human()
    db.upsert_source("CELLxGENE", "https://cellxgene.cziscience.com/datasets",
                     "Curation API https://api.cellxgene.cziscience.com/curation/v1/datasets ; H5AD asset per dataset "
                     "(obsm['spatial'] + counts in raw.X/X)", len(rows), "harvested deterministically")
    for r in rows:
        db.record_dataset({
            "source": "CELLxGENE", "accession": r["accession"], "title": r["title"], "url": r["explorer_url"],
            "publication": r["doi"] or "", "technology": r["assay"], "tissue": r["tissue"], "disease": r["disease"],
            "n_samples": 1, "verdict": "ready" if r["h5ad_url"] else "no_processed_data",
            "notes": f"collection {r['collection']} ({r['collection_id']}); cells={r['cell_count']}",
            "links": [f"https://cellxgene.cziscience.com/collections/{r['collection_id']}"],
            "samples": [{"sample_id": r["accession"], "files": [
                {"url": r["h5ad_url"], "role": "matrix+coords", "fmt": "h5ad", "size_bytes": r["h5ad_bytes"]}]}]
            if r["h5ad_url"] else []}, curated_by="harvest")
    print(f"CELLxGENE: {len(rows)} human spatial datasets recorded")


def status():
    for q in ["SELECT name, est_human_datasets, access_method FROM sources ORDER BY name",
              "SELECT source, status, COUNT(*) n FROM candidates GROUP BY 1,2",
              "SELECT source, verdict, COUNT(*) n FROM datasets GROUP BY 1,2 ORDER BY 1,3 DESC",
              "SELECT role, dl_status, COUNT(*) n, ROUND(SUM(size_bytes)/1e9,1) gb FROM files GROUP BY 1,2",
              "SELECT ROUND(SUM(usd),2) usd, SUM(input) inp, SUM(output) outp FROM usage"]:
        print("\n" + q)
        for r in db.query(q):
            print("  ", {k: (str(v)[:90] if isinstance(v, str) else v) for k, v in r.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    ap.add_argument("lanes", nargs="*")
    ap.add_argument("-n", type=int, default=None, help="max candidates to curate")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--budget", type=float, default=50.0, help="global $ cap across all agent calls")
    ap.add_argument("--source", default=None)
    a = ap.parse_args()

    if a.cmd == "harvest":
        harvest()
    elif a.cmd == "scout":
        from hsa import roles
        lanes = a.lanes or list(roles.LANES)
        with ThreadPoolExecutor(min(a.workers, len(lanes))) as ex:
            futs = {ex.submit(roles.scout, l, a.budget): l for l in lanes}
            for fu in as_completed(futs):
                try:
                    final, usage, rec = fu.result()
                    print(f"\n=== {futs[fu]} === recorded {len(rec)}\n{final}\n{usage}", flush=True)
                except Exception as e:
                    print(f"\n=== {futs[fu]} FAILED: {e}", flush=True)
    elif a.cmd == "curate":
        from hsa import roles
        sql = "SELECT * FROM candidates WHERE status='pending'"
        if a.source:
            sql += f" AND source='{a.source}'"
        cands = db.query(sql + " ORDER BY added" + (f" LIMIT {a.n}" if a.n else ""))
        print(f"curating {len(cands)}; spent so far ${db.total_usd():.2f}", flush=True)
        with ThreadPoolExecutor(a.workers) as ex:
            futs = {ex.submit(roles.curate, c, a.budget): c for c in cands}
            for i, fu in enumerate(as_completed(futs)):
                c = futs[fu]
                try:
                    final, _ = fu.result()
                    print(f"[{i+1}/{len(cands)}] {c['accession']}: {final.strip()[:150]}", flush=True)
                except Exception as e:
                    print(f"[{i+1}/{len(cands)}] {c['accession']} ERROR {type(e).__name__}: {e}", flush=True)
                    if "budget" in str(e):
                        ex.shutdown(cancel_futures=True)
                        break
        print(f"total spent ${db.total_usd():.2f}")
    elif a.cmd == "download":
        from hsa import download
        download.run(workers=a.workers, source=a.source, limit=a.n)
    elif a.cmd == "bundles":
        from hsa import download
        download.run_bundles(workers=a.workers, source=a.source, limit=a.n)
    elif a.cmd == "status":
        status()
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
