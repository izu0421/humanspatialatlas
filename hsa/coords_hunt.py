"""Ask the model where a deposit's spatial coordinates are, when the loader could not find them.

1,202 samples across 307 studies load their expression matrix but have no coordinates. The
loader searches for file-name patterns it knows; a depositor who named the file something else,
put the coordinates inside the matrix object, or described them only in the paper is invisible
to it. That is a reading problem, not a parsing problem, so it goes to the model: show it every
file in the deposit plus the dataset text, and ask which file holds the positions and how to read
it.

The answer is a HYPOTHESIS, not a fact. It names a file and a way to read it, which the loader
can then try and either succeed or fail on -- so a wrong answer costs one load attempt and is
caught, rather than silently producing a sample with invented positions.
"""
from __future__ import annotations

import json
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from . import db
from .config import MODEL, PRICE, RUNS
from .metadata import CLIENT, page_text

logger = logging.getLogger(__name__)

SCHEMA = {
    "type": "object", "additionalProperties": False,
    "required": ["accession", "verdict", "coord_file", "how_to_read", "confidence", "evidence"],
    "properties": {
        "accession": {"type": "string"},
        "verdict": {"type": "string",
                    "enum": ["file_in_deposit", "inside_matrix_file", "in_another_accession",
                             "paper_only", "not_available"],
                    "description": "where the coordinates are, if anywhere"},
        "coord_file": {"type": "string",
                       "description": "the exact filename from the list that holds positions, "
                                      "or the key/slot inside the matrix file; '' if none"},
        "how_to_read": {"type": "string",
                        "description": "concretely: which columns are x and y, whether there is a "
                                       "header, which column is the barcode, any scale factor or "
                                       "units mentioned; '' if unknown"},
        "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
        "evidence": {"type": "string", "description": "what you based this on, briefly"}}}

SYSTEM = """You locate spatial coordinates for a spatial-transcriptomics deposit whose coordinate \
file our loader could not identify.

You are given: the dataset title and description, its publication if known, and the FULL list of \
files in the deposit with their roles and sizes. Decide where each spot's or cell's x/y position \
is, and say exactly how to read it.

Coordinates hide in these places, in roughly this order of likelihood:
  - a file whose name does not match the usual patterns: 'positions', 'coords', 'centroids', \
'spatial', 'locations', 'barcode_pos', 'tissue_positions', 'metadata', 'obs', 'cell_info', \
'spot_data', or a bare 'metadata.csv' holding x/y among other columns
  - INSIDE the matrix file: an .h5ad's obsm['spatial'] or obs columns; a Seurat object's images \
slot; a .loom column attribute; an .mtx deposit's barcode names encoding position as \
'12x34' or 'A1_12_34'
  - the barcode itself: Slide-seq puck barcodes and some Stereo-seq deposits encode x_y directly
  - a separate accession or a supplementary file of the PAPER rather than the data deposit
  - genuinely absent: some depositors upload only the expression matrix

Rules. Name a file from the list EXACTLY as given; do not invent filenames. If the positions are \
inside the matrix file, set coord_file to the slot ('obsm/spatial', 'obs: x,y') and verdict to \
inside_matrix_file. If the coordinates exist only in a paper figure or supplementary PDF and not \
as data, say paper_only -- that is useful to know and not a failure. If nothing in the deposit \
can hold positions, say not_available rather than guessing: a wrong file name costs a load \
attempt, but an invented one that happens to parse would put fictional positions into the atlas."""


def _files_for(source, accession, limit=60):
    rows = db.query("SELECT url, role, fmt, size_bytes, local_path FROM files "
                    "WHERE source=? AND accession=?", (source, accession))
    out = []
    for r in rows[:limit]:
        name = str(r["local_path"] or r["url"] or "").split("/")[-1]
        out.append({"file": name[:90], "role": r["role"], "fmt": r["fmt"],
                    "MB": round((r["size_bytes"] or 0) / 1e6, 1)})
    return out, len(rows)


def one(ds, run):
    files, n_total = _files_for(ds["source"], ds["accession"])
    ctx = {"accession": ds["accession"], "source": ds["source"],
           "title": ds.get("title"), "publication": ds.get("publication"),
           "url": ds.get("url"), "n_files_total": n_total, "files": files}
    if ds["source"] != "GEO":
        ctx["page_text"] = page_text(ds, n=4000)
    with CLIENT.messages.stream(
            model=MODEL, max_tokens=4000, system=SYSTEM,
            messages=[{"role": "user", "content": json.dumps(ctx, default=str)[:120000]}],
            output_config={"effort": "low",
                           "format": {"type": "json_schema", "schema": SCHEMA}}) as st:
        msg = st.get_final_message()
    u = msg.usage
    step = {"in": u.input_tokens or 0, "out": u.output_tokens or 0,
            "cache_read": u.cache_read_input_tokens or 0,
            "cache_write": u.cache_creation_input_tokens or 0}
    db.log_usage(run, step, (step["in"] * PRICE["in"] + step["out"] * PRICE["out"]
                             + step["cache_read"] * PRICE["cache_read"]
                             + step["cache_write"] * PRICE["cache_write"]) / 1e6)
    r = json.loads(msg.content[0].text)
    r["accession"] = ds["accession"]
    r["source"] = ds["source"]
    return r


def run(limit=None, workers=4, budget=20.0):
    q = pd.read_csv(RUNS / "sample_qc.csv", low_memory=False).drop_duplicates(
        ["source", "accession", "sample_id"], keep="last")
    nc = q[q.status.astype(str).str.startswith("no_coords")]
    counts = nc.groupby(["source", "accession"]).size().sort_values(ascending=False)
    ds_rows = {(r["source"], r["accession"]): r
               for r in db.query("SELECT source, accession, title, url, publication FROM datasets")}
    todo = [ds_rows[k] for k in counts.index if k in ds_rows]
    if limit:
        todo = todo[:limit]
    logger.info("%d studies with no coordinates (%d samples)", len(todo), int(counts.sum()))
    out = []
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(one, d, f"coords_{d['accession']}"[:120]): d for d in todo}
        for i, f in enumerate(as_completed(futs), 1):
            d = futs[f]
            try:
                r = f.result(); out.append(r)
                logger.info("[%d/%d] %s -> %s (%s) %s", i, len(todo), d["accession"][:34],
                            r["verdict"], r["confidence"], str(r["coord_file"])[:40])
            except Exception as e:
                logger.warning("[%d/%d] %s FAILED %s", i, len(todo), d["accession"][:34],
                               str(e)[:80])
            if i % 20 == 0:
                pd.DataFrame(out).to_csv(RUNS / "coords_hunt.csv", index=False)
            spent = db.query("SELECT COALESCE(SUM(usd),0) s FROM usage "
                             "WHERE run LIKE 'coords_%'")[0]["s"]
            if spent > budget:
                logger.warning("budget reached (%.2f)", spent)
                break
    d = pd.DataFrame(out)
    d.to_csv(RUNS / "coords_hunt.csv", index=False)
    if len(d):
        logger.info("verdicts: %s", d.verdict.value_counts().to_dict())
    return d
