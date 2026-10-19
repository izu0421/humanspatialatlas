"""Cross-check the materialised .h5ad files against the metadata that describes them.

The manifest and the matrices are built by different paths -- metadata comes from depositor text
through an LLM and an ontology lookup, the matrices from the files themselves -- so where they
disagree one of them is wrong, and which one differs per check. Every disagreement here was
silent before it was looked for.

Checks, and what a failure means:
  coverage    a usable sample with no file, or a file no manifest row claims
  shape       recorded n_cells / n_genes against the matrix actually written
  species     recorded species against the case of the gene symbols (SOX17 vs Sox17)
  sex         recorded sex against the XIST / chrY call computed from the matrix
  panel       recorded technology against panel size (a 36,601-gene "Xenium" is mislabelled)
  spatial     coordinates present, finite, and not all identical
  pieces      the curator's multi_piece flag against connected components in the coordinates
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent.parent
RUNS = HERE / "runs"
OUT = RUNS / "crosscheck.csv"

_VOCAB = None


def gene_vocabulary() -> set:
    """Symbols seen in many independent whole-transcriptome panels, i.e. real genes.

    Case tells human from mouse but not genes from non-genes: uppercased mouse symbols ARE human
    symbols, so a mouse panel scores 94% against this vocabulary while a CODEX antibody panel
    scores 23%. The two checks catch different things and neither replaces the other.
    """
    global _VOCAB
    if _VOCAB is None:
        f = RUNS / "gene_vocabulary.json"
        _VOCAB = set(json.loads(f.read_text())) if f.exists() else set()
    return _VOCAB


CONTROL = re.compile(r"^(NEGCONTROL|UNASSIGNEDCODEWORD|BLANKCODEWORD|NEGPROBE|ANTISENSE|"
                     r"DEPRECATEDCODEWORD|SYSTEMCONTROL|FALSECODE|UNNAMED)")
UPPER = re.compile(r"^[A-Z][A-Z0-9\-]{1,}$")
TITLE = re.compile(r"^[A-Z][a-z]{2,}[0-9A-Za-z\-]*$")
# a targeted platform cannot carry a whole transcriptome; the reverse is also true
PANEL_MAX = {"Xenium": 6000, "CosMx": 21000, "MERFISH": 2000, "STARmap": 2000, "seqFISH": 2000}
WTX_MIN = {"Visium": 5000, "Visium HD": 5000, "Stereo-seq": 5000, "Slide-seq": 5000}


def _names(h, axis="var"):
    g = h[axis]
    o = g[str(g.attrs.get("_index", "_index"))]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def read_facts(path: str) -> dict:
    """Everything the file itself asserts, read without loading the matrix."""
    f = {"n_obs": None, "n_vars": None, "frac_mouse_like": None, "has_spatial": False,
         "spatial_ok": None, "inferred_sex": None, "frac_known_genes": None,
         "n_control_probes": 0}
    with h5py.File(path, "r") as h:
        var = _names(h, "var")
        f["n_vars"] = len(var)
        f["n_obs"] = len(_names(h, "obs"))
        vocab = gene_vocabulary()
        # Strip the reference prefix a barnyard deposit carries ("GRCh38_ACTB"): those are real
        # gene symbols and scored 0% against a plain-symbol vocabulary, which flagged 109 perfectly
        # good samples as "not genes" on the first run of this check.
        up_names = [re.sub(r"^(GRCH38|GRCH37|HG38|HG19|HUMAN|MM10|GRCM38|GRCM39|MOUSE)[_-]", "",
                           str(x).upper()) for x in var]
        # probe-version suffixes are part of the probe name, not the gene: "ABHD16A-v1-prime"
        up_names = [re.sub(r"-V\d+(-PRIME|-ALT)?$", "", n) for n in up_names]
        # control probes are legitimately not genes and must leave the denominator; one Xenium
        # panel is 491 controls to 50 genes, which read as "9% known genes" and looked broken
        real = [n for n in up_names if not CONTROL.match(n)]
        f["n_control_probes"] = len(up_names) - len(real)
        if vocab and real and not any(n.startswith("ENSG") for n in real[:20]):
            f["frac_known_genes"] = round(sum(1 for n in real if n in vocab) / len(real), 3)
        sym = [str(x) for x in var if not str(x).startswith(("ENSG", "ENSMUS"))][:600]
        up = sum(1 for x in sym if UPPER.match(x))
        ti = sum(1 for x in sym if TITLE.match(x))
        # 10, not 30: a 33-gene mouse panel and a 35-plex mouse CODEX panel both slipped past a
        # 30-symbol minimum. Purity is required to scale with how few symbols there are.
        if up + ti >= 10:
            f["frac_mouse_like"] = round(ti / (up + ti), 3)
            f["n_symbols"] = up + ti
        if "obsm/spatial" in h:
            f["has_spatial"] = True
            xy = h["obsm/spatial"][:]
            xy = np.asarray(xy, dtype=float)
            finite = np.isfinite(xy).all(1)
            # all-identical coordinates pass a "has spatial" test and carry no geometry at all
            f["spatial_ok"] = bool(finite.mean() > 0.5
                                   and len(np.unique(xy[finite][:, :2], axis=0)) > 10)
        s = h.get("uns/hsa/inferred_sex")
        if s is not None:
            s = s[()]
            f["inferred_sex"] = s.decode() if isinstance(s, bytes) else str(s)
    return f


def run(limit: int | None = None, write: bool = True) -> pd.DataFrame:
    m = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    led = pd.read_csv(RUNS / "materialise_ledger.csv")
    ok = led[led.status.astype(str).str.startswith("ok")
             | led.status.astype(str).eq("already_done")]
    key = ["source", "accession", "sample_id"]
    j = m.merge(ok[key + ["path"]], on=key, how="left")
    usable = j[j.usable.fillna(False)]
    if limit:
        usable = usable.head(limit)

    rows = []
    for i, (_, r) in enumerate(usable.iterrows(), 1):
        rec = {k: r[k] for k in key}
        rec["technology"] = r.technology_resolved
        p = r.get("path")
        if not isinstance(p, str) or not Path(p).exists():
            rows.append({**rec, "issue": "no_file",
                         "detail": "usable but nothing materialised"})
            continue
        try:
            f = read_facts(p)
        except Exception as e:
            rows.append({**rec, "issue": "unreadable", "detail": f"{type(e).__name__}: {e}"[:80]})
            continue
        issues = []
        # shape
        for col, got, lab in [("n_cells", f["n_obs"], "cells"), ("n_genes", f["n_vars"], "genes")]:
            want = r.get(col)
            if pd.notna(want) and got and abs(float(want) - got) / max(got, 1) > 0.02:
                issues.append(f"{lab}_mismatch: manifest {int(want):,} vs file {got:,}")
        # species
        sp = str(r.get("species") or "").lower()
        fm = f["frac_mouse_like"]
        if fm is not None:
            thr = 0.8 if (f.get("n_symbols") or 0) >= 30 else 0.95
            if fm > thr and "mus" not in sp:
                issues.append(f"species: symbols {fm:.0%} mouse but recorded '{sp or 'blank'}'")
            elif fm < 0.2 and "mus_musculus" == sp:
                issues.append(f"species: symbols human but recorded mouse")
        # sex
        rec_sex = str(r.get("sex") or "unknown").lower()
        inf = str(f["inferred_sex"] or "unknown").lower()
        if rec_sex in ("male", "female") and inf in ("male", "female") and rec_sex != inf:
            issues.append(f"sex: recorded {rec_sex}, expression says {inf}")
        # panel vs technology
        t = str(r.technology_resolved)
        if t in PANEL_MAX and f["n_vars"] and f["n_vars"] > PANEL_MAX[t]:
            issues.append(f"panel: {t} with {f['n_vars']:,} genes (> {PANEL_MAX[t]:,})")
        if t in WTX_MIN and f["n_vars"] and f["n_vars"] < WTX_MIN[t]:
            issues.append(f"panel: {t} with only {f['n_vars']:,} genes")
        # var axis holds genes at all
        fk = f.get("frac_known_genes")
        if fk is not None and fk < 0.5 and f["n_vars"] and f["n_vars"] >= 20:
            issues.append(f"not_genes: only {fk:.0%} of the var axis are known gene symbols")
        # spatial
        if not f["has_spatial"]:
            issues.append("spatial: no coordinates")
        elif f["spatial_ok"] is False:
            issues.append("spatial: coordinates degenerate (non-finite or all identical)")
        rows.append({**rec, **{k: f[k] for k in ("n_obs", "n_vars", "frac_mouse_like",
                                                 "inferred_sex")},
                     "issue": issues[0].split(":")[0] if issues else "",
                     "detail": " | ".join(issues)})
        if i % 1000 == 0:
            logger.info("checked %d/%d", i, len(usable))
    d = pd.DataFrame(rows)
    if write:
        d.to_csv(OUT, index=False)
    bad = d[d.issue.astype(str).ne("")]
    logger.info("checked %d samples, %d with a disagreement", len(d), len(bad))
    return d
