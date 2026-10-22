"""Read-only consistency review of every materialised sample.

The atlas is assembled from ~1,500 deposits in a dozen formats, and each loader fix so far has
been found by one sample behaving oddly rather than by asking the whole corpus the same question.
This asks the same question of every file, which is the only way to see that 212 samples index by
Ensembl while 7,500 index by symbol, or that one obs column exists in 94% of files and is
therefore not safe to rely on.

Nothing here writes. Severity is the point of the output:

  error    the file cannot be used as it claims to be -- no coordinates, no genes, a var axis of
           positional integers, an obs index that is not cell identity
  warn     usable but will surprise someone -- a gene dialect that is not the corpus norm, a
           coordinate unit that is not microns, a matrix that is not counts
  info     recorded for completeness

A check is only worth having if it can fail, so each one states what it would catch.
"""
from __future__ import annotations

import logging
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

import h5py
import numpy as np
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)
OUT = RUNS / "review_samples.csv"
REPORT = RUNS / "review_report.md"

ENSEMBL = re.compile(r"^ENS[A-Z]*G\d{6,}")
SYMBOL = re.compile(r"^[A-Za-z][A-Za-z0-9\-\._]{0,20}$")
POSITIONAL = re.compile(r"^\d+$")


def _names(h, axis: str) -> list[str]:
    g = h.get(axis)
    if g is None:
        return []
    idx = g.attrs.get("_index", "_index")
    idx = idx.decode() if isinstance(idx, bytes) else str(idx)
    d = g.get(idx)
    if d is None:
        return []
    if isinstance(d, h5py.Group):                       # categorical or nullable-string
        d = d.get("categories") if "categories" in d else d.get("values")
        if d is None:
            return []
    out = d[: min(5000, d.shape[0])]
    return [x.decode() if isinstance(x, bytes) else str(x) for x in out]


def _cols(h, axis: str) -> list[str]:
    g = h.get(axis)
    if g is None:
        return []
    order = g.attrs.get("column-order", [])
    out = [x.decode() if isinstance(x, bytes) else str(x) for x in order]
    return out or [k for k in g.keys() if not str(k).startswith("_")]


def review_one(path: str) -> dict:
    r: dict = {"path": str(path), "errors": [], "warnings": []}
    with h5py.File(path, "r") as h:
        # --- shape
        obs_names, var_names = _names(h, "obs"), _names(h, "var")
        X = h.get("X")
        if X is None:
            r["errors"].append("no X")
            return r
        if isinstance(X, h5py.Group):
            r["X_format"] = "sparse"
            r["X_dtype"] = str(X["data"].dtype)
            sh = X.attrs.get("shape")
            r["n_obs"], r["n_vars"] = (int(sh[0]), int(sh[1])) if sh is not None else (None, None)
        else:
            r["X_format"] = "dense"
            r["X_dtype"] = str(X.dtype)
            r["n_obs"], r["n_vars"] = int(X.shape[0]), int(X.shape[1])

        # --- coordinates: present, right shape, finite, and actually varying
        sp = h.get("obsm/spatial")
        if sp is None:
            r["errors"].append("no obsm/spatial")
            r["coord_ok"] = False
        else:
            xy = np.asarray(sp[: min(50_000, sp.shape[0])], dtype=float)
            r["coord_dims"] = int(xy.shape[1]) if xy.ndim == 2 else 0
            if xy.ndim != 2 or xy.shape[1] < 2:
                r["errors"].append(f"obsm/spatial has {r['coord_dims']} columns")
                r["coord_ok"] = False
            else:
                fin = np.isfinite(xy[:, :2]).all(1)
                r["coord_frac_finite"] = round(float(fin.mean()), 4)
                r["coord_n_unique"] = int(np.unique(xy[fin][:, :2], axis=0).shape[0]) \
                    if fin.any() else 0
                if fin.mean() < 0.5:
                    r["errors"].append(f"coordinates {100*(1-fin.mean()):.0f}% non-finite")
                elif r["coord_n_unique"] < max(10, 0.05 * max(1, fin.sum())):
                    r["errors"].append("coordinates degenerate (few distinct positions)")
                r["coord_ok"] = not r["errors"]
                if r["n_obs"] and sp.shape[0] != r["n_obs"]:
                    r["errors"].append(f"obsm/spatial rows {sp.shape[0]} != n_obs {r['n_obs']}")

        # --- gene axis: identity, dialect, uniqueness
        if not var_names:
            r["errors"].append("no var names")
        else:
            n = len(var_names)
            ens = sum(1 for v in var_names if ENSEMBL.match(v))
            pos = sum(1 for v in var_names if POSITIONAL.match(v))
            sym = sum(1 for v in var_names if SYMBOL.match(v) and not ENSEMBL.match(v))
            r["gene_frac_ensembl"] = round(ens / n, 3)
            r["gene_frac_positional"] = round(pos / n, 3)
            r["gene_frac_symbol"] = round(sym / n, 3)
            r["gene_dialect"] = ("ensembl" if ens / n > .5 else
                                 "positional" if pos / n > .5 else
                                 "symbol" if sym / n > .5 else "mixed")
            r["gene_has_version"] = bool(any("." in v and ENSEMBL.match(v) for v in var_names))
            r["gene_dup_frac"] = round(1 - len(set(var_names)) / n, 4)
            # a var axis of row numbers is not gene identity, and joins to nothing
            if pos / n > 0.5:
                r["errors"].append("var axis is positional integers, not gene identity")
            if r["gene_dup_frac"] > 0.01:
                r["warnings"].append(f"{100*r['gene_dup_frac']:.0f}% duplicate gene names")
            if r["gene_dialect"] == "mixed":
                r["warnings"].append("mixed gene identifier dialects on one axis")
            r["var_cols"] = ";".join(_cols(h, "var")[:12])

        # --- obs axis: identity and columns
        if obs_names:
            m = len(obs_names)
            r["obs_frac_positional"] = round(
                sum(1 for v in obs_names if POSITIONAL.match(v)) / m, 3)
            r["obs_dup_frac"] = round(1 - len(set(obs_names)) / m, 4)
            if r["obs_frac_positional"] > 0.5:
                r["warnings"].append("obs index is positional, not barcodes or cell ids")
            if r["obs_dup_frac"] > 0.01:
                r["warnings"].append(f"{100*r['obs_dup_frac']:.0f}% duplicate obs names")
        r["obs_cols"] = ";".join(_cols(h, "obs")[:20])

        # --- what the pipeline recorded about itself
        u = h.get("uns/hsa")
        if u is None:
            r["warnings"].append("no uns/hsa block")
        else:
            def g(k):
                if k not in u:
                    return None
                v = u[k][()]
                if isinstance(v, bytes):
                    return v.decode()
                # h5py returns an Empty object for a null scalar, and it is unhashable -- it
                # passed straight through into a value_counts() and broke the report after the
                # per-sample pass had already finished.
                if isinstance(v, h5py.Empty) or v is None:
                    return None
                return v if np.isscalar(v) else str(v)
            for k in ("hsa_schema", "matrix_kind", "coord_coord_units", "technology_resolved",
                      "celltype_model", "coord_source"):
                r[k] = g(k)
            if r.get("matrix_kind") not in (None, "raw_counts"):
                r["warnings"].append(f"matrix is {r['matrix_kind']}, not counts")
            if r.get("coord_coord_units") not in (None, "microns"):
                r["warnings"].append(f"coordinates in {r['coord_coord_units']}")
    r["n_errors"], r["n_warnings"] = len(r["errors"]), len(r["warnings"])
    r["errors"] = " | ".join(r["errors"])
    r["warnings"] = " | ".join(r["warnings"])
    return r


def _safe(path):
    try:
        return review_one(path)
    except Exception as e:
        return {"path": path, "errors": f"unreadable: {type(e).__name__}: {str(e)[:90]}",
                "warnings": "", "n_errors": 1, "n_warnings": 0}


def run(workers: int = 8, limit: int | None = None) -> pd.DataFrame:
    m = pd.read_csv(RUNS / "sample_manifest.csv", low_memory=False)
    m = m[m.h5ad_path.notna()]
    if limit:
        m = m.head(limit)
    key = ["source", "accession", "sample_id", "technology_resolved", "tissue"]
    meta = m[key + ["h5ad_path"]].rename(columns={"h5ad_path": "path"})
    logger.info("reviewing %d samples with %d workers", len(meta), workers)

    rows = []
    with ProcessPoolExecutor(workers) as ex:
        futs = {ex.submit(_safe, p): p for p in meta.path}
        for i, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if i % 1000 == 0:
                logger.info("  %d/%d", i, len(futs))
    d = pd.DataFrame(rows).merge(meta, on="path", how="left")
    d.to_csv(OUT, index=False)
    logger.info("wrote %s", OUT)
    return d


# ---------------------------------------------------------------- corpus-level consistency
# A per-sample check sees one file and cannot tell a corpus-wide convention from an exception.
# These ask what is NORMAL across the atlas and name what departs from it, which is how a
# minority dialect (212 Ensembl-indexed samples among 7,500 symbol-indexed ones) becomes visible.

UNIVERSAL = 0.99            # present in at least this fraction = safe to rely on
COMMON = 0.80


def consistency(d: pd.DataFrame) -> dict:
    """What the corpus agrees on, and which samples depart from it."""
    ok = d[d.n_errors.eq(0)] if "n_errors" in d else d
    n = max(1, len(ok))
    out: dict = {"n_reviewed": len(d), "n_usable": len(ok)}

    # obs columns: which are safe to rely on, and which exist in only part of the corpus
    cnt = Counter()
    for s in ok.obs_cols.fillna(""):
        cnt.update([c for c in str(s).split(";") if c])
    cols = pd.DataFrame({"column": list(cnt), "n": list(cnt.values())})
    cols["frac"] = cols.n / n
    out["obs_universal"] = sorted(cols.loc[cols.frac >= UNIVERSAL, "column"])
    out["obs_partial"] = cols[(cols.frac < UNIVERSAL) & (cols.frac >= 0.01)] \
        .sort_values("frac", ascending=False).head(25).to_dict("records")

    for field, label in (("gene_dialect", "gene_dialect"), ("X_dtype", "X_dtype"),
                         ("X_format", "X_format"), ("matrix_kind", "matrix_kind"),
                         ("coord_coord_units", "coord_units"), ("hsa_schema", "hsa_schema")):
        if field in ok.columns:
            vc = ok[field].fillna("(absent)").value_counts()
            out[label] = vc.to_dict()
            # the minority convention is the one that breaks a downstream join
            if len(vc) > 1 and vc.iloc[0] / n < UNIVERSAL:
                out.setdefault("_split", []).append(
                    {"field": label, "majority": str(vc.index[0]),
                     "majority_frac": round(vc.iloc[0] / n, 3),
                     "minority": {str(k): int(v) for k, v in vc.iloc[1:6].items()}})
    return out


def report(d: pd.DataFrame | None = None, path=REPORT) -> str:
    import json
    if d is None:
        d = pd.read_csv(OUT, low_memory=False)
    c = consistency(d)
    L = ["# Atlas consistency review", "",
         f"{c['n_reviewed']:,} samples reviewed, {c['n_usable']:,} with no structural error.", ""]

    err = d[d.n_errors > 0]
    L += ["## Structural errors", ""]
    if not len(err):
        L.append("None.")
    else:
        ec = Counter()
        for s in err.errors.fillna(""):
            for p in str(s).split("|"):
                if p.strip():
                    ec[p.strip()[:70]] += 1
        L.append(f"{len(err):,} samples ({100*len(err)/max(1,len(d)):.1f}%).")
        L.append("")
        L += [f"- **{v}** — {k}" for k, v in ec.most_common(15)]
    L.append("")

    warn = d[d.n_warnings > 0]
    L += ["## Warnings", "", f"{len(warn):,} samples carry at least one.", ""]
    wc = Counter()
    for s in warn.warnings.fillna(""):
        for p in str(s).split("|"):
            if p.strip():
                wc[p.strip()[:70]] += 1
    L += [f"- **{v}** — {k}" for k, v in wc.most_common(15)]
    L.append("")

    L += ["## What the corpus agrees on", ""]
    for f in ("gene_dialect", "X_dtype", "X_format", "matrix_kind", "coord_units", "hsa_schema"):
        if f in c:
            inner = ", ".join(f"{k} {v:,}" for k, v in list(c[f].items())[:6])
            L.append(f"- **{f}**: {inner}")
    L.append("")

    if c.get("_split"):
        L += ["### Fields the corpus does NOT agree on", "",
              "A minority convention is what breaks a downstream join, because code written "
              "against the majority passes every test until it meets one of these.", ""]
        for s in c["_split"]:
            L.append(f"- **{s['field']}** — {s['majority_frac']:.0%} are `{s['majority']}`; "
                     f"the rest: {s['minority']}")
        L.append("")

    L += ["## obs columns", "",
          f"Present in >= {UNIVERSAL:.0%} of usable samples (safe to rely on):", ""]
    L.append("`" + "`, `".join(c["obs_universal"]) + "`" if c["obs_universal"] else "_none_")
    L += ["", "Present in only part of the corpus (a join on these silently drops samples):", ""]
    for r in c["obs_partial"][:15]:
        L.append(f"- `{r['column']}` — {r['frac']:.1%} ({r['n']:,})")
    L.append("")

    txt = "\n".join(L)
    path.write_text(txt)
    (RUNS / "review_consistency.json").write_text(json.dumps(c, indent=1, default=str))
    return txt


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    import sys
    d = run(workers=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
    print(report(d))
