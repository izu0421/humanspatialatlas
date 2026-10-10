"""One row per sample, joining everything HSA knows about it.

The pieces were scattered: biological labels in sample_meta, QC in sample_qc.csv, panel identity
in panel_registry.csv, bin size and coordinate scale in grid_geometry.csv, the corrected platform
in technology_resolved.csv, study-level context in the datasets table. A consumer of a standardised
.h5ad should not have to re-join six files, so the manifest is built once and embedded in each.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

from . import db
from .config import RUNS

logger = logging.getLogger(__name__)
KEY = ["source", "accession", "sample_id"]
MANIFEST = RUNS / "sample_manifest.csv"

# what a downstream user most often filters on, in the order a reader would want it
CORE = ["source", "accession", "sample_id", "technology", "technology_resolved", "unit", "usable",
        "standardised", "cell_called", "n_cells", "n_genes", "species",
        "tissue", "tissue_ontology_term_id", "disease_state", "disease_level1", "disease_level2",
        "disease_ontology_term_id_level1", "sex", "age_group", "age_num", "age_unit",
        "development_stage", "development_stage_ontology_term_id", "self_reported_ethnicity",
        "self_reported_ethnicity_ontology_term_id", "donor_id", "condition"]


def _read(p: Path, cols=None):
    if not p.exists():
        return pd.DataFrame(columns=KEY + (cols or []))
    d = pd.read_csv(p)
    for c in KEY:
        if c in d.columns:
            d[c] = d[c].astype(str)
    return d


def build() -> pd.DataFrame:
    """Join every per-sample source into one table and write runs/sample_manifest.csv."""
    qc = _read(RUNS / "sample_qc.csv")
    m = qc.copy()

    sm = pd.DataFrame(db.query("SELECT * FROM sample_meta"))
    if len(sm):
        sm = sm.drop(columns=[c for c in ("llm_json", "updated") if c in sm.columns])
        for c in KEY:
            sm[c] = sm[c].astype(str)
        m = m.merge(sm, on=KEY, how="left", suffixes=("", "_meta"))

    ds = pd.DataFrame(db.query("SELECT source, accession, title, url, publication, technology "
                               "AS dataset_technology, verdict FROM datasets"))
    if len(ds):
        for c in ("source", "accession"):
            ds[c] = ds[c].astype(str)
        m = m.merge(ds, on=["source", "accession"], how="left")

    for path, cols, pre in [
        (RUNS / "panel_registry.csv", ["panel_id", "n_genes", "reader"], "panel_"),
        (RUNS / "grid_geometry.csv", ["bin_um", "capture_mm", "grid_rows", "verdict"], "grid_"),
        (RUNS / "technology_resolved.csv", ["resolved", "evidence"], "tech_"),
        (RUNS / "coordinate_scale.csv", ["um_per_unit", "coord_units"], "coord_"),
    ]:
        d = _read(path, cols)
        keep = [c for c in cols if c in d.columns]
        if not keep:
            continue
        d = d[KEY + keep].rename(columns={c: pre + c for c in keep})
        m = m.merge(d, on=KEY, how="left")

    # --- processing state: what exists for this sample, not just what is known about it
    led = RUNS / "materialise_ledger.csv"
    if led.exists():
        L = pd.read_csv(led)
        if len(L):
            for c in KEY:
                L[c] = L[c].astype(str)
            L = L[L.status.isin(["ok", "already_done"])].drop_duplicates(KEY, keep="last")
            cols = [c for c in ["path", "bytes_out", "inferred_sex", "n_obs", "n_vars"]
                    if c in L.columns]
            m = m.merge(L[KEY + cols].rename(columns={
                "path": "h5ad_path", "bytes_out": "h5ad_bytes",
                "n_obs": "h5ad_n_obs", "n_vars": "h5ad_n_vars"}), on=KEY, how="left")
    m["standardised"] = m.get("h5ad_path").notna() if "h5ad_path" in m.columns else False

    # cell-level Visium HD is a separate derived product keyed by its own tag
    try:
        from .visium_hd_cells import HDOUT, _tag
        if HDOUT.exists():
            present = {f.stem for f in HDOUT.glob("*.h5ad")}
            m["hd_cells_path"] = [
                str(HDOUT / f"{_tag(a, s_)}.h5ad") if _tag(a, s_) in present else None
                for a, s_ in zip(m.accession, m.sample_id)]
        else:
            m["hd_cells_path"] = None
    except Exception:
        m["hd_cells_path"] = None
    m["cell_called"] = m.hd_cells_path.notna()

    m["technology_resolved"] = m.get("tech_resolved").where(
        m.get("tech_resolved").notna() & ~m.get("tech_resolved").eq("unknown"), m.technology) \
        if "tech_resolved" in m.columns else m.technology
    # one honest flag for "can a downstream user open this and get coordinates"
    m["usable"] = m.status.astype(str).str.startswith("ok")
    m["unit"] = pd.Series("cell", index=m.index).where(
        ~m.technology_resolved.astype(str).isin(
            ["Visium", "Visium HD", "Slide-seq", "Stereo-seq", "DBiT-seq", "GeoMx", "ST (legacy)"]),
        "spot/bin")
    cols = [c for c in CORE if c in m.columns] + [c for c in m.columns if c not in CORE]
    m = m[cols]
    MANIFEST.parent.mkdir(exist_ok=True, parents=True)
    m.to_csv(MANIFEST, index=False)
    # the master sheet: same content, published alongside the other export tables
    from .config import ROOT
    exp = ROOT / "exports"
    exp.mkdir(exist_ok=True, parents=True)
    m.to_csv(exp / "hsa_master_sheet.csv", index=False)
    logger.info("manifest: %d samples x %d fields -> %s", len(m), m.shape[1], MANIFEST)
    return m


def for_sample(source: str, accession: str, sample_id: str) -> dict:
    """The manifest row for one sample, as a JSON-safe dict for embedding in .h5ad uns."""
    if not MANIFEST.exists():
        build()
    m = pd.read_csv(MANIFEST)
    for c in KEY:
        m[c] = m[c].astype(str)
    r = m[(m.source == str(source)) & (m.accession == str(accession))
          & (m.sample_id == str(sample_id))]
    if not len(r):
        return {}
    d = r.iloc[0].to_dict()
    return {k: (None if pd.isna(v) else (v.item() if hasattr(v, "item") else v))
            for k, v in d.items()}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    m = build()
    print(f"{len(m)} samples, {m.shape[1]} fields")
    print("completeness of the core fields:")
    for c in CORE:
        if c in m.columns:
            known = (~m[c].isin(["unknown", "", None]) & m[c].notna()).sum()
            print(f"  {c:<42} {known:>6,} / {len(m):,}  ({100*known/len(m):.0f}%)")
