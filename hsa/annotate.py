"""Organ-specific cell-type annotation for the cell-resolution part of the atlas.

Why a model per organ AND per panel, rather than one pre-trained CellTypist model:

  * The published CellTypist models are trained on whole transcriptomes. A Xenium section carries
    477 genes and a MERFISH section 315, so applying a 20,000-gene model means imputing >97% of
    its input as zero. The model then scores on absence, and absence in a targeted panel means
    "not on the panel", not "not expressed" -- the single most common silent failure in this
    corpus. Training on the panel's own genes keeps every feature a measured one.
  * Cell types are organ-specific. A pan-tissue model spends its capacity separating hepatocytes
    from neurons; an organ model spends it on the distinctions that are actually in the section.

Reference data is CELLxGENE Census: uniformly processed, ontology-labelled, and queryable per
tissue, which is what makes an organ-by-organ reference practical at all. Labels are kept at
`cell_type` (the ontology term) so annotations are comparable across studies.

Each trained model is cached by (organ, panel hash), because a panel is shared across many
sections -- 1,261 Xenium sections use a handful of distinct panels, so the training cost is paid
per panel, not per section.
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .config import RUNS

logger = logging.getLogger(__name__)
MODELS = RUNS.parent / "celltype_models"
OUT = RUNS / "annotation_log.csv"

CELL_TECH = ("Xenium", "CosMx", "MERFISH", "Atera", "cell-level WTA", "STARmap", "seqFISH")
MIN_GENES = 50            # below this a panel cannot support a multi-class model
MIN_REF_CELLS = 500       # per organ
# 40,000 balanced cells, not 120,000. The model has at most a few hundred features (a panel),
# so the extra rows buy very little and cost a lot: the 97,816-cell lung fetch died part way
# through pulling its matrix. A bounded fetch that completes is worth more than a large one
# that does not.
MAX_REF_CELLS = 40_000    # per organ, balanced across types
MIN_CELLS_PER_TYPE = 50


def panel_key(genes: list[str]) -> str:
    """Stable id for a gene panel, so one model serves every section that shares it."""
    return hashlib.md5("|".join(sorted(set(map(str, genes)))).encode()).hexdigest()[:12]


# Census uses UBERON-backed tissue_general; the atlas stores free-text tissue. Map the organs
# that actually carry cell-resolution sections here, and fall back to a pan-tissue model when a
# sample's organ has no usable reference.
ORGAN_MAP = {
    "lung": "lung", "brain": "brain", "dorsolateral_prefrontal_cortex": "brain",
    "cerebral_cortex": "brain", "prefrontal_cortex": "brain",
    "skin_of_body": "skin of body", "zone_of_skin": "skin of body",
    "breast": "breast", "kidney": "kidney", "pancreas": "pancreas",
    "colon": "large intestine", "rectum": "large intestine", "large_intestine": "large intestine",
    "small_intestine": "small intestine", "liver": "liver", "heart": "heart",
    "bone_marrow": "bone marrow", "lymph_node": "lymph node", "spleen": "spleen",
    "prostate_gland": "prostate gland", "uterus": "uterus", "ovary": "ovary",
    "stomach": "stomach", "esophagus": "esophagus", "thymus": "thymus",
    "placenta": "placenta", "testis": "testis", "bladder_organ": "urinary bladder",
    "synovial_membrane_of_synovial_joint": "blood", "tonsil": "tonsil",
    # Sub-regions map to the organ Census indexes them under. Left out, a third of the
    # cell-resolution corpus had no reference at all -- 33 carotid artery sections, 22 frontal
    # cortex, 20 lung parenchyma -- not because the reference is missing but because the atlas
    # records a finer anatomical term than Census's tissue_general.
    "frontal_cortex": "brain", "occipital_cortex": "brain", "parietal_cortex": "brain",
    "temporal_cortex": "brain", "primary_visual_cortex": "brain", "cerebellum": "brain",
    "ammon_s_horn": "brain", "striatum": "brain", "caudate_nucleus": "brain",
    "middle_temporal_gyrus": "brain", "inferior_parietal_cortex": "brain",
    "right_frontal_lobe": "brain", "hippocampal_formation": "brain",
    "lung_parenchyma": "lung", "pleura": "lung",
    "sigmoid_colon": "large intestine", "ascending_colon": "large intestine",
    "colonic_mucosa": "large intestine", "descending_colon": "large intestine",
    "transverse_colon": "large intestine",
    "ileum": "small intestine", "terminal_ileum": "small intestine",
    "duodenum": "small intestine", "jejunum": "small intestine", "intestine": "small intestine",
    "heart_left_ventricle": "heart", "heart_right_ventricle": "heart",
    "cardiac_atrium": "heart", "cardiac_ventricle": "heart",
    "carotid_artery_segment": "blood vessel", "ascending_aorta": "blood vessel",
    "coronary_artery": "blood vessel", "aorta": "blood vessel", "artery": "blood vessel",
    "urinary_bladder": "urinary bladder", "cortex_of_kidney": "kidney",
    "caudate_lobe_of_liver": "liver",
    "palatine_tonsil": "tonsil", "oropharynx": "tonsil",
    "skeletal_muscle_tissue": "skeletal muscle tissue",
    "subcutaneous_adipose_tissue": "adipose tissue", "adipose_tissue": "adipose tissue",
    "mucosa_of_vagina": "vagina", "vagina": "vagina",
    "sural_nerve": "nerve", "peripheral_nerve": "nerve",
}


def census_reference(organ: str, genes: list[str], census_version: str = "stable",
                     max_cells: int = MAX_REF_CELLS, min_per_type: int = MIN_CELLS_PER_TYPE,
                     seed: int = 0):
    """Per-organ reference restricted to `genes`, sampled before the matrix is fetched.

    Done in two passes on purpose. The obs pass is cheap and tells us which cells exist and what
    they are labelled; the matrix pass is the expensive one. Lung alone holds 3.1M Census cells
    across 253 types, so fetching the matrix first and subsampling after would move two orders of
    magnitude more data than the model will ever use. Sampling the indices first, then asking for
    exactly those rows, keeps each organ to a bounded fetch.

    The sample is balanced across cell types rather than drawn at random: Census is dominated by
    whatever assay happened to be run most, and an unbalanced reference teaches the model the
    prior rather than the biology.
    """
    import cellxgene_census

    rng = np.random.default_rng(seed)
    obs_filter = (f"tissue_general == '{organ}' and is_primary_data == True "
                  "and disease == 'normal'")
    with cellxgene_census.open_soma(census_version=census_version) as census:
        obs = cellxgene_census.get_obs(
            census, "Homo sapiens", value_filter=obs_filter,
            column_names=["soma_joinid", "cell_type", "assay", "donor_id"])
        if not len(obs):
            raise ValueError(f"no Census cells for organ={organ}")
        vc = obs.cell_type.value_counts()
        keep = vc[vc >= min_per_type].index
        obs = obs[obs.cell_type.isin(keep)]
        if not len(obs):
            raise ValueError(f"no Census type in {organ} reaches {min_per_type} cells")
        per = max(min_per_type, int(max_cells / max(1, obs.cell_type.nunique())))
        take = []
        for _, sub in obs.groupby("cell_type", observed=True):
            ids = sub.soma_joinid.to_numpy()
            take.append(rng.choice(ids, per, replace=False) if len(ids) > per else ids)
        coords = np.sort(np.concatenate(take))
        logger.info("%s: %d Census cells, %d types -> sampling %d",
                    organ, len(obs), obs.cell_type.nunique(), len(coords))

        a = cellxgene_census.get_anndata(
            census, organism="Homo sapiens",
            obs_coords=coords,
            var_value_filter=f"feature_name in {list(map(str, genes))}",
            column_names={"obs": ["cell_type", "assay", "donor_id"],
                          "var": ["feature_name", "feature_id"]},
        )
    a.var_names = a.var["feature_name"].astype(str)
    a.var_names_make_unique()
    a.obs["cell_type"] = a.obs["cell_type"].astype(str)
    return a


def _balance(a, rng, max_cells=MAX_REF_CELLS, min_per_type=MIN_CELLS_PER_TYPE):
    """Drop types too rare to learn, then cap the common ones so one type cannot dominate."""
    vc = a.obs.cell_type.value_counts()
    keep_types = vc[vc >= min_per_type].index
    a = a[a.obs.cell_type.isin(keep_types)].copy()
    if a.n_obs <= max_cells:
        return a
    per = max(min_per_type, int(max_cells / max(1, a.obs.cell_type.nunique())))
    idx = []
    for t, sub in a.obs.groupby("cell_type", observed=True):
        take = sub.index.to_numpy()
        if len(take) > per:
            take = rng.choice(take, per, replace=False)
        idx.extend(take)
    return a[a.obs_names.isin(set(idx))].copy()


def train(organ: str, genes: list[str], seed: int = 0, force: bool = False):
    """Train (or load) a CellTypist model for one organ restricted to one panel."""
    import celltypist
    import scanpy as sc

    MODELS.mkdir(parents=True, exist_ok=True)
    key = panel_key(genes)
    path = MODELS / f"{organ.replace(' ', '_')}__{key}.pkl"
    meta_path = path.with_suffix(".json")
    if path.exists() and not force:
        return celltypist.models.Model.load(str(path)), json.loads(meta_path.read_text())

    rng = np.random.default_rng(seed)
    ref = census_reference(organ, genes, seed=seed)
    ref = _balance(ref, rng)   # a second pass: types can still fall under the floor after the fetch
    present = [g for g in genes if g in set(ref.var_names)]
    if len(present) < MIN_GENES:
        raise ValueError(f"only {len(present)} of {len(genes)} panel genes in the {organ} "
                         "reference")
    ref = ref[:, present].copy()
    # CellTypist expects log1p of library-size-normalised counts at 10,000
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)
    model = celltypist.train(ref, labels="cell_type", n_jobs=8, feature_selection=False,
                             check_expression=False)
    model.write(str(path))
    meta = {"organ": organ, "panel_key": key, "n_genes_requested": len(genes),
            "n_genes_used": len(present), "n_ref_cells": int(ref.n_obs),
            "n_types": int(ref.obs.cell_type.nunique()),
            "types": sorted(map(str, ref.obs.cell_type.unique()))}
    meta_path.write_text(json.dumps(meta, indent=1))
    logger.info("trained %s: %d genes, %d cells, %d types", path.name, len(present),
                ref.n_obs, meta["n_types"])
    return model, meta


def annotate_one(path: str, organ: str, majority_voting: bool = True, seed: int = 0) -> dict:
    """Annotate one section in place: writes obs.cell_type / cell_type_conf / cell_type_model."""
    import anndata as ad
    import celltypist
    import scanpy as sc

    a = ad.read_h5ad(path)
    a.var_names_make_unique()
    genes = [str(g) for g in a.var_names]
    model, meta = train(organ, genes)

    use = [g for g in model.features if g in set(a.var_names)]
    sub = a[:, use].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    pred = celltypist.annotate(sub, model=model, majority_voting=majority_voting)
    lab = pred.predicted_labels
    col = "majority_voting" if majority_voting and "majority_voting" in lab.columns \
        else "predicted_labels"
    calls = lab[col].astype(str).to_numpy()
    conf = (pred.probability_matrix.max(axis=1).to_numpy()
            if hasattr(pred, "probability_matrix") else np.full(len(calls), np.nan))

    import h5py
    with h5py.File(path, "r+") as h:
        g = h.require_group("obs")
        for name, vals, dt in (("cell_type", calls, h5py.string_dtype("utf-8")),
                               ("cell_type_conf", conf.astype("float32"), "float32")):
            if name in g:
                del g[name]
            g.create_dataset(name, data=vals, dtype=dt)
        u = h.require_group("uns/hsa")
        for k, v in (("celltype_model", f"{organ}__{meta['panel_key']}"),
                     ("celltype_n_genes_used", len(use)),
                     ("celltype_n_types", meta["n_types"])):
            if k in u:
                del u[k]
            if isinstance(v, str):
                u.create_dataset(k, data=v, dtype=h5py.string_dtype("utf-8"))
            else:
                u.create_dataset(k, data=v)
    vc = pd.Series(calls).value_counts()
    return {"path": path, "organ": organ, "n_cells": int(a.n_obs), "n_genes_used": len(use),
            "n_types_called": int(vc.size), "top_type": str(vc.index[0]),
            "top_frac": float(vc.iloc[0] / len(calls)),
            "median_conf": float(np.nanmedian(conf)), "model": meta["panel_key"]}
