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
REFS = RUNS.parent / "celltype_refs"
OUT = RUNS / "annotation_log.csv"

CELL_TECH = ("Xenium", "CosMx", "MERFISH", "Atera", "cell-level WTA", "STARmap", "seqFISH")
MIN_GENES = 50            # below this a panel cannot support a multi-class model
MIN_REF_CELLS = 500       # per organ
# 40,000 balanced cells, not 120,000. The model has at most a few hundred features (a panel),
# so the extra rows buy very little and cost a lot: the 97,816-cell lung fetch died part way
# through pulling its matrix. A bounded fetch that completes is worth more than a large one
# that does not.
MAX_REF_CELLS = 25_000    # per organ, balanced across types; the fetch is the slow step
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
    "synovial_membrane_of_synovial_joint": "blood",
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
    # Census calls this "vasculature"; "blood vessel" returns nothing and the 47 carotid, aorta
    # and coronary sections would have had no reference at all.
    "carotid_artery_segment": "vasculature", "ascending_aorta": "vasculature",
    "coronary_artery": "vasculature", "aorta": "vasculature", "artery": "vasculature",
    "urinary_bladder": "urinary bladder", "cortex_of_kidney": "kidney",
    "caudate_lobe_of_liver": "liver",
    # Census names these differently from the atlas, and all three returned nothing under the
    # obvious term: "musculature" not "skeletal muscle tissue", "reproductive system" not
    # "vagina", and tonsil has no term of its own -- it is secondary lymphoid tissue, so the
    # lymph node reference is the closest honest match and is noted as a substitution.
    "skeletal_muscle_tissue": "musculature", "muscle_tissue": "musculature",
    "mucosa_of_vagina": "reproductive system", "vagina": "reproductive system",
    "tonsil": "lymph node", "palatine_tonsil": "lymph node", "oropharynx": "lymph node",
    "subcutaneous_adipose_tissue": "adipose tissue", "adipose_tissue": "adipose tissue",
    # Peripheral nerve has no tissue_general in Census (only "central nervous system", which is
    # a different tissue). 8 sections are left without a reference rather than annotated from a
    # CNS model that has never seen a Schwann cell.
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

        # Genes are selected by integer coordinate, not by a `feature_name in [...]` filter.
        # A 480-term filter string is handed to the tiledbsoma C++ layer and the fetch died
        # there with no Python traceback -- the silent-exit signature of a crash below the
        # interpreter, not an exception. Resolving the names to indices once keeps the query
        # trivial and makes the gene mapping inspectable.
        var = cellxgene_census.get_var(census, "Homo sapiens",
                                       column_names=["soma_joinid", "feature_name"])
        want = set(map(str, genes))
        hit = var[var.feature_name.astype(str).isin(want)]
        if not len(hit):
            raise ValueError(f"none of the {len(want)} panel genes are in Census")
        var_coords = np.sort(hit.soma_joinid.to_numpy())
        logger.info("%s: %d of %d panel genes resolved in Census", organ, len(hit), len(want))

        # Fetched in blocks of cells: one request for tens of thousands of rows is the step that
        # kept failing, and a failure part way through a single request leaves nothing.
        import anndata as ad_
        blocks, B = [], 5_000
        for i in range(0, len(coords), B):
            blocks.append(cellxgene_census.get_anndata(
                census, organism="Homo sapiens",
                obs_coords=coords[i:i + B], var_coords=var_coords,
                column_names={"obs": ["cell_type", "assay", "donor_id"],
                              "var": ["feature_name", "feature_id"]}))
            logger.info("  fetched %d/%d cells", min(i + B, len(coords)), len(coords))
        a = ad_.concat(blocks, join="outer", merge="first") if len(blocks) > 1 else blocks[0]
        if "feature_name" not in a.var.columns:
            a.var = blocks[0].var.reindex(a.var_names)
    a.var_names = a.var["feature_name"].astype(str)
    # anndata refuses to write a var index whose name collides with a column of different values,
    # so the column goes once it has been promoted to the index. The fetch itself had succeeded;
    # only the write failed, which is the worst place to lose a 20-minute network call.
    a.var = a.var.drop(columns=["feature_name"])
    a.var.index.name = None
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


def cached_reference(organ: str, genes: list[str], seed: int = 0, force: bool = False):
    """Fetch an organ's reference once and reuse it for every panel in that organ.

    A Census fetch is the slow step by a wide margin -- block times for one lung reference ranged
    from 24 seconds to 13 minutes, all network -- and an organ carries several panels that mostly
    overlap. Fetching per (organ, panel) pays that cost once per panel; fetching the UNION of the
    organ's panel genes once and subsetting locally pays it once per organ.

    The cache is keyed on the organ and checked for coverage, so a later panel with genes the
    cached reference does not hold triggers one re-fetch of the widened union rather than
    silently training on a subset of its own features.
    """
    import anndata as ad

    REFS.mkdir(parents=True, exist_ok=True)
    f = REFS / f"{organ.replace(' ', '_')}.h5ad"
    want = {str(g) for g in genes}
    if f.exists() and not force:
        a = ad.read_h5ad(f)
        # A cache file that opens is not a cache file that is usable. A write that failed after
        # the obs axis left 20,000-cell references with ZERO genes on disk -- readable, the right
        # n_obs, and worthless. Anything downstream would have trained on nothing.
        if a.n_vars == 0 or a.n_obs == 0:
            logger.warning("%s: cached reference is empty %s, refetching", organ, a.shape)
            f.unlink(missing_ok=True)
        else:
            # Census does not carry every panel gene, so a gene missing from the cache is only a
            # reason to refetch when the cache was built from a NARROWER request than this one --
            # which is why the request itself is recorded at build time.
            built_from = {str(x) for x in a.uns.get("requested_genes", [])}
            if not (want - built_from):
                return a
            logger.info("%s: cache built from %d genes, now want %d; widening", organ,
                        len(built_from), len(want))
            want |= built_from
    a = census_reference(organ, sorted(want), seed=seed)
    a.uns["requested_genes"] = sorted(want)
    if a.n_vars == 0 or a.n_obs == 0:
        raise ValueError(f"{organ}: fetched reference is empty {a.shape}; refusing to cache")
    a.write_h5ad(f, compression="gzip")
    check = ad.read_h5ad(f, backed="r")          # read back: the write is the step that failed
    if check.n_vars != a.n_vars or check.n_obs != a.n_obs:
        f.unlink(missing_ok=True)
        raise ValueError(f"{organ}: wrote {a.shape} but read back {check.shape}")
    logger.info("%s: cached reference %s -> %s", organ, a.shape, f.name)
    return a


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
    ref = cached_reference(organ, genes, seed=seed)
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
