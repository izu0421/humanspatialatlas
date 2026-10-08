"""Cell-level Visium HD by nuclei segmentation of the H&E image (bin2cell).

Space Ranger 4's own segmentation needs FASTQs *and* an H&E image, neither of which HSA holds, so
re-running it is not possible. bin2cell reaches the same place from what we do have: 2 um bins plus
the full-resolution image.

Images are 5-6 GB each and /data has ~2 TB free, so this streams one sample at a time:
    locate image -> download -> segment -> write cell-level .h5ad -> DELETE image
Peak extra disk is one image.
"""
from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from . import db, tools_api as T
from .config import DATA

logger = logging.getLogger(__name__)
HDOUT = DATA.parent / "visium_hd_cells"
IMG_TMP = DATA.parent / "tmp_images"
# full-resolution H&E only: downsampled tif/png cannot resolve nuclei
# any TIFF-family image; SIZE is the real discriminator, since names vary wildly between submitters
FULLRES = re.compile(r"\.(ome\.tiff?|tiff?|btf|qptiff|svs|ndpi)(\.gz)?$", re.I)
# exclude CytAssist images (~10-30 MB) and fiducial/thumbnail jpgs: those cannot resolve nuclei
TOO_SMALL = 200e6


def eligible_samples() -> pd.DataFrame:
    """HD samples that have 2 um bins + scalefactors, i.e. everything bin2cell needs but the image."""
    from .export import norm_tech
    d = pd.DataFrame(db.query("SELECT source, accession, technology FROM datasets "
                              "WHERE verdict IN ('ready','bundle_only')"))
    d["t"] = d.technology.map(norm_tech)
    hd = d[d.t == "Visium HD"][["source", "accession"]]
    f = pd.DataFrame(db.query("SELECT source, accession, sample_id, local_path FROM files "
                              "WHERE dl_status='done' AND local_path IS NOT NULL")).merge(hd, on=["source", "accession"])
    f["two"] = f.local_path.str.contains("square_002um")
    f["sf"] = f.local_path.str.contains("scalefactors")
    g = f.groupby(["source", "accession", "sample_id"]).agg(two=("two", "any"), sf=("sf", "any")).reset_index()
    return g[g.two & g.sf].drop(columns=["two", "sf"])


def _probe(url: str) -> dict | None:
    """HEAD a candidate image; keep it only if it is real and big enough to resolve nuclei."""
    try:
        h = T.S.head(url, timeout=60, allow_redirects=True)
    except Exception:
        return None
    if h.status_code != 200:
        return None
    size = int(h.headers.get("Content-Length", 0) or 0)
    if size < TOO_SMALL:
        return None
    return {"name": url.rsplit("/", 1)[-1], "url": url, "size_bytes": size}


def _tenx_candidates(accession: str) -> list[str]:
    """10x publishes the full-res H&E next to the other outputs but does not link it on the page.

    Derive `<dir>/<Name>_tissue_image.btf` from any file URL we already recorded for this dataset.
    The S3 origin is used because the CDN stalls on large range reads.
    """
    rows = db.query("SELECT url FROM files WHERE accession=? AND url LIKE 'http%' LIMIT 40", (accession,))
    out = []
    for r in rows:
        u = r["url"].split("#")[0].split("?")[0]
        if "10xgenomics.com" not in u and "10x.files" not in u:
            continue
        u = u.replace("https://cf.10xgenomics.com/", "https://s3-us-west-2.amazonaws.com/10x.files/")
        d, base = u.rsplit("/", 1)
        name = d.rsplit("/", 1)[-1]
        for suf in ("_tissue_image.btf", "_tissue_image.ome.tif", "_image.btf"):
            cand = f"{d}/{name}{suf}"
            if cand not in out:
                out.append(cand)
    return out


def find_image_url(source: str, accession: str, sample_id: str) -> dict | None:
    """Locate a full-resolution image for one sample, per source. Returns None if none is public."""
    try:
        if source == "10x Datasets":
            for c in _tenx_candidates(accession):
                if (hit := _probe(c)):
                    return hit
            return None

        if source == "GEO":
            m = re.search(r"(GSM\d+)", sample_id)
            files = []
            if m:
                files += T.geo_all_sample_files([m.group(1)]).get(m.group(1), [])
            if not files or not any(f.get("role") == "image" for f in files):
                files += T.geo_series(accession).get("series_suppl", [])     # some submit at series level
        elif source == "Zenodo":
            files = T.zenodo_files(accession)["files"]
        elif source == "figshare":
            files = T.figshare_files(accession)["files"]
        elif source == "Mendeley Data":
            files = T.mendeley_files(accession)["files"]
        else:
            return None
    except Exception as e:
        logger.warning("image lookup failed %s/%s: %s", accession, sample_id, str(e)[:80])
        return None

    cands = [f for f in files if FULLRES.search(f.get("name", ""))
             and (f.get("size_bytes") or 0) > TOO_SMALL]
    # prefer an image whose name shares the sample's own token, else the largest
    tok = re.sub(r"[^A-Za-z0-9]", "", str(sample_id).split("/")[0]).lower()
    named = [f for f in cands if tok and tok in re.sub(r"[^A-Za-z0-9]", "", f["name"]).lower()]
    pool = named or cands
    return max(pool, key=lambda f: f.get("size_bytes") or 0) if pool else None


def all_hd_samples() -> pd.DataFrame:
    from .export import norm_tech
    d = pd.DataFrame(db.query("SELECT source, accession, technology FROM datasets "
                              "WHERE verdict IN ('ready','bundle_only')"))
    d["t"] = d.technology.map(norm_tech)
    hd = d[d.t == "Visium HD"][["source", "accession"]]
    f = pd.DataFrame(db.query("SELECT DISTINCT source, accession, sample_id FROM files "
                              "WHERE dl_status='done' AND local_path IS NOT NULL"))
    return f.merge(hd, on=["source", "accession"])


def survey(limit: int | None = None, everything: bool = True) -> pd.DataFrame:
    """What fraction of HD samples have a reachable full-res image, and how big?"""
    el = all_hd_samples() if everything else eligible_samples()
    if limit:
        el = el.head(limit)
    rows = []
    for _, r in el.iterrows():
        img = find_image_url(r.source, r.accession, r.sample_id)
        rows.append({**r.to_dict(), "image": (img or {}).get("name"),
                     "image_gb": round((img or {}).get("size_bytes", 0) / 1e9, 2) if img else 0.0,
                     "image_url": (img or {}).get("url")})
    return pd.DataFrame(rows)


def _stage_spaceranger(source: str, accession: str, sample_id: str, stage: Path) -> Path | None:
    """Symlink this sample's scattered files into the Space Ranger layout bin2cell expects.

    HSA stores files wherever the depositor put them; b2c.read_visium wants
        <dir>/filtered_feature_bc_matrix.h5
        <dir>/spatial/{tissue_positions*, scalefactors_json.json}
    Symlinks cost no disk and leave the originals untouched.
    """
    paths = [Path(r["local_path"]) for r in db.query(
        "SELECT local_path FROM files WHERE source=? AND accession=? AND sample_id=? "
        "AND dl_status='done' AND local_path IS NOT NULL", (source, accession, sample_id))
        if Path(r["local_path"]).exists()]
    two = [p for p in paths if "square_002um" in str(p)]
    if not two:
        return None
    mat = (next((p for p in two if p.name == "filtered_feature_bc_matrix.h5"), None)
           or next((p for p in two if "filtered_feature_bc_matrix" in p.name and p.suffix == ".h5"), None)
           or next((p for p in two if p.name.endswith(".h5") and "raw" not in p.name and "probe" not in p.name), None))
    pos = next((p for p in two if p.name.startswith("tissue_positions")), None)
    sf = next((p for p in two if "scalefactors" in p.name), None)
    if not (pos and sf):
        return None

    stage.mkdir(parents=True, exist_ok=True)
    (stage / "spatial").mkdir(exist_ok=True)
    def link(src: Path, dst: Path):
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src.resolve())

    if mat:
        link(mat, stage / "filtered_feature_bc_matrix.h5")
    else:                                        # mtx trio instead of an .h5
        d = stage / "filtered_feature_bc_matrix"
        d.mkdir(exist_ok=True)
        for pat, name in [("matrix.mtx", "matrix.mtx.gz"), ("barcodes.tsv", "barcodes.tsv.gz"),
                          ("features.tsv", "features.tsv.gz")]:
            src = next((p for p in two if pat in p.name and "raw" not in str(p)), None)
            if src is None:
                return None
            link(src, d / name)
    link(pos, stage / "spatial" / pos.name)
    link(sf, stage / "spatial" / "scalefactors_json.json")
    for extra in two:                            # hires/lowres images if the depositor included them
        if "hires_image" in extra.name or "lowres_image" in extra.name:
            link(extra, stage / "spatial" / extra.name)
    return stage


def _prepare_image(img: Path) -> Path:
    """Make the image readable by bin2cell (which goes through cv2.imread).

    GEO deposits gzipped TIFFs and 10x ships BigTIFF (.btf); cv2 reads neither, and returns an empty
    array rather than raising, which surfaces much later as an opaque assertion inside cvtColor.
    """
    import gzip as _gzip
    if img.suffix == ".gz":
        out = img.with_suffix("")
        if not out.exists():
            with _gzip.open(img, "rb") as fi, open(out, "wb") as fo:
                shutil.copyfileobj(fi, fo, 1 << 24)
        img.unlink()
        img = out
    import cv2
    if cv2.imread(str(img), cv2.IMREAD_UNCHANGED) is None:   # BigTIFF and friends
        import tifffile
        arr = tifffile.imread(img)
        out = img.with_suffix(".png") if max(arr.shape[:2]) < 60000 else img.with_suffix(".tiff")
        tifffile.imwrite(out, arr, bigtiff=False) if out.suffix == ".tiff" else cv2.imwrite(str(out), arr)
        img.unlink()
        img = out
    return img


def _attach_positions(adata, stage: Path, img_path: Path) -> None:
    """Join tissue positions and scalefactors onto an AnnData read with load_images=False."""
    import json as _json
    pos = next(iter((stage / "spatial").glob("tissue_positions*")), None)
    if pos is None:
        raise FileNotFoundError("no tissue_positions in staged spatial/")
    df = pd.read_parquet(pos) if pos.suffix == ".parquet" else pd.read_csv(
        pos, header=0 if pos.name == "tissue_positions.csv" else None,
        names=None if pos.name == "tissue_positions.csv" else
        ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    df = df.set_index("barcode")
    df.index = df.index.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    common = adata.obs_names.intersection(df.index)
    if len(common) < 0.5 * adata.n_obs:
        raise ValueError(f"positions match only {len(common)}/{adata.n_obs} barcodes")
    adata._inplace_subset_obs(adata.obs_names.isin(common))
    for c in ["in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]:
        adata.obs[c] = df.loc[adata.obs_names, c].to_numpy()
    adata.obsm["spatial"] = adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
    sf_path = stage / "spatial" / "scalefactors_json.json"
    sf = _json.loads(sf_path.read_text()) if sf_path.exists() else {}
    lib = next(iter(adata.uns.get("spatial", {})), "library")
    adata.uns.setdefault("spatial", {}).setdefault(lib, {})
    adata.uns["spatial"][lib]["scalefactors"] = sf
    adata.uns["spatial"][lib].setdefault("images", {})
    adata.uns["spatial"][lib]["metadata"] = {"source_image_path": str(img_path)}


def segment_one(source: str, accession: str, sample_id: str, image_url: str,
                mpp: float = 0.3, prob_thresh: float = 0.01, keep_image: bool = False) -> dict:
    """Download image -> bin2cell -> cell-level .h5ad -> DELETE image. Peak extra disk: one image."""
    import bin2cell as b2c

    rec = {"source": source, "accession": accession, "sample_id": sample_id, "status": "start"}
    HDOUT.mkdir(parents=True, exist_ok=True)
    IMG_TMP.mkdir(parents=True, exist_ok=True)
    tag = re.sub(r"[^A-Za-z0-9]+", "_", f"{accession}__{sample_id}")[:120]
    out = HDOUT / f"{tag}.h5ad"
    if out.exists():
        return {**rec, "status": "already_done", "path": str(out)}

    stage = IMG_TMP / f"stage_{tag}"
    img_path = IMG_TMP / f"{tag}__{Path(image_url.split('?')[0]).name}"
    scaled = IMG_TMP / f"{tag}_scaled.tiff"
    labels = IMG_TMP / f"{tag}_labels.npz"
    try:
        if _stage_spaceranger(source, accession, sample_id, stage) is None:
            return {**rec, "status": "cannot_stage_spaceranger_layout"}

        if not img_path.exists():
            with T.S.get(image_url, stream=True, timeout=900) as r:
                if r.status_code != 200:
                    return {**rec, "status": f"image_http_{r.status_code}"}
                with open(img_path, "wb") as fh:
                    for chunk in r.iter_content(1 << 23):
                        fh.write(chunk)
        rec["image_gb"] = round(img_path.stat().st_size / 1e9, 2)
        img_path = _prepare_image(img_path)       # gunzip / BigTIFF -> something cv2 can read

        # load_images=False is required (the SR hires/lowres pngs were never downloaded), but that
        # code path also skips the tissue-position join, so attach positions + scalefactors here.
        adata = b2c.read_visium(stage, source_image_path=img_path, load_images=False)
        adata.var_names_make_unique()
        _attach_positions(adata, stage, img_path)
        adata = adata[:, adata.X.sum(0).A1 > 0].copy() if hasattr(adata.X, "A1") else adata
        adata.obs["n_counts"] = np.asarray(adata.X.sum(1)).ravel()   # destripe expects this
        b2c.destripe(adata)                                   # remove the HD row/column striping
        b2c.scaled_he_image(adata, mpp=mpp, save_path=str(scaled))
        b2c.stardist(image_path=str(scaled), labels_npz_path=str(labels),
                     stardist_model="2D_versatile_he", prob_thresh=prob_thresh)
        b2c.insert_labels(adata, labels_npz_path=str(labels), basis="spatial",
                          spatial_key="spatial_cropped_150_buffer", mpp=mpp, labels_key="labels_he")
        b2c.expand_labels(adata, labels_key="labels_he", expanded_labels_key="labels_he_expanded")
        cells = b2c.bin_to_cell(adata, labels_key="labels_he_expanded",
                                spatial_keys=["spatial", "spatial_cropped_150_buffer"])
        cells.uns["hsa"] = {"source": source, "accession": accession, "sample_id": sample_id,
                            "technology": "Visium HD", "unit": "cell",
                            "cell_calling": f"bin2cell stardist 2D_versatile_he mpp={mpp} p={prob_thresh}",
                            "image": Path(image_url).name, "image_url": image_url}
        cells.write_h5ad(out, compression="gzip")
        rec.update(status="ok", n_cells=int(cells.n_obs), n_genes=int(cells.n_vars),
                   n_bins=int(adata.n_obs), path=str(out))
    except Exception as e:
        rec["status"] = f"error: {type(e).__name__}: {e}"[:220]
    finally:
        if not keep_image:
            for f in list(IMG_TMP.glob(f"{tag}*")) + [scaled, labels]:
                if f.is_file():
                    f.unlink()
            if stage.exists():
                shutil.rmtree(stage, ignore_errors=True)
    return rec


def segment_all(csv: str = "runs/hd_ready_for_b2c.csv", limit: int | None = None) -> pd.DataFrame:
    """Stream every ready sample: one image on disk at a time."""
    todo = pd.read_csv(csv)
    todo = todo[todo.image_gb > 0]
    if limit:
        todo = todo.head(limit)
    res, outcsv = [], Path("runs/bin2cell_results.csv")
    for i, (_, r) in enumerate(todo.iterrows(), 1):
        logger.info("[%d/%d] %s / %s (%.1f GB image)", i, len(todo), r.accession, r.sample_id, r.image_gb)
        rec = segment_one(r.source, r.accession, str(r.sample_id), r.image_url)
        logger.info("    -> %s", rec.get("status"))
        res.append(rec)
        pd.DataFrame(res).to_csv(outcsv, index=False)
    return pd.DataFrame(res)
