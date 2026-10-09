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
import hashlib
import re
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import db, tools_api as T
from .config import DATA, RUNS

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
    # 10x names the 2 um output directory square_002um, but GEO depositors use their own
    # conventions (slide2_2um_filtered_feature_bc_matrix.h5), so try progressively looser patterns
    # and fall back to every file for the sample -- build_queue() has already established that this
    # sample IS 2 um from its bin grid, so there is nothing to disambiguate when nothing matches.
    two = []
    for pat in (r"square_0*2um", r"0{2,}2um", r"[^0-9]2um", r"2um"):
        two = [p for p in paths if re.search(pat, str(p), re.I)]
        if two:
            break
    if not two:
        other = [p for p in paths if re.search(r"0*(8|16|32)um", str(p), re.I)]
        two = [p for p in paths if p not in other]
    if not two:
        return None
    mat = (next((p for p in two if p.name == "filtered_feature_bc_matrix.h5"), None)
           or next((p for p in two if "filtered_feature_bc_matrix" in p.name and p.suffix == ".h5"), None)
           or next((p for p in two if p.name.endswith(".h5") and "raw" not in p.name and "probe" not in p.name), None))
    # substring, not startswith: depositors prefix these (slide2_2um_tissue_positions.parquet.gz)
    pos = next((p for p in two if "tissue_positions" in p.name), None)
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
    # read_visium reads these by name and cannot handle gzip, so normalise the names and
    # decompress where needed rather than symlinking a .gz under an uncompressed name.
    def place(src: Path, dst: Path):
        if src.name.endswith(".gz"):
            import gzip as _gz
            with _gz.open(src, "rb") as fi, open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo, 1 << 24)
        else:
            link(src, dst)

    stem = pos.name[:-3] if pos.name.endswith(".gz") else pos.name
    pos_name = ("tissue_positions.parquet" if stem.endswith(".parquet")
                else "tissue_positions_list.csv" if "list" in stem
                else "tissue_positions.csv")
    place(pos, stage / "spatial" / pos_name)
    place(sf, stage / "spatial" / "scalefactors_json.json")
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


def _attach_positions(adata, stage: Path, img_path: Path | None) -> None:
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
    adata.uns["spatial"][lib]["metadata"] = {"source_image_path":
                                             None if img_path is None else str(img_path)}
    # Space Ranger ships downsampled hires/lowres overviews; HSA never downloaded them, so build
    # them from the source image using the scalefactors the depositor did provide. Skipped entirely
    # when there is no image: the gene-expression segmentation path works from array coordinates
    # and per-bin counts, so it never reads an overview.
    imgs = adata.uns["spatial"][lib].setdefault("images", {})
    if img_path is not None and "hires" not in imgs:
        import cv2
        src = cv2.imread(str(img_path), cv2.IMREAD_UNCHANGED)
        if src is None:
            import tifffile
            src = tifffile.imread(img_path)
        if src.ndim == 2:
            src = cv2.cvtColor(src, cv2.COLOR_GRAY2RGB)
        elif src.shape[-1] == 4:
            src = src[..., :3]
        for key, sk in (("hires", "tissue_hires_scalef"), ("lowres", "tissue_lowres_scalef")):
            f = float(sf.get(sk, 0.1))
            imgs[key] = cv2.resize(src, (max(1, int(src.shape[1] * f)), max(1, int(src.shape[0] * f))),
                                   interpolation=cv2.INTER_AREA)
        del src


def _tag(accession: str, sample_id: str) -> str:
    """Filesystem tag for one sample.

    The slug alone is not safe: 10x accessions plus sample ids run well past any truncation limit,
    and three distinct lung-cancer post-Xenium samples truncated to the SAME 120 characters, so two
    of them were skipped as "already_done" against the first one's output. The hash of the full,
    untruncated identity makes the name injective again.
    """
    slug = re.sub(r"[^A-Za-z0-9]+", "_", f"{accession}__{sample_id}")[:110].strip("_")
    h = hashlib.md5(f"{accession}__{sample_id}".encode()).hexdigest()[:8]
    return f"{slug}_{h}"


def _fetch_image(url: str, dest: Path, min_bps: float = 1.5e6, grace_s: float = 600.0,
                 max_s: float = 4 * 3600.0):
    """Stream a tissue image, giving up on a connection that trickles.

    requests' `timeout` with stream=True is a per-read deadline, not a total one, so a server that
    delivers a few kilobytes at a time keeps the transfer alive forever: one 9 GB image arrived at
    ~100 kB/s and held the batch for 17 hours without ever tripping the timeout. Enforce an actual
    throughput floor and a hard ceiling, and resume with a Range request on the next attempt.
    """
    start = dest.stat().st_size if dest.exists() else 0
    headers = {"Range": f"bytes={start}-"} if start else {}
    t0 = time.time()
    with T.S.get(url, stream=True, timeout=300, headers=headers) as r:
        if r.status_code not in (200, 206):
            return f"image_http_{r.status_code}", t0
        n = 0
        with open(dest, "ab" if start and r.status_code == 206 else "wb") as fh:
            for chunk in r.iter_content(1 << 23):
                fh.write(chunk)
                n += len(chunk)
                el = time.time() - t0
                if el > max_s:
                    return "image_download_exceeded_time_budget", t0
                if el > grace_s and n / el < min_bps:
                    logger.warning("image trickling at %.0f kB/s after %.0f s, abandoning (resumable)",
                                n / el / 1e3, el)
                    return "image_download_too_slow", t0
    return n, t0


def _stardist_params(image_path: Path) -> dict:
    """Tiling parameters that survive StarDist's own grid adjustment.

    stardist.big.cover() asserts `min_overlap + 2*context < block_size <= size`, and
    predict_instances_big first rounds each of those UP to a multiple of the model's grid (16).
    So passing block_size = size is not enough: 3350 was promoted to 3360 and then failed the
    assertion against a 3350 px image. Everything is therefore rounded DOWN to a multiple of 32,
    which is already grid-divisible and so passes through the adjustment unchanged.

    The gene-expression grid images make this necessary -- a 2 um sample rasterised at 2 um per
    pixel is 3350 px, below bin2cell's default 4096 block.
    """
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    with Image.open(image_path) as im:
        size = min(im.size)
    G = 32
    if size >= 4096 + 384:
        return {"block_size": 4096, "min_overlap": 128, "context": 128}
    block = max(G, (size // G) * G)
    context = max(G, ((block // 8) // G) * G)
    overlap = max(G, ((block // 16) // G) * G)
    while overlap + 2 * context >= block and context > G:
        context -= G
        overlap = max(G, overlap - G)
    if overlap + 2 * context >= block:              # image too small to tile at all
        raise ValueError(f"image {size}px too small for StarDist tiling")
    return {"block_size": block, "min_overlap": overlap, "context": context}


def _read_staged(stage: Path, img_path: Path | None):
    """AnnData from a staged Space Ranger layout, whether it holds an .h5 or an mtx trio.

    b2c.read_visium (like sc.read_visium) reads only filtered_feature_bc_matrix.h5, so the mtx
    fallback in _stage_spaceranger staged files that nothing downstream could open -- the sample
    failed with a missing-.h5 error after its image had already been downloaded. sc.read_10x_mtx
    reads the staged trio directly; spatial metadata is attached afterwards either way.
    """
    import bin2cell as b2c
    import scanpy as sc
    if (stage / "filtered_feature_bc_matrix.h5").exists():
        return b2c.read_visium(stage, source_image_path=img_path, load_images=False)
    d = stage / "filtered_feature_bc_matrix"
    if d.is_dir():
        return sc.read_10x_mtx(d)
    raise FileNotFoundError("staged layout has neither filtered_feature_bc_matrix.h5 nor an mtx dir")


def segment_one(source: str, accession: str, sample_id: str, image_url: str | None = None,
                mpp: float = 0.3, prob_thresh: float = 0.01, keep_image: bool = False,
                gex_mpp: float = 2.0, gex_prob_thresh: float = 0.05) -> dict:
    """Download image -> bin2cell -> cell-level .h5ad -> DELETE image. Peak extra disk: one image.

    With no image_url, segments on gene-expression density instead: b2c.grid_image() rasterises
    per-bin total counts and StarDist's fluorescence model finds nuclei in that. This is the only
    route for the 47 of 69 two-micron samples that have no public full-resolution H&E, so it is
    worth having even though H&E segmentation is the better of the two. Where an image does exist
    both are run and salvage_secondary_labels() fills H&E gaps with GEX calls, which is what the
    bin2cell authors' own workflow does.
    """
    import bin2cell as b2c

    rec = {"source": source, "accession": accession, "sample_id": sample_id, "status": "start"}
    HDOUT.mkdir(parents=True, exist_ok=True)
    IMG_TMP.mkdir(parents=True, exist_ok=True)
    tag = _tag(accession, sample_id)
    out = HDOUT / f"{tag}.h5ad"
    legacy_slug = re.sub(r"[^A-Za-z0-9]+", "_", f"{accession}__{sample_id}")[:120]
    legacy = HDOUT / f"{legacy_slug}.h5ad"
    if not out.exists() and legacy.exists():
        legacy.rename(out)                       # migrate off the collision-prone name
    if out.exists():
        return {**rec, "status": "already_done", "path": str(out)}

    stage = IMG_TMP / f"stage_{tag}"
    use_he = bool(image_url)
    img_path = (IMG_TMP / f"{tag}__{Path(image_url.split('?')[0]).name}") if use_he else None
    scaled = IMG_TMP / f"{tag}_scaled.tiff"
    labels = IMG_TMP / f"{tag}_labels.npz"
    gex_img = IMG_TMP / f"{tag}_gex.tiff"
    gex_labels = IMG_TMP / f"{tag}_gex_labels.npz"
    try:
        if _stage_spaceranger(source, accession, sample_id, stage) is None:
            return {**rec, "status": "cannot_stage_spaceranger_layout"}

        if use_he and not img_path.exists():
            part = img_path.with_suffix(img_path.suffix + ".part")
            got, t0 = _fetch_image(image_url, part)
            if isinstance(got, str):
                # A slow or missing image is no longer fatal: fall back to segmenting on gene
                # expression. Returning here forfeited the sample entirely, which is the wrong
                # trade when the GEX route is available and costs no download at all.
                logger.warning("%s: %s -> falling back to GEX-only segmentation", sample_id, got)
                rec["image_status"] = got
                use_he, img_path = False, None
            else:
                part.rename(img_path)
        if use_he:
            rec["image_gb"] = round(img_path.stat().st_size / 1e9, 2)
            img_path = _prepare_image(img_path)   # gunzip / BigTIFF -> something cv2 can read

        # load_images=False is required (the SR hires/lowres pngs were never downloaded), but that
        # code path also skips the tissue-position join, so attach positions + scalefactors here.
        adata = _read_staged(stage, img_path)
        adata.var_names_make_unique()
        _attach_positions(adata, stage, img_path)
        adata = adata[:, adata.X.sum(0).A1 > 0].copy() if hasattr(adata.X, "A1") else adata
        # Drop empty bins BEFORE destriping: at 2 um most bins are empty, and destripe divides by a
        # per-row/column quantile that is then 0, producing inf that only surfaces much later.
        adata.obs["n_counts"] = np.asarray(adata.X.sum(1)).ravel()
        adata = adata[adata.obs["n_counts"] > 0].copy()
        b2c.destripe(adata)
        for key in ("destripe_factor", "n_counts_adjusted"):
            if key in adata.obs:
                v = adata.obs[key].to_numpy(dtype=float)
                adata.obs[key] = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)                                   # remove the HD row/column striping
        primary = None
        if use_he:
            b2c.scaled_he_image(adata, mpp=mpp, save_path=str(scaled))
            b2c.stardist(image_path=str(scaled), labels_npz_path=str(labels),
                         stardist_model="2D_versatile_he", prob_thresh=prob_thresh,
                         **_stardist_params(scaled))
            b2c.insert_labels(adata, labels_npz_path=str(labels), basis="spatial",
                              spatial_key="spatial_cropped_150_buffer", mpp=mpp,
                              labels_key="labels_he")
            primary = "labels_he"

        # gene-expression nuclei: rasterise per-bin counts, then StarDist's fluorescence model
        b2c.grid_image(adata, "n_counts_adjusted" if "n_counts_adjusted" in adata.obs
                       else "n_counts", mpp=gex_mpp, sigma=5, save_path=str(gex_img))
        b2c.stardist(image_path=str(gex_img), labels_npz_path=str(gex_labels),
                     stardist_model="2D_versatile_fluo", prob_thresh=gex_prob_thresh,
                     **_stardist_params(gex_img))
        b2c.insert_labels(adata, labels_npz_path=str(gex_labels), basis="array",
                          mpp=gex_mpp, labels_key="labels_gex")

        # salvage_secondary_labels() expects the primary labels already EXPANDED -- its own default
        # is labels_he_expanded -- so expand first and salvage into the final key, rather than
        # salvaging raw labels and expanding the union afterwards.
        if primary is None:
            b2c.expand_labels(adata, labels_key="labels_gex",
                              expanded_labels_key="labels_expanded")
            method = f"stardist 2D_versatile_fluo mpp={gex_mpp} p={gex_prob_thresh} (GEX only)"
            label_key = "labels_gex"
        else:
            b2c.expand_labels(adata, labels_key="labels_he",
                              expanded_labels_key="labels_he_expanded")
            b2c.salvage_secondary_labels(adata, primary_label="labels_he_expanded",
                                         secondary_label="labels_gex",
                                         labels_key="labels_expanded")
            label_key = "labels_joint"
            method = (f"stardist 2D_versatile_he mpp={mpp} p={prob_thresh} + GEX salvage "
                      f"(2D_versatile_fluo mpp={gex_mpp})")
        spatial_keys = ["spatial", "spatial_cropped_150_buffer"] if use_he else ["spatial"]
        cells = b2c.bin_to_cell(adata, labels_key="labels_expanded", spatial_keys=spatial_keys)
        cells.uns["hsa"] = {"source": source, "accession": accession, "sample_id": sample_id,
                            "technology": "Visium HD", "unit": "cell",
                            "cell_calling": f"bin2cell {method}",
                            "image": Path(image_url).name if use_he else None,
                            "image_url": image_url}
        cells.write_h5ad(out, compression="gzip")
        rec.update(status="ok", n_cells=int(cells.n_obs), n_genes=int(cells.n_vars),
                   n_bins=int(adata.n_obs), path=str(out), label_source=label_key)
    except Exception as e:
        import traceback
        rec["status"] = f"error: {type(e).__name__}: {e}"[:220]
        rec["traceback"] = traceback.format_exc()[-1200:]
        logger.error("%s / %s failed:\n%s", accession, sample_id, rec["traceback"])
    finally:
        if not keep_image:
            for f in list(IMG_TMP.glob(f"{tag}*")) + [scaled, labels, gex_img, gex_labels]:
                if f.is_file():
                    f.unlink()
            if stage.exists():
                shutil.rmtree(stage, ignore_errors=True)
    return rec


def build_queue() -> pd.DataFrame:
    """Every genuine Visium HD sample at 2 um bins, with its image if one exists.

    Supersedes the image-first eligibility, which required a full-resolution H&E on the same row
    and so admitted only 22 of the 69 two-micron samples. Since segment_one() can segment on gene
    expression alone, having no image is no longer disqualifying -- it only decides which label
    source is used. 8 and 16 um samples are deliberately excluded: an 8 um bin is already roughly
    cell-sized, so "cell calling" on it would not mean anything.
    """
    res = pd.read_csv(RUNS / "technology_resolved.csv")
    key = ["source", "accession", "sample_id"]
    res[key] = res[key].astype(str)
    two = res[(res.resolved == "Visium HD") & (res.bin_um == 2.0)].copy()
    surv = RUNS / "hd_image_survey_all.csv"
    if surv.exists():
        img = pd.read_csv(surv)
        img[key] = img[key].astype(str)
        two = two.merge(img[key + ["image_gb", "image_url"]], on=key, how="left")
    else:
        two["image_gb"], two["image_url"] = float("nan"), None
    two["label_source"] = np.where(two.image_url.notna(), "he+gex", "gex")
    out = two[key + ["bin_um", "n_cells", "image_gb", "image_url", "label_source"]]
    out.to_csv(RUNS / "hd_b2c_queue.csv", index=False)
    return out


def segment_all(csv: str = "runs/hd_b2c_queue.csv", limit: int | None = None,
                require_image: bool = False) -> pd.DataFrame:
    """Stream every queued sample: one image on disk at a time."""
    if not Path(csv).exists():
        build_queue()
    todo = pd.read_csv(csv)
    if require_image:
        todo = todo[todo.image_gb > 0]
    todo = todo.sort_values("image_url", na_position="last")   # H&E samples first
    if limit:
        todo = todo.head(limit)
    res, outcsv = [], Path("runs/bin2cell_results.csv")
    for i, (_, r) in enumerate(todo.iterrows(), 1):
        has_img = isinstance(r.get("image_url"), str) and r["image_url"].strip()
        logger.info("[%d/%d] %s / %s (%s)", i, len(todo), r.accession, r.sample_id,
                    f"{r.image_gb:.1f} GB H&E" if has_img else "no image -> GEX segmentation")
        rec = segment_one(r.source, r.accession, str(r.sample_id),
                          r["image_url"] if has_img else None)
        logger.info("    -> %s", rec.get("status"))
        res.append(rec)
        pd.DataFrame(res).to_csv(outcsv, index=False)
    return pd.DataFrame(res)
