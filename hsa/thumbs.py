"""Spot-scatter thumbnails: a picture of each section, from its own coordinates.

A catalogue of 4,535 datasets described only in bars and counts gives no sense of
what is in it. Every standardised sample already carries the coordinates it was
measured on, so the cheapest honest illustration is the measurement itself: the
spots or cells plotted where they actually sit, coloured by cell type where the
section has one.

WHY THESE ARE EMBEDDED RATHER THAN SERVED AS FILES. The dashboard is a single
AES-encrypted page behind a password. A thumbnail written as a file next to it
would be fetched without that password, so the shape of every section in the
atlas would be public while the catalogue describing it was not. They are
therefore data URIs inside the encrypted payload, which bounds how many there can
be: about 1.4 KB each after base64, so a few hundred, not 7,705.

WHAT IS DRAWN, AND WHAT IS NOT. Points are subsampled to at most 8,000, so a
625k-cell Xenium section and a 2k-spot Visium section are drawn with comparable
ink; a thumbnail says where the measurements are, not how many. Aspect ratio is
preserved and the scale is per section, so these are shapes, not a size
comparison. Sections whose coordinates are degenerate -- a single point, a line,
all NaN -- are skipped rather than drawn as a dot in the middle, because a
thumbnail that looks like a failed section and one that IS a failed section have
to be distinguishable.

Usage:
    micromamba run -n hsa python -I -m hsa.thumbs --limit 240
"""
from __future__ import annotations

import argparse
import base64
import collections
import io
import json
import logging
import time
from pathlib import Path

import numpy as np

log = logging.getLogger("thumbs")
ROOT = Path(__file__).resolve().parent.parent
STD = ROOT / "standardised"
CACHE = ROOT / "exports_private" / "thumbs.json"

# the project palette; the first colours go to the commonest cell types
BLUE, RED, GREY = (33, 102, 172), (178, 24, 43), (154, 154, 162)
PAL = [BLUE, RED, (26, 122, 106), (201, 135, 40), (120, 94, 170), (90, 90, 98)]
SIZE = 96
MAX_PTS = 8000


def render(xy, labels=None, size: int = SIZE, pad: int = 3) -> bytes | None:
    """One section as a transparent PNG, or None if its coordinates say nothing."""
    from PIL import Image

    xy = np.asarray(xy, dtype=float)
    if xy.ndim != 2 or xy.shape[1] < 2 or len(xy) < 8:
        return None
    xy = xy[:, :2]
    ok = np.isfinite(xy).all(axis=1)
    xy = xy[ok]
    if len(xy) < 8:
        return None
    if labels is not None:
        labels = np.asarray(labels)[ok]
    if len(xy) > MAX_PTS:
        pick = np.random.default_rng(0).choice(len(xy), MAX_PTS, replace=False)
        xy = xy[pick]
        if labels is not None:
            labels = labels[pick]

    lo, hi = xy.min(axis=0), xy.max(axis=0)
    extent = hi - lo
    span = float(np.max(extent))
    # a section collapsed onto a point or a line is not a shape
    if not np.isfinite(span) or span <= 0 or float(np.min(extent)) <= span * 1e-4:
        return None

    # Drawn at SS times the final size and then resampled down. A capture array is a
    # regular lattice, and plotting one point per pixel onto a 96px raster beats that
    # lattice against the pixel grid: the first version striped every Visium section
    # with horizontal moire bands that are a sampling artefact and not the tissue.
    # Supersampling also makes density legible, because a pixel that several points
    # land in comes out darker instead of simply being overwritten.
    SS = 3
    big = size * SS
    q = ((xy - (lo + hi) / 2) / span * (big - 2 * pad * SS) + big / 2).astype(int)
    q[:, 1] = big - 1 - q[:, 1]                        # image y grows downward
    img = np.zeros((big, big, 4), np.uint8)
    if labels is None:
        cols = np.tile(np.array(BLUE + (255,), np.uint8), (len(q), 1))
    else:
        top = [t for t, _ in collections.Counter(labels).most_common(len(PAL))]
        idx = {t: i for i, t in enumerate(top)}
        cols = np.array([PAL[idx[t]] + (255,) if t in idx else GREY + (255,)
                         for t in labels], np.uint8)
    yy = np.clip(q[:, 1], 0, big - 1)
    xx = np.clip(q[:, 0], 0, big - 1)
    # a point is a small disc at this scale, so the lattice reads as spots
    for dy in range(SS - 1):
        for dx in range(SS - 1):
            img[np.clip(yy + dy, 0, big - 1), np.clip(xx + dx, 0, big - 1)] = cols

    im = Image.fromarray(img, "RGBA").resize((size, size), Image.LANCZOS)
    im = im.quantize(colors=48, method=Image.FASTOCTREE)
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=True)
    return buf.getvalue()


def one(path: Path) -> dict | None:
    import anndata as ad

    try:
        a = ad.read_h5ad(path, backed="r")
    except Exception as e:                                  # noqa: BLE001
        log.debug("%s: unreadable (%s)", path.name, e)
        return None
    if "spatial" not in a.obsm:
        return None
    try:
        xy = np.asarray(a.obsm["spatial"])
    except Exception:                                       # noqa: BLE001
        return None
    labels = None
    if "cell_type" in a.obs.columns:
        try:
            labels = a.obs["cell_type"].astype(str).values
        except Exception:                                   # noqa: BLE001
            labels = None
    png = render(xy, labels)
    if png is None:
        return None

    def col(name):
        try:
            v = a.obs[name].astype(str)
            return str(v.iloc[0]) if len(v) else ""
        except Exception:                                   # noqa: BLE001
            return ""

    return {
        "png": "data:image/png;base64," + base64.b64encode(png).decode(),
        "n": int(a.shape[0]),
        "tissue": col("tissue"),
        "tech": col("technology_resolved"),
        "disease": col("disease_state"),
        "accession": col("accession"),
        "source": col("source"),
        "typed": bool(labels is not None),
    }


def pick(limit: int) -> list[Path]:
    """A spread, not the top of a sorted list.

    Taking the biggest N sections would return one technology and three tissues.
    Instead the files are grouped by (technology, tissue) as the filename and size
    allow, and the groups are visited round-robin, largest file first within each,
    so the gallery spans the catalogue before it goes deep on any part of it.
    """
    # Smallest first, deliberately. Every section is subsampled to MAX_PTS before it
    # is drawn, so a 7 GB file produces exactly the same picture as a 7 MB one and
    # costs a hundred times more to open. Sorting the other way spent fifteen minutes
    # on the largest files in the catalogue and produced nothing.
    files = sorted(STD.glob("*.h5ad"), key=lambda p: p.stat().st_size)
    groups: dict[str, list[Path]] = {}
    for f in files:
        # the source prefix is the cheapest grouping available without opening
        # every file; the per-file metadata read happens later and only for the
        # files that are actually drawn
        key = f.name.split("_")[0]
        groups.setdefault(key, []).append(f)
    out, i = [], 0
    while len(out) < limit:
        added = False
        for key in sorted(groups):
            if i < len(groups[key]):
                out.append(groups[key][i])
                added = True
                if len(out) >= limit:
                    break
        if not added:
            break
        i += 1
    return out


def build(limit: int = 240, force: bool = False) -> dict:
    """Thumbnails keyed by the standardised file's stem, cached between runs."""
    cache = {}
    if CACHE.exists() and not force:
        try:
            cache = json.loads(CACHE.read_text())
        except Exception:                                   # noqa: BLE001
            cache = {}
    want = pick(limit)
    todo = [f for f in want if f.stem not in cache]
    if todo:
        log.info("rendering %d thumbnails (%d already cached)", len(todo), len(want) - len(todo))
        t0, made, skipped = time.time(), 0, 0
        CACHE.parent.mkdir(parents=True, exist_ok=True)

        def flush():
            tmp = CACHE.with_suffix(".json.tmp")
            tmp.write_text(json.dumps(cache))
            tmp.replace(CACHE)

        for n, f in enumerate(todo, 1):
            r = one(f)
            if r:
                cache[f.stem] = r
                made += 1
            else:
                # a negative result is cached too, so a section with no usable
                # coordinates is not re-opened on every export
                cache[f.stem] = None
                skipped += 1
            # Written as it goes, not at the end. The first run of this was killed by a
            # timeout at fifteen minutes and lost every thumbnail it had rendered,
            # because the only write was after the loop.
            if n % 20 == 0:
                flush()
                log.info("  %d/%d  (%.1f s, cached)", n, len(todo), time.time() - t0)
        flush()
        log.info("  %d drawn, %d skipped for unusable coordinates, %.1f s",
                 made, skipped, time.time() - t0)
    out = {f.stem: cache[f.stem] for f in want if cache.get(f.stem)}
    kb = sum(len(v["png"]) for v in out.values()) / 1024
    log.info("%d thumbnails, %.0f KB of data URIs", len(out), kb)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=240)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--sheet", help="write a contact sheet PNG here to eyeball them")
    a = ap.parse_args()
    t = build(a.limit, a.force)
    if a.sheet:
        from PIL import Image
        items = list(t.values())[:64]
        cols = 8
        rows = (len(items) + cols - 1) // cols
        sheet = Image.new("RGBA", (cols * SIZE, rows * SIZE), (255, 255, 255, 255))
        for i, v in enumerate(items):
            raw = base64.b64decode(v["png"].split(",", 1)[1])
            sheet.paste(Image.open(io.BytesIO(raw)).convert("RGBA"),
                        ((i % cols) * SIZE, (i // cols) * SIZE))
        sheet.save(a.sheet)
        print("wrote", a.sheet)
