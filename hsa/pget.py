"""Parallel range-request downloader.

NCBI throttles a single HTTP connection to about 70 KB/s, which makes a 2 TB transfer a
331-day job; eight concurrent range requests to the same file sustain ~26 MB/s, bringing it to
about a day. Each chunk is written to its own part file and concatenated only once every part
has arrived, so an interrupted transfer resumes per chunk rather than restarting.
"""
from __future__ import annotations

import logging
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

logger = logging.getLogger(__name__)


def size_of(url: str, timeout: int = 60) -> int | None:
    r = subprocess.run(["curl", "-sfIL", "--max-time", str(timeout), url],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            return int(line.split(":", 1)[1].strip())
    return None


def _chunk(url: str, out: Path, start: int, end: int, tries: int = 4) -> int:
    for attempt in range(tries):
        have = out.stat().st_size if out.exists() else 0
        want = end - start + 1
        if have >= want:
            return have
        r = subprocess.run(["curl", "-sfL", "--max-time", "3600",
                            "-r", f"{start + have}-{end}", "-o", str(out),
                            *(["-C", "-"] if have else []), url],
                           capture_output=True)
        if r.returncode == 0 and out.exists() and out.stat().st_size >= want:
            return out.stat().st_size
        time.sleep(2 ** attempt)
    return out.stat().st_size if out.exists() else 0


def get(url: str, dest: str | Path, conns: int = 8, chunk_mb: int = 256) -> Path:
    dest = Path(dest)
    if dest.exists():
        return dest
    total = size_of(url)
    if not total:
        subprocess.run(["curl", "-sfL", "-o", str(dest), url], check=True)
        return dest
    part_dir = dest.with_suffix(dest.suffix + ".parts")
    part_dir.mkdir(parents=True, exist_ok=True)
    step = chunk_mb * 1024 * 1024
    spans = [(i, s, min(s + step - 1, total - 1))
             for i, s in enumerate(range(0, total, step))]
    t0 = time.time()
    with ThreadPoolExecutor(conns) as ex:
        futs = {ex.submit(_chunk, url, part_dir / f"{i:05d}", s, e): (i, s, e)
                for i, s, e in spans}
        done = 0
        for f in as_completed(futs):
            f.result(); done += 1
            if done % max(1, len(spans) // 8) == 0:
                got = sum(p.stat().st_size for p in part_dir.glob("*"))
                logger.info("  %s: %.0f%% (%.1f MB/s)", dest.name, got / total * 100,
                            got / 1e6 / max(time.time() - t0, 1))
    got = sum(p.stat().st_size for p in part_dir.glob("*"))
    if got < total:
        raise IOError(f"{dest.name}: got {got:,} of {total:,} bytes")
    tmp = dest.with_suffix(dest.suffix + ".joining")
    with open(tmp, "wb") as o:
        for p in sorted(part_dir.glob("*")):
            with open(p, "rb") as f:
                while (b := f.read(1 << 22)):
                    o.write(b)
    tmp.rename(dest)
    for p in part_dir.glob("*"):
        p.unlink()
    part_dir.rmdir()
    logger.info("%s: %.1f GB in %.0fs (%.1f MB/s)", dest.name, total / 1e9,
                time.time() - t0, total / 1e6 / max(time.time() - t0, 1))
    return dest
