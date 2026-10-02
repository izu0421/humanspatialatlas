"""Deterministic downloader: only matrix / coords / matrix+coords files. Resumable, size-capped."""
import os
import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import db
from .config import DATA, MAX_FILE_GB, MIN_FREE_TB
from .tools_api import S, ZIP_SEP, _throttle, zip_extract

GET_ROLES = ("matrix", "coords", "matrix+coords")


def _safe(s):
    return re.sub(r"[^A-Za-z0-9._\-]+", "_", s)[:150]


def _dest(f):
    name = (f["url"].split(ZIP_SEP, 1)[1] if ZIP_SEP in f["url"] else f["url"].split("?")[0]).rstrip("/").rsplit("/", 1)[-1]
    if "cellxgene" in f["url"] and not name.endswith(".h5ad"):
        name += ".h5ad"
    return DATA / _safe(f["source"]) / _safe(f["accession"]) / _safe(f["sample_id"]) / _safe(name)


def _one(f):
    dest = _dest(f)
    if dest.exists():
        return "done", str(dest), "exists"
    if shutil.disk_usage(DATA).free < MIN_FREE_TB * 1e12:
        return "pending", None, "disk floor reached"
    size = f["size_bytes"]
    if size is None and ZIP_SEP not in f["url"]:
        try:
            h = S.head(f["url"].split(ZIP_SEP)[0], timeout=60, allow_redirects=True)
            size = int(h.headers.get("Content-Length", 0)) or None
        except Exception:
            size = None
    if size and size > MAX_FILE_GB * 1e9:
        return "skipped_large", None, f"{size/1e9:.1f} GB"
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    if ZIP_SEP in f["url"]:            # member of a remote zip: range-read only that member
        zip_extract(f["url"], part)
        os.replace(part, dest)
        return "done", str(dest), f"{dest.stat().st_size/1e6:.1f} MB (zip member)"
    have = part.stat().st_size if part.exists() else 0
    hdr = {"Range": f"bytes={have}-"} if have else {}
    _throttle(f["url"])
    with S.get(f["url"], stream=True, timeout=120, headers=hdr) as r:
        if r.status_code not in (200, 206):
            return "failed", None, f"HTTP {r.status_code}"
        mode = "ab" if r.status_code == 206 else "wb"
        with open(part, mode) as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
    os.replace(part, dest)
    return "done", str(dest), f"{dest.stat().st_size/1e6:.1f} MB"


def run(workers=8, source=None, limit=None):
    db.execute("UPDATE files SET dl_status='deferred' WHERE role='bundle' AND dl_status='pending'")
    sql = "SELECT * FROM files WHERE dl_status='pending' AND role IN ('matrix','coords','matrix+coords')"
    args = []
    if source:
        sql += " AND source=?"
        args.append(source)
    todo = db.query(sql + (f" LIMIT {int(limit)}" if limit else ""), args)
    print(f"{len(todo)} files to fetch")
    ok = 0
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(_one, f): f for f in todo}
        for i, fu in enumerate(as_completed(futs)):
            f = futs[fu]
            try:
                status, path, msg = fu.result()
            except Exception as e:
                status, path, msg = "failed", None, str(e)[:300]
            db.execute("UPDATE files SET dl_status=?, local_path=?, dl_msg=? WHERE url=?",
                       (status, path, msg, f["url"]))
            ok += status == "done"
            if i % 50 == 0:
                print(f"[{i+1}/{len(todo)}] {status} {f['accession']} {f['url'].rsplit('/',1)[-1]} {msg}", flush=True)
    print(f"done: {ok}/{len(todo)}")


# ---------------------------------------------------------------- bundles (tar / tar.gz / zip)
import tarfile

from .tools_api import classify, zip_list

MAX_BUNDLE_GB = 60.0
TAR_SEP = "#tarmember="


def _member_sample(name, default):
    m = re.search(r"(GSM\d+)", name)
    return m.group(1) if m else default


def _add_member(b, member_url, name, role, path, size):
    db.execute("INSERT OR REPLACE INTO files (source, accession, sample_id, url, role, fmt, size_bytes, local_path, "
               "dl_status, dl_msg) VALUES (?,?,?,?,?,?,?,?,?,?)",
               (b["source"], b["accession"], _member_sample(name, b["sample_id"]), member_url, role,
                name.rsplit(".", 1)[-1], size, str(path) if path else None,
                "done" if path else "pending", f"from bundle {b['url'].rsplit('/', 1)[-1]}"))


def _bundle(b):
    url = b["url"]
    base = DATA / _safe(b["source"]) / _safe(b["accession"])
    if url.lower().endswith(".zip"):
        files = [f for f in zip_list(url)["files"] if f["role"] in GET_ROLES]
        for f in files:
            dest = base / _safe(_member_sample(f["name"], b["sample_id"])) / _safe(f["name"].rsplit("/", 1)[-1])
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists():
                zip_extract(f["url"], dest)
            _add_member(b, f["url"], f["name"], f["role"], dest, f["size_bytes"])
        return "done", None, f"zip: {len(files)} members"
    # tar / tar.gz: one streaming pass, keep only matrix/coords members (nested .tar.gz members are kept whole)
    n = 0
    with S.get(url, stream=True, timeout=300) as r:
        if r.status_code != 200:
            return "failed", None, f"HTTP {r.status_code}"
        size = int(r.headers.get("Content-Length", 0))
        if size > MAX_BUNDLE_GB * 1e9:
            return "skipped_large", None, f"{size/1e9:.1f} GB bundle"
        r.raw.decode_content = True
        with tarfile.open(fileobj=r.raw, mode="r|*") as tf:
            for m in tf:
                if not m.isfile():
                    continue
                role = classify(m.name)
                if role not in GET_ROLES:
                    continue
                dest = base / _safe(_member_sample(m.name, b["sample_id"])) / _safe(m.name.rsplit("/", 1)[-1])
                dest.parent.mkdir(parents=True, exist_ok=True)
                with tf.extractfile(m) as src, open(dest, "wb") as out:
                    shutil.copyfileobj(src, out, 1 << 22)
                _add_member(b, f"{url}{TAR_SEP}{m.name}", m.name, role, dest, m.size)
                n += 1
    return "done", None, f"tar: {n} members kept"


def run_bundles(workers=4, source=None, limit=None):
    sql = "SELECT * FROM files WHERE role='bundle' AND dl_status IN ('deferred','pending')"
    args = []
    if source:
        sql += " AND source=?"
        args.append(source)
    todo = db.query(sql + (f" LIMIT {int(limit)}" if limit else ""), args)
    print(f"{len(todo)} bundles", flush=True)
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(_bundle, b): b for b in todo}
        for fu in as_completed(futs):
            b = futs[fu]
            try:
                status, _, msg = fu.result()
            except Exception as e:
                status, msg = "failed", f"{type(e).__name__}: {e}"[:300]
            db.execute("UPDATE files SET dl_status=?, dl_msg=? WHERE url=?", (status, msg, b["url"]))
            print(f"{status} {b['accession']} {b['url'].rsplit('/',1)[-1]} {msg}", flush=True)
