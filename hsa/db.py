"""SQLite catalogue: sources -> candidates -> datasets -> samples -> files."""
import json
import sqlite3
import threading
import time

from .config import DB_PATH

SCHEMA = """
CREATE TABLE IF NOT EXISTS sources (
  name TEXT PRIMARY KEY, url TEXT, access_method TEXT, est_human_datasets INTEGER,
  notes TEXT, updated REAL);
CREATE TABLE IF NOT EXISTS candidates (
  source TEXT, accession TEXT, title TEXT, url TEXT, reason TEXT,
  status TEXT DEFAULT 'pending',            -- pending | curated | failed
  added REAL, PRIMARY KEY (source, accession));
CREATE TABLE IF NOT EXISTS datasets (
  source TEXT, accession TEXT, title TEXT, url TEXT, publication TEXT,
  technology TEXT, tissue TEXT, disease TEXT, n_samples INTEGER,
  verdict TEXT,      -- ready | bundle_only | restricted | not_human | not_spatial | no_processed_data
  notes TEXT, links TEXT, curated_by TEXT, updated REAL,
  PRIMARY KEY (source, accession));
CREATE TABLE IF NOT EXISTS files (
  source TEXT, accession TEXT, sample_id TEXT, url TEXT, role TEXT,  -- matrix | coords | matrix+coords | bundle
  fmt TEXT, size_bytes INTEGER, local_path TEXT,
  dl_status TEXT DEFAULT 'pending',          -- pending | done | skipped_large | failed | deferred
  dl_msg TEXT, PRIMARY KEY (url));
CREATE TABLE IF NOT EXISTS usage (
  run TEXT, ts REAL, input INTEGER, output INTEGER, cache_read INTEGER, cache_write INTEGER, usd REAL);
"""

_lock = threading.Lock()


def connect():
    con = sqlite3.connect(DB_PATH, timeout=60, check_same_thread=False)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA journal_mode=WAL")
    con.executescript(SCHEMA)
    return con


CON = connect()


def execute(sql, args=()):
    with _lock:
        cur = CON.execute(sql, args)
        CON.commit()
        return cur


def query(sql, args=()):
    with _lock:
        return [dict(r) for r in CON.execute(sql, args).fetchall()]


def upsert_source(name, url, access_method, est, notes):
    execute("INSERT OR REPLACE INTO sources VALUES (?,?,?,?,?,?)",
            (name, url, access_method, est, notes, time.time()))


def add_candidates(source, items, reason):
    n = 0
    for it in items:
        cur = execute("INSERT OR IGNORE INTO candidates (source, accession, title, url, reason, added) "
                      "VALUES (?,?,?,?,?,?)",
                      (source, it["accession"], it.get("title", ""), it.get("url", ""), reason, time.time()))
        n += cur.rowcount
    return n


def record_dataset(d, curated_by):
    execute("INSERT OR REPLACE INTO datasets VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (d["source"], d["accession"], d.get("title", ""), d.get("url", ""), d.get("publication", ""),
             d.get("technology", ""), d.get("tissue", ""), d.get("disease", ""), d.get("n_samples", 0),
             d["verdict"], d.get("notes", ""), json.dumps(d.get("links", [])), curated_by, time.time()))
    # a re-record replaces the not-yet-downloaded file list
    execute("DELETE FROM files WHERE source=? AND accession=? AND dl_status='pending'", (d["source"], d["accession"]))
    n = 0
    for s in d.get("samples", []):
        for f in s.get("files", []):
            execute("INSERT OR IGNORE INTO files (source, accession, sample_id, url, role, fmt, size_bytes) "
                    "VALUES (?,?,?,?,?,?,?)",
                    (d["source"], d["accession"], s["sample_id"], f["url"], f["role"], f.get("fmt", ""),
                     f.get("size_bytes") or None))
            n += 1
    execute("UPDATE candidates SET status='curated' WHERE source=? AND accession=?",
            (d["source"], d["accession"]))
    return n


def log_usage(run, u, usd):
    execute("INSERT INTO usage VALUES (?,?,?,?,?,?,?)",
            (run, time.time(), u["in"], u["out"], u["cache_read"], u["cache_write"], usd))


def total_usd():
    """Agent (curation/scout) spend only; metadata harmonisation ('meta_*' runs) has its own cap."""
    return query("SELECT COALESCE(SUM(usd),0) AS s FROM usage WHERE run NOT LIKE 'meta_%'")[0]["s"]
