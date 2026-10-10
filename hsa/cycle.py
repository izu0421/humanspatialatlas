"""One maintenance cycle: find what is new, take it, make it match the rest, report.

Run on a schedule, this is what makes the Atlas track the literature rather than being a
snapshot. Each stage is already resumable and each writes its own ledger, so a cycle that dies
halfway costs time and not state.

Deliberately conservative about two things:
  * money -- the agent stages stop at a per-cycle budget AND an absolute lifetime cap, so a
    scheduled job cannot quietly drain the account if a lane starts looping.
  * disk -- materialisation refuses to start below a floor, because deleting raw files to make
    room for their replacements is only safe while there is somewhere to put the replacements.
"""
from __future__ import annotations

import logging
import os
import subprocess
import sys
import time
from pathlib import Path

from . import db
from .config import DATA, ROOT, RUNS

logger = logging.getLogger(__name__)
LOGDIR = RUNS / "cycles"


def _free_gb() -> float:
    st = os.statvfs(str(DATA))
    return st.f_bavail * st.f_frsize / 1e9


def _run(label: str, args: list[str], log, timeout_s: int) -> int:
    log.write(f"\n===== {label} :: {' '.join(args)}\n")
    log.flush()
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, *args], cwd=str(ROOT), stdout=log, stderr=log,
                           timeout=timeout_s)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        log.write(f"----- {label} TIMED OUT after {timeout_s}s\n")
        rc = -1
    log.write(f"----- {label} rc={rc} in {time.time()-t0:.0f}s\n")
    log.flush()
    return rc


def cycle(budget_usd: float = 40.0, lifetime_cap_usd: float = 500.0,
          min_free_gb: float = 150.0, scout: bool = True) -> dict:
    """Scout -> curate -> download -> QC -> harmonise -> materialise -> export."""
    LOGDIR.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M")
    out = {"started": stamp, "spend_before": db.total_usd(), "free_gb_before": round(_free_gb())}

    if out["spend_before"] >= lifetime_cap_usd:
        out["skipped"] = f"lifetime cap reached (${out['spend_before']:.2f})"
        return out

    with open(LOGDIR / f"cycle_{stamp}.log", "w") as log:
        log.write(f"HSA cycle {stamp}\nspend ${out['spend_before']:.2f} / {lifetime_cap_usd}\n"
                  f"free {out['free_gb_before']} GB\n")
        # the per-cycle budget is the smaller of the allowance and the headroom to the cap
        b = min(budget_usd, max(0.0, lifetime_cap_usd - out["spend_before"]))
        if scout and b > 1:
            _run("scout", ["run_hsa.py", "scout", "--budget", str(b), "--workers", "6"], log, 7200)
        if db.total_usd() < lifetime_cap_usd:
            _run("curate", ["run_hsa.py", "curate", "--budget", str(lifetime_cap_usd),
                            "--workers", "6"], log, 10800)
        _run("download", ["run_hsa.py", "download", "--workers", "6"], log, 21600)
        _run("qc", ["-m", "hsa.standardise"], log, 21600)
        _run("panels", ["-m", "hsa.panels"], log, 3600)
        if _free_gb() >= min_free_gb:
            _run("materialise", ["-c",
                 "from hsa import materialise as M; M.convert_all(delete_raw=True, min_free_gb=%d)"
                 % min_free_gb], log, 43200)
        else:
            log.write(f"\n===== materialise SKIPPED: only {_free_gb():.0f} GB free\n")
        _run("export", ["-m", "hsa.export"], log, 1800)

    out.update(spend_after=db.total_usd(), free_gb_after=round(_free_gb()),
               log=str(LOGDIR / f"cycle_{stamp}.log"))
    out["spent"] = round(out["spend_after"] - out["spend_before"], 2)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=float, default=40.0, help="$ for this cycle's agent stages")
    ap.add_argument("--cap", type=float, default=500.0, help="absolute lifetime $ cap")
    ap.add_argument("--min-free-gb", type=float, default=150.0)
    ap.add_argument("--no-scout", action="store_true")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    r = cycle(a.budget, a.cap, a.min_free_gb, scout=not a.no_scout)
    print(r)
