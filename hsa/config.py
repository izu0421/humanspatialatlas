"""HSA (Human Spatial Atlas) configuration."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # .../space_agent/hsa
DATA = ROOT / "data"                                 # downloaded matrices
RUNS = ROOT / "runs"                                 # agent transcripts, logs
DB_PATH = ROOT / "hsa_catalogue.sqlite"
KEY_FILE = ROOT.parent / "api.txt"

MODEL = "claude-sonnet-5-5"
# $ per 1M tokens (Sonnet 5.5): input, output, cache read, cache write (5 min)
PRICE = {"in": 2.00, "out": 10.00, "cache_read": 0.20, "cache_write": 2.50}

# download guards
MAX_FILE_GB = 5.0          # skip single files larger than this
MIN_FREE_TB = 1.0          # stop downloading if /data free space drops below this
HTTP_UA = "HSA-human-spatial-atlas/0.1 (academic; yzy21@cam.ac.uk)"


def load_key() -> None:
    if not os.environ.get("ANTHROPIC_API_KEY"):
        os.environ["ANTHROPIC_API_KEY"] = KEY_FILE.read_text().strip()
