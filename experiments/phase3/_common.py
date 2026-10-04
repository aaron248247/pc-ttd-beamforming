import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
RES = ROOT / "results" / "phase3"
RAW = RES / "raw"                       # raw plane-field dumps (large, not tracked by git)
RES.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

import multiprocessing as mp_
import numpy as np


def pool_map(fn, items, n=8):
    with mp_.get_context("spawn").Pool(min(n, len(items))) as p:
        return p.map(fn, items, chunksize=1)
