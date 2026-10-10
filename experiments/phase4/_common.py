import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
RES = ROOT / "results" / "phase4"
RAW = RES / "raw"                       # raw plane-field dumps (large, not tracked by git)
RES.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

import multiprocessing as mp_
import numpy as np

# Pre-declared in docs/phase2_5_criteria.md, Revision 4 (2026-10-10)
N0, DN, NCH = 8, 1, 4
LENGTHS = [N0 + m * DN for m in range(NCH)]
F0 = 0.2495
TAU_PP = 3.9323                          # Phase 3 D0 median per-period delay (a/c), validated +-10% band, res 30
DTAU_TARGET = DN * TAU_PP
VALID = (0.231, 0.268)


def pool_map(fn, items, n=8):
    with mp_.get_context("spawn").Pool(min(n, len(items))) as p:
        return p.map(fn, items, chunksize=1)
