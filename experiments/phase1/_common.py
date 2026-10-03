import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
RES = ROOT / "results" / "phase1"
RES.mkdir(parents=True, exist_ok=True)

import multiprocessing as mp_
import numpy as np


def pool_map(fn, items, n=10):
    with mp_.get_context("spawn").Pool(min(n, len(items))) as p:
        return p.map(fn, items, chunksize=1)


def cfg_dict(cfg):
    d = asdict(cfg)
    d["freqs"] = [float(min(cfg.freqs)), float(max(cfg.freqs)), len(cfg.freqs)]
    return d
