import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
P4 = ROOT / "results" / "phase4"
RES = ROOT / "results" / "phase5"
FIG = ROOT / "figures" / "phase5"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

import numpy as np

# Pre-declared in docs/phase2_5_criteria.md, Revision 6
THETA0 = 30.0
CASES = ["fixed", "ideal", "raw", "trim"]
LABEL = {"fixed": "fixed phase", "ideal": "ideal TTD", "raw": "raw PhC", "trim": "phase-intercept-trimmed PhC",
         "raw_unit": "raw PhC, |H| = 1 (diagnostic)", "trim_unit": "trimmed PhC, |H| = 1 (diagnostic)"}
COL = {"fixed": "#d62728", "ideal": "#000000", "raw": "#ff7f0e", "trim": "#1f77b4", "raw_unit": "#ff7f0e", "trim_unit": "#1f77b4"}
LS = {"fixed": "-", "ideal": "--", "raw": "-", "trim": "-", "raw_unit": ":", "trim_unit": ":"}
