import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
P4 = ROOT / "results" / "phase4"
P5 = ROOT / "results" / "phase5"
RES = ROOT / "results" / "phase6"
FIG = ROOT / "figures" / "phase6"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

import numpy as np

# Pre-declared in docs/phase2_5_criteria.md, Revision 7
THETA0 = 30.0
MSUBS = [1, 2, 4, 8]
NQ = 4
CASES = ["fixed", "ideal_elem", "ideal_sub", "raw", "trim", "trim_p4"]
DIAG = ["raw_unit", "trim_unit"]
LABEL = {"fixed": "fixed phase", "ideal_elem": "ideal element-level TTD", "ideal_sub": "ideal subarray TTD",
         "raw": "raw PhC subarray TTD", "trim": "phase-intercept-trimmed PhC (re-trimmed to required step)",
         "trim_p4": "Phase 4 H_trim as is (trimmed for 3.9323 a)", "raw_unit": "raw PhC, |H| = 1 (diagnostic)", "trim_unit": "re-trimmed PhC, |H| = 1 (diagnostic)"}
SHORT = {"fixed": "fixed phase", "ideal_elem": "ideal elem. TTD", "ideal_sub": "ideal subarray TTD", "raw": "raw PhC", "trim": "trimmed PhC", "trim_p4": "P4 H_trim as is"}
COL = {"fixed": "#d62728", "ideal_elem": "#000000", "ideal_sub": "#7f7f7f", "raw": "#ff7f0e", "trim": "#1f77b4", "trim_p4": "#2ca02c", "raw_unit": "#ff7f0e", "trim_unit": "#1f77b4"}
LS = {"fixed": "-", "ideal_elem": "--", "ideal_sub": "-.", "raw": "-", "trim": "-", "trim_p4": ":", "raw_unit": ":", "trim_unit": ":"}
MCOL = {1: "#9467bd", 2: "#8c564b", 4: "#1f77b4", 8: "#e377c2"}
