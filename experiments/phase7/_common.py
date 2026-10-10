import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
RESULTS = ROOT / "results"
RES = RESULTS / "final"
FIG = ROOT / "figures" / "final"
DOCS = ROOT / "docs"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

import json
import numpy as np


def jload(rel):
    return json.load(open(RESULTS / rel))


def nload(rel):
    return np.load(RESULTS / rel, allow_pickle=True)


# Authoritative commits carried forward from the repository history (checked by tests/test_phase7.py against git when available)
COMMITS = {"phase6_final_head": "3c10780", "phase6_gate": "914628b", "phase6_criteria": "65f1ab3", "phase5_gate": "959e4c9", "phase5_criteria": "6553a1d"}

# Result files that Phase 7 must not modify (Phases 4-6)
PROTECTED = ["phase4/phase4_handoff.npz", "phase4/phase4_handoff.json", "phase4/4b_res30.npz", "phase4/4b_res30.json", "phase4/4b_res40.npz", "phase4/4b_res40.json",
             "phase5/5a_array.npz", "phase5/5a_array.json", "phase5/phase5_handoff.npz", "phase5/phase5_handoff.json",
             "phase6/6a_subarrays.npz", "phase6/6a_subarrays.json", "phase6/6b_robustness.npz", "phase6/6b_robustness.json",
             "phase6/phase6_handoff.npz", "phase6/phase6_handoff.json"]
