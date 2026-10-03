import os
os.environ.setdefault("OMP_NUM_THREADS", "1")      # one thread per worker; parallelism is across designs
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
RES = ROOT / "results" / "phase0p5"

# declared study parameters (docs/phase0p5_criteria.md)
EPS = 12.0
M = {"triangular": 5, "square": 4}
F1 = [0.8, 0.9, 1.0, 1.1, 1.2]
S = [-0.10, -0.05, 0.0, 0.05, 0.10]
SWEEP_RES = 24
SWEEP_TOL = 1e-7
BULK_RES = 32
