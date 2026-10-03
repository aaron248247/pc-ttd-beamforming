"""Record interpreter and dependency versions to environment/env_versions.json."""
import json, platform, subprocess, sys, io, contextlib
import numpy, scipy, matplotlib, pandas, pytest
import meep
from meep import mpb  # noqa: F401  (fails loudly if MPB is missing)

info = {
    "python": sys.version,
    "platform": platform.platform(),
    "numpy": numpy.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
    "pandas": pandas.__version__, "pytest": pytest.__version__, "meep": meep.__version__,
    "meep_with_mpi": bool(meep.with_mpi()), "mpb_available": True,
}
try:
    info["git_head"] = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
except Exception:
    info["git_head"] = None
json.dump(info, open("environment/env_versions.json", "w"), indent=2)
print(json.dumps(info, indent=2))
