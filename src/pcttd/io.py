"""Result I/O with provenance. Every saved result records code version, library versions and a timestamp."""
from __future__ import annotations

import datetime
import json
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def _git(*args):
    try:
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


def provenance() -> dict:
    import meep, scipy
    dirty = _git("status", "--porcelain")
    return dict(
        timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        git_head=_git("rev-parse", "HEAD"), git_dirty=bool(dirty),
        meep=meep.__version__, numpy=np.__version__, scipy=scipy.__version__,
    )


class _Enc(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def save_json(path, obj, with_provenance=True):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if with_provenance and isinstance(obj, dict):
        obj = {**obj, "provenance": provenance()}
    path.write_text(json.dumps(obj, indent=1, cls=_Enc))


def load_json(path):
    return json.loads(Path(path).read_text())
