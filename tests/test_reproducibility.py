"""Reproducibility: identical configurations give identical results, and fresh runs reproduce the saved results. Marked ``meep``."""
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("meep")
from pcttd import fdtd1d as F1  # noqa: E402
from pcttd import fdtd2d as F2  # noqa: E402

pytestmark = pytest.mark.meep
RES = Path(__file__).resolve().parents[1] / "results" / "phase1"


def test_1d_run_is_repeatable_bit_for_bit():
    cfg = F1.Cfg1D(res=20)
    a, b = F1.run_1d(cfg), F1.run_1d(cfg)
    assert np.array_equal(a["E"], b["E"]) and a["t_end"] == b["t_end"]


@pytest.mark.skipif(not (RES / "1a_base_curves.npz").exists(), reason="saved 1A result not present")
def test_fresh_1d_vacuum_run_reproduces_saved_1a_curves():
    z = np.load(RES / "1a_base_curves.npz")
    r = F1.run_1d(F1.Cfg1D())                                   # the 1A base configuration
    assert np.array_equal(r["f"], z["f"])
    assert np.max(np.abs(r["E"] - z["E"])) <= 1e-12 * np.max(np.abs(z["E"]))


def test_2d_run_is_repeatable_bit_for_bit():
    cfg = F2.Cfg2D(N=4, ny=24, decay_tol=1e-6)
    a, b = F2.run_2d(cfg), F2.run_2d(cfg)
    assert np.array_equal(a["T"], b["T"]) and np.array_equal(a["R"], b["R"])


@pytest.mark.skipif(not (RES / "1d_pc2d.json").exists(), reason="saved 1D result not present")
def test_fresh_2d_run_reproduces_saved_pc2d_curve():
    saved = json.loads((RES / "1d_pc2d.json").read_text())
    row = [r for r in saved["length"] if r["N"] == 2][0]        # cheapest saved case: N = 2, ny = 40, default settings
    cfg = replace(F2.Cfg2D(N=2, ny=40))
    r = F2.run_2d(cfg)
    assert np.allclose(r["T"], row["curve"]["T"], rtol=1e-10, atol=1e-14)
    assert np.allclose(r["R"], row["curve"]["R"], rtol=1e-10, atol=1e-14)
