"""Compact 2-D FDTD check of the finite triangular air-hole crystal against the MPB directional band structure (Phase 1D).

Marked ``meep``. The full convergence studies are in results/phase1; this is the regression guard.
MPB (results/phase1/1d_mpb_directional.json): y-even gap [0.20705, 0.2890]; all-parity gap top 0.2744.
"""
import numpy as np
import pytest

pytest.importorskip("meep")
from pcttd import fdtd2d as F2  # noqa: E402

pytestmark = pytest.mark.meep
EVEN_GAP = (0.20705, 0.2890)
ALL_PARITY_TOP = 0.2744


@pytest.fixture(scope="module")
def run():
    cfg = F2.Cfg2D(N=8, ny=24, decay_tol=1e-8)
    return cfg, F2.run_2d(cfg)


def test_frequency_axis_is_the_requested_one(run):
    cfg, r = run
    assert np.allclose(r["f"], np.array(cfg.freqs), atol=1e-12)


@pytest.mark.parametrize("ny", [24, 32, 40])
def test_cell_height_is_exactly_whole_pixels_so_lattice_is_not_strained(ny):
    # regression: an integer resolution made MEEP round the periodic height sqrt(3) a to 1.75 a (1% strain). With
    # resolution = ny / sqrt(3) the allocated pixel count across y is exactly ny, i.e. the height is exactly sqrt(3) a.
    cfg = F2.Cfg2D(ny=ny)
    sim = F2._build(cfg, True)
    sim.init_sim()
    assert sim.structure.gv.ny() == ny
    assert sim.structure.gv.nx() == round(sim.cell_size.x * cfg.res)           # x length is also a whole number of pixels


def test_stop_band_matches_mpb_within_documented_tolerance(run):
    _, r = run
    f, T = r["f"], r["T"]
    inside = (f >= 0.22) & (f <= 0.27)
    assert T[inside].max() < 1e-2                               # measured 1.5e-3 at N = 8 (termination resonance), all resolutions
    # a pass band exists below the gap: some frequency transmits almost fully
    assert T[(f >= 0.12) & (f <= 0.20)].max() > 0.9             # measured 0.99
    # above the even gap top the crystal transmits again
    assert T[(f >= 0.30) & (f <= 0.36)].max() > 0.5             # measured > 0.9
    # symmetry selection: stop band continues through the odd-parity-only window (all-parity gap would end at 0.2744)
    win = (f >= ALL_PARITY_TOP + 0.001) & (f <= EVEN_GAP[1] - 0.001)
    assert T[win].max() < 1e-2                                  # measured 1.5e-4 at N = 8


def test_energy_conservation_in_the_gap(run):
    _, r = run
    f = r["f"]
    gap = (f >= 0.22) & (f <= 0.27)
    assert np.max(np.abs(r["T"] + r["R"] - 1)[gap]) < 5e-3      # measured 8.5e-5 (N=8, ny=24, tol 1e-8)
