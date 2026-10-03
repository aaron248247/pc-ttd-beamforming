"""MPB-dependent verification against closed-form results (homogeneous media). Marked ``mpb``."""
import numpy as np
import pytest

pytest.importorskip("meep")
from pcttd import mpb_bands as MB  # noqa: E402
from pcttd import geometry as G  # noqa: E402

pytestmark = pytest.mark.mpb


@pytest.fixture(scope="module")
def uniform_eps12():
    # triangular host (eps=12), no holes, M=3 supercell -> a plain eps=12 slab-free medium
    H = 3 * G.SQRT3
    return MB.defect_bands("triangular", 3, 0.3, res=16, nk=11, _override=(H, [])), H


def test_group_index_equals_sqrt_eps(uniform_eps12):
    d, _ = uniform_eps12
    # lowest band at k>0: plane wave with ky=0: w = k/sqrt(eps); v_g = 1/sqrt(eps)
    vg = d["vg"][1:, 0]
    assert np.allclose(1.0 / np.abs(vg), np.sqrt(12.0), rtol=1e-4)
    w = d["freq"][1:, 0]
    assert np.allclose(w, d["k"][1:] / np.sqrt(12.0), rtol=1e-4)


def test_lowest_band_is_even_and_degenerate_pair_has_both_parities(uniform_eps12):
    d, H = uniform_eps12
    assert (d["parity"][1:, 0] == +1).all()          # ky = 0 plane wave -> Phi even about y=0
    # next level: ky = +-1/H degenerate (cos, sin) -> one even, one odd
    ik = 5
    assert d["freq"][ik, 1] == pytest.approx(d["freq"][ik, 2], rel=1e-5)
    assert sorted(d["parity"][ik, 1:3]) == [-1, 1]
    expected = np.hypot(d["k"][ik], 1.0 / H) / np.sqrt(12.0)
    assert d["freq"][ik, 1] == pytest.approx(expected, rel=1e-4)


def test_strip_energy_fraction_uniform_field(uniform_eps12):
    d, H = uniform_eps12
    # uniform ky=0 wave -> energy fraction equals area fraction, up to the grid quantum (2 rows of ny = H*res pixels)
    quantum = 2.0 / (H * 16)
    expected = (2 * 1.5 * G.row_pitch("triangular")) / H
    assert d["strip"][3, 0] == pytest.approx(expected, abs=quantum)
    edge_expected = (2 * G.row_pitch("triangular")) / H
    assert d["edge"][3, 0] == pytest.approx(edge_expected, abs=quantum)


def test_vacuum_group_index_is_one():
    H = 9.0
    d = MB.defect_bands("square", 4, 0.2, res=16, nk=9, _override=(H, []))
    assert np.allclose(1.0 / np.abs(d["vg"][1:, 0]), 1.0, rtol=1e-4)
    assert (d["parity"] != 0).all()          # every band has a defined parity (guard for vanishing Phi)


def test_mpb_group_velocity_matches_finite_difference_in_real_structure():
    d = MB.defect_bands("triangular", 3, 0.3, res=16, nk=41)
    # lowest band, away from k=0 and the zone edge
    fd = np.gradient(d["freq"][:, 0], d["k"])
    sel = slice(4, 36)
    assert np.allclose(fd[sel], d["vg"][sel, 0], rtol=0.03, atol=2e-3)


def test_symmetric_structure_is_mirror_symmetric_on_the_grid():
    d = MB.defect_bands("triangular", 3, 0.3, res=16, nk=3)
    assert d["eps_mirror_asym"] < 1e-6
