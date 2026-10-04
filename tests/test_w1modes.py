"""``pcttd.w1modes``: band inversion (pure) and MPB checks of the guided band (marked ``mpb``)."""
import numpy as np
import pytest

pytest.importorskip("meep")
from pcttd import w1modes as WM  # noqa: E402


def _fake_bands():
    k = np.linspace(0, 0.5, 51)
    f = 0.28 - 0.06 * np.sin(np.pi * k) ** 1.0 * 0 - 0.06 * k * 2 + 0.02 * k**2   # monotone decreasing
    vg = np.gradient(f, k)
    vg[0] = vg[-1] = 0.0                                                          # flat ends as at band edges
    return dict(k=k, freq=f[:, None], vg=vg[:, None])


def test_kappa_of_f_inverts_a_monotone_decreasing_band():
    b = _fake_bands()
    kk = np.array([0.1, 0.2, 0.3, 0.4])
    f = np.interp(kk, b["k"], b["freq"][:, 0])
    assert np.allclose(WM.kappa_of_f(b, 1, f), kk, atol=2e-4)


def test_kappa_of_f_rejects_out_of_range_and_non_monotone_bands():
    b = _fake_bands()
    with pytest.raises(ValueError):
        WM.kappa_of_f(b, 1, [0.5])
    b2 = dict(b); b2["freq"] = b["freq"].copy(); b2["freq"][20:30, 0] += 0.05
    b2["vg"] = np.full_like(b["vg"], -0.1)
    with pytest.raises(ValueError):
        WM.kappa_of_f(b2, 1, [0.25])


@pytest.mark.mpb
def test_guided_band_decreasing_with_consistent_group_velocity_and_even_profile():
    d = WM.even_bands(M=6, res=16, nk=41, nb=9)
    f, vg, k = d["freq"][:, 7], d["vg"][:, 7], d["k"]
    mid = slice(8, 33)
    assert np.all(np.diff(f[mid]) < 0) and np.all(vg[mid] < 0)
    # Hellmann-Feynman v_g agrees with the centered finite difference of the band
    fd = np.gradient(f, k)
    assert np.allclose(vg[mid], fd[mid], rtol=0.05)
    P = WM.forward_profiles([0.25], 8, M=6, res=16)
    # Hz of the Hz-even class is symmetric about y = 0 (up to the sampling convention of the grid)
    hz = P["Hz"][0]
    assert np.allclose(np.abs(hz), np.abs(hz[::-1]), atol=0.02 * np.abs(hz).max()) or np.allclose(np.abs(hz[1:]), np.abs(hz[1:][::-1]), atol=0.02 * np.abs(hz).max())
    assert P["vg"][0] > 0 and f[mid].min() < P["freq"][0] < f[mid].max()
