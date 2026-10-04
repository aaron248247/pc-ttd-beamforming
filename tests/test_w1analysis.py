import numpy as np
import pytest

from pcttd import extract as X
from pcttd import guided as GD
from pcttd import w1analysis as WA


def test_linear_length_fit_recovers_slope_intercept_and_flags_nonlinearity():
    N = np.array([6, 12, 24])
    s = np.array([0.8, 1.9]); c = np.array([0.3, -0.2])
    P = c[None, :] + s[None, :] * N[:, None]
    sl, ic, rms = WA.linear_length_fit(N, P)
    assert np.allclose(sl, s) and np.allclose(ic, c) and np.allclose(rms, 0, atol=1e-12)
    P[1] += 0.05                                         # a bump in the middle length is not removable by a line
    assert np.all(WA.linear_length_fit(N, P)[2] > 0.01)


def test_cascade_theta_recovers_phase_of_added_section_with_mismatched_interfaces():
    # device = interface * propagation(N) * interface^-1, lossless: S11, S21 from the T product
    rng = np.random.default_rng(3)
    R = 0.3; S11i = np.sqrt(R) * np.exp(0.7j); S21i = np.sqrt(1 - R) * np.exp(-0.4j)
    Ti = GD.t_matrix(S11i, S21i, GD.s22_from_unitarity(S11i, S21i))
    k = np.linspace(0.9, 1.4, 11)
    S11, S21 = {}, {}
    for N in (4, 9):
        T = np.array([Ti @ GD.propagate_T(kk, N) @ np.linalg.inv(Ti) for kk in k])
        s11, s21, _, _ = GD.s_from_t(T)
        S11[N], S21[N] = s11, s21
    th, mod = WA.cascade_theta(S11, S21, 4, 9)
    assert np.allclose(th, np.angle(np.exp(1j * 5 * k)), atol=1e-9) and np.allclose(mod, 1)


def test_ng_from_theta_with_meep_sign_and_subsample():
    f = np.linspace(0.2, 0.3, 401)
    ng0 = 4.0
    theta = X.omega(f) * ng0
    assert np.allclose(WA.ng_from_theta(f, theta)[2:-2], ng0)
    f2, (t2,) = WA.subsample(f, [theta], 4)
    assert np.allclose(WA.ng_from_theta(f2, t2)[2:-2], ng0) and len(f2) == 101


def test_ng_vs_reference_stats():
    r = WA.ng_vs_reference(np.arange(3), np.array([4.0, 4.4, 3.8]), np.array([4.0, 4.0, 4.0]))
    assert np.allclose(r["rel"], [0, 0.1, -0.05]) and r["max_abs_rel"] == pytest.approx(0.1)


def test_pulse_centroid_delay_of_shifted_gaussian():
    t = np.linspace(0, 100, 4001)
    g = lambda t0: np.exp(-((t - t0) / 5) ** 2) * np.cos(2 * (t - t0))
    assert WA.pulse_centroid_delay(t, g(30), g(41.5)) == pytest.approx(11.5, abs=1e-3)


def test_unitarity_defect_zero_for_lossless():
    assert np.allclose(WA.unitarity_defect(np.sqrt(0.2), np.sqrt(0.8)), 0)
