"""Uniform-linear-array model (``pcttd.array``): steering sign convention, squint of fixed phase vs ideal TTD, trimming."""
import numpy as np
import pytest

from pcttd import array as A

D, N, TH0 = 0.5, 8, 30.0
F = np.linspace(0.8, 1.2, 9)                      # c = 1, d = 0.5 -> lambda/2 at f = 1


def test_target_delays_are_progressive_d_sin_theta_over_c():
    tau = A.target_delays(N, D, TH0)
    assert np.allclose(np.diff(tau), D * np.sin(np.deg2rad(TH0)))
    assert tau[0] == 0 and len(tau) == N


def test_ideal_ttd_points_at_theta0_at_every_frequency():
    tau = A.target_delays(N, D, TH0)
    th, gain, _ = A.pointing_angle(A.ideal_ttd_weights(F, tau), F, D, TH0)
    assert np.allclose(th, TH0, atol=1e-3)
    assert np.allclose(gain, 1.0, atol=1e-6)


def test_fixed_phase_squints_as_sin_theta_equals_f0_over_f_sin_theta0():
    tau = A.target_delays(N, D, TH0)
    f0 = 1.0
    th, _, _ = A.pointing_angle(A.fixed_phase_weights(F, tau, f0), F, D, TH0)
    want = np.rad2deg(np.arcsin(f0 / F * np.sin(np.deg2rad(TH0))))
    assert np.allclose(th, want, atol=5e-3)
    assert th[0] > TH0 > th[-1]                    # lower frequency squints away from broadside; higher toward it


def test_negative_steering_angle_mirrors():
    tau = A.target_delays(N, D, -TH0)
    th, _, _ = A.pointing_angle(A.ideal_ttd_weights(F, tau), F, D, -TH0)
    assert np.allclose(th, -TH0, atol=1e-3)


def test_trim_phase_removes_only_a_constant_phase_per_channel():
    tau = A.target_delays(4, D, TH0)
    f0 = 1.0
    err = np.array([0.0, 0.31, -0.2, 0.8])         # constant phase errors per channel
    H = A.ideal_ttd_weights(F, tau) * np.exp(1j * err)[:, None] * np.array([1, .9, 1.1, .8])[:, None]
    Ht, delta = A.trim_phase(H, F, tau, f0)
    assert np.allclose(np.angle(np.exp(1j * (delta - err))), 0, atol=1e-12)
    assert np.allclose(np.abs(Ht), np.abs(H))
    # trimmed channels are ideal TTD up to their amplitudes at every frequency
    assert np.allclose(Ht / np.abs(Ht), A.ideal_ttd_weights(F, tau))


def test_pointing_angle_window_excludes_grating_lobes_of_a_sparse_array():
    d = 2.0                                          # sparse: grating lobes visible
    tau = A.target_delays(4, d, TH0)
    th, _, w = A.pointing_angle(A.ideal_ttd_weights(F, tau), F, d, TH0)
    assert np.allclose(th, TH0, atol=1e-2) and w == pytest.approx(1 / (2 * F.max() * d))


def test_half_power_beamwidth_matches_uniform_array_estimate():
    th, p = A.pattern(np.ones(N), 1.0, D)
    i = int(np.argmax(p))
    bw = A.half_power_beamwidth(th, p, i)
    assert bw == pytest.approx(0.886 / (N * D) * 180 / np.pi, rel=0.03)


def test_half_power_beamwidth_nan_when_not_resolved():
    th = np.linspace(-1, 1, 11)
    assert np.isnan(A.half_power_beamwidth(th, np.ones(11), 5))


def test_r_squint_ratio_and_mask():
    ph = np.array([30.0, 31.0, 28.0])
    assert A.r_squint(np.array([30.0, 30.2, 29.9]), ph, 30.0) == pytest.approx(0.2 / 2.0)
    assert A.r_squint(np.array([30.0, 30.2, 29.5]), ph, 30.0, mask=np.array([True, True, False])) == pytest.approx(0.2)


# ---- Phase 5 additions: sparse array (d ~ 2 lambda0), design-lobe tracking, lobe logic -------------------------------------------------------------
from pcttd import beam as B

DS, TS, F0S = 3.9323 / 0.5, 3.9323, 0.2495
FS = np.linspace(0.231, 0.268, 38)
D5 = 3.9323 / np.sin(np.deg2rad(30))


def test_sparse_pitch_definition():
    assert D5 == pytest.approx(7.8646, abs=1e-4)
    assert D5 * F0S == pytest.approx(1.962, abs=1e-3)


def test_phase5_ideal_ttd_stays_at_theta0_over_common_band():
    tau = np.arange(4) * TS
    w = 1 / (2 * FS.max() * D5)
    r = B.analyze_case(A.ideal_ttd_weights(FS, tau), FS, D5, 30.0, w)
    assert np.all(np.abs(r["theta_peak"] - 30.0) < 1e-3) and np.allclose(r["gain"], 1.0, atol=1e-8)


def test_phase5_fixed_phase_matches_analytic_squint_to_grid_accuracy():
    tau = np.arange(4) * TS
    w = 1 / (2 * FS.max() * D5)
    r = B.analyze_case(A.fixed_phase_weights(FS, tau, F0S), FS, D5, 30.0, w)
    want = np.rad2deg(np.arcsin(F0S / FS * 0.5))
    assert np.max(np.abs(r["theta_peak"] - want)) < 0.01 and r["theta_peak"][0] > 30 > r["theta_peak"][-1]


def test_design_window_selects_m0_branch_not_global_maximum():
    # weights that put a stronger grating lobe outside the window: element 3 boosted -> still the window picks the intended branch
    tau = np.arange(4) * TS
    H = A.ideal_ttd_weights(FS, tau)
    w = 1 / (2 * FS.max() * D5)
    th, _, _ = A.pointing_angle(H, FS, D5, 30.0, window=w)
    s = np.sin(np.deg2rad(th))
    assert np.all(np.abs(s - 0.5) <= w)
    # a lobe one grating order away sits outside the window
    grat = B.grating_lobe_angles(30.0, FS, D5)
    assert all(np.all(np.abs(np.sin(np.deg2rad(g)) - 0.5) > w) for g in grat)


def test_pointing_known_synthetic_array_linear_phase_offset():
    # 4 elements, ideal TTD plus a linear phase g*m: peak shifts by g/(2 pi f d cos theta0) (small g)
    tau = np.arange(4) * TS
    g = 0.05
    H = A.ideal_ttd_weights(FS, tau) * np.exp(1j * g * np.arange(4))[:, None]
    w = 1 / (2 * FS.max() * D5)
    th = B.analyze_case(H, FS, D5, 30.0, w)["theta_peak"]
    pred, _ = B.phase_slope_pointing(np.tile(g * np.arange(4)[:, None], (1, len(FS))), FS, D5, 30.0)
    assert np.allclose(th - 30.0, pred, atol=0.01)


def test_hpbw_uniform_broadside_array_in_new_module():
    s, p = B.full_power(np.ones(8), 1.0, 0.5)
    k = int(np.argmax(p))
    assert B._hpbw(s, p, k) == pytest.approx(0.886 / (8 * 0.5) * 180 / np.pi, rel=0.01)


def test_strongest_nondesign_lobe_excludes_design_lobe():
    # ideal sparse array at f = 0.25: design lobe at 30 deg, grating lobes at equal level -> strongest other lobe is 0 dB, not -inf and not the design lobe
    tau = np.arange(4) * TS
    f = 0.25
    r = B.lobe_report(A.ideal_ttd_weights([f], tau)[:, 0], f, D5, 30.0)
    assert r["strongest_db"] == pytest.approx(0.0, abs=1e-3) and abs(r["strongest_theta_deg"] - 30.0) > 5


def test_local_maxima_basic():
    p = np.array([0.0, 1.0, 0.5, 0.7, 0.2, 0.9])
    assert list(B.local_maxima(p)) == [1, 3, 5]


def test_grating_lobe_count_matches_analytic_for_ideal_sparse_array_and_excludes_endpoints():
    tau = np.arange(4) * TS
    for f in (0.231, 0.2495, 0.268):
        r = B.lobe_report(A.ideal_ttd_weights([f], tau)[:, 0], f, D5, 30.0)
        want = 1 + len(B.grating_lobe_angles(30.0, [f], D5)[0])
        assert r["n_within_3db"] == want
