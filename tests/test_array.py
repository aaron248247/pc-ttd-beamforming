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
