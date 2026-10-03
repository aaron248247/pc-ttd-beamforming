import numpy as np
import pytest

from pcttd import extract as X


F = np.linspace(0.15, 0.45, 301)          # df = 0.001 -> unambiguous delay < 500


@pytest.mark.parametrize("tau", [0.0, 1.0, 10.0, 123.4])
def test_group_delay_of_pure_delay_recovers_tau_with_meep_sign(tau):
    # MEEP convention: H = exp(+i omega tau) for a causal delay tau
    H = np.exp(1j * X.omega(F) * tau)
    t_phase, ok = X.group_delay_phase(F, H, sign=+1)
    assert ok
    assert np.allclose(t_phase[2:-2], tau, atol=1e-6)
    assert np.allclose(X.group_delay_ratio(F, H, sign=+1), tau, atol=1e-6)


def test_wrong_sign_gives_negative_delay_and_is_detected_by_sign():
    H = np.exp(1j * X.omega(F) * 10.0)
    assert np.all(X.group_delay_phase(F, H, sign=-1)[0][2:-2] < 0)


def test_opposite_time_convention_needs_opposite_sign():
    # exp(+i omega t) convention: causal delay appears as H = exp(-i omega tau); then sign=-1 recovers +tau
    H = np.exp(-1j * X.omega(F) * 7.0)
    assert np.allclose(X.group_delay_phase(F, H, sign=-1)[0][2:-2], 7.0, atol=1e-6)


def test_unwrapping_handles_many_wraps_and_amplitude_variation():
    tau = 200.0                                    # ~60 full 2*pi wraps over the band
    amp = 1.0 + 0.5 * np.sin(40 * F)
    H = amp * np.exp(1j * (X.omega(F) * tau + 0.3))
    t, ok = X.group_delay_phase(F, H)
    assert ok and np.allclose(t[2:-2], tau, atol=1e-6)


def test_sampling_guard_rejects_grids_too_coarse_for_the_expected_delay():
    coarse = np.linspace(0.15, 0.45, 31)           # df = 0.01 -> unambiguous delay < 50
    assert X.max_unambiguous_delay(0.01) == pytest.approx(50.0)
    with pytest.raises(ValueError):
        X.assert_sampling_adequate(coarse, tau_max=200.0)
    assert X.assert_sampling_adequate(coarse, tau_max=30.0)
    assert X.assert_sampling_adequate(F, tau_max=300.0)


def test_aliasing_is_silent_without_the_guard_documented_limitation():
    # tau = 200 with df = 0.01 is a phase step of exactly 4 pi per sample: it aliases to "no phase change".
    # The wrapped-step check cannot see this, which is why the a-priori guard above is mandatory.
    coarse = np.linspace(0.15, 0.45, 31)
    H = np.exp(1j * X.omega(coarse) * 200.0)
    t, ok = X.group_delay_phase(coarse, H)
    assert ok                                       # not flagged
    assert not np.allclose(t[2:-2], 200.0, atol=1.0)   # ... yet wrong
    assert np.allclose(t[2:-2], 0.0, atol=1e-6)         # it reads as zero delay


def test_unwrapping_inside_the_limit():
    df = 0.001
    tau = 0.85 * X.max_unambiguous_delay(df)        # 425 < 500; phase step 0.85*pi
    f = np.arange(0.2, 0.3, df)
    H = np.exp(1j * X.omega(f) * tau)
    t, ok = X.group_delay_phase(f, H)
    assert ok and np.allclose(t[2:-2], tau, atol=1e-5)


def test_ratio_method_has_half_the_delay_limit_of_unwrapping():
    df = 0.001
    f = np.arange(0.2, 0.3, df)
    ok_tau = 0.85 * X.max_unambiguous_delay(df) / 2           # 212.5: inside the 2-sample-stencil limit
    assert np.allclose(X.group_delay_ratio(f, np.exp(1j * X.omega(f) * ok_tau))[1:-1], ok_tau, atol=1e-5)
    bad_tau = 0.85 * X.max_unambiguous_delay(df)              # 425: fine for unwrapping, aliased for the ratio method
    assert not np.allclose(X.group_delay_ratio(f, np.exp(1j * X.omega(f) * bad_tau))[1:-1], bad_tau, rtol=0.1)
    assert X.assert_sampling_adequate(f, 0.8 * ok_tau, span=2)
    with pytest.raises(ValueError):
        X.assert_sampling_adequate(f, bad_tau, span=2)


def test_past_the_limit_the_unwrapped_result_is_wrong():
    df = 0.001
    tau = 1.3 * X.max_unambiguous_delay(df)         # step 1.3 pi -> aliases to -0.7 pi
    f = np.arange(0.2, 0.3, df)
    H = np.exp(1j * X.omega(f) * tau)
    t, _ = X.group_delay_phase(f, H)
    assert not np.allclose(t[2:-2], tau, rtol=0.1)
    with pytest.raises(ValueError):
        X.assert_sampling_adequate(f, tau)


def test_masking_low_amplitude_samples_does_not_corrupt_neighbors():
    tau = 25.0
    H = np.exp(1j * X.omega(F) * tau)
    H[100:110] = 1e-12 * np.exp(1j * 5.0)            # deep null with garbage phase
    t, ok = X.group_delay_phase(F, H, mask_below=1e-6)
    good = np.isfinite(t)
    assert not good[100:110].any()
    assert np.allclose(t[good][5:-5][np.abs(F[good][5:-5] - F[105]) > 0.01], tau, atol=1e-4)


def test_phase_and_complex_derivative_agree_for_dispersive_response():
    # smooth dispersive H: phase = a*(w-w0)^2/2 + tau0*w -> tau(w) = tau0 + a*(w-w0)
    w = X.omega(F)
    w0, a, tau0 = X.omega(0.3), 150.0, 30.0       # |tau| <= ~171 < 250 (ratio-method limit at df = 0.001)
    H = np.exp(1j * (0.5 * a * (w - w0) ** 2 + tau0 * w))
    exact = tau0 + a * (w - w0)
    for t in (X.group_delay_phase(F, H)[0], X.group_delay_ratio(F, H)):
        assert np.allclose(t[2:-2], exact[2:-2], atol=1e-6)


@pytest.mark.parametrize("tau", [5.0, 20.0, 60.0])
def test_derivative_method_error_follows_its_predicted_scaling(tau):
    # truncation error of the central-difference derivative of exp(i w tau): tau*(dw*tau)^2/6
    H = np.exp(1j * X.omega(F) * tau)
    dw = X.omega(F[1]) - X.omega(F[0])
    err = abs(X.group_delay_derivative(F, H)[150] - tau)
    pred = tau * (dw * tau) ** 2 / 6.0
    assert err == pytest.approx(pred, rel=0.05)


def test_delay_vs_length_fit():
    L = np.array([5.0, 10.0, 20.0, 40.0])
    slope, icpt, res = X.fit_delay_vs_length(L, 2.0 * L + 0.01)
    assert slope == pytest.approx(2.0) and icpt == pytest.approx(0.01) and res < 1e-12


def test_observed_order_recovers_known_power_law():
    h = np.array([1 / 20, 1 / 40, 1 / 80, 1 / 160])
    for p in (1.0, 2.0, 3.0):
        order, _ = X.observed_order(h, 7.0 * h ** p)
        assert order == pytest.approx(p, abs=1e-9)
    assert np.isnan(X.observed_order([0.1], [1.0])[0])


def test_yee_prediction_reduces_to_ideal_delay_at_fine_resolution():
    f = np.array([0.1, 0.3])
    for n in (1.0, 2.0):
        dx = 1 / 2000
        ph, tau = X.yee_1d_delay_phase(f, 10.0, n, dx, 0.5 * dx / 1.0)
        assert np.allclose(tau, 10.0 * n, rtol=1e-6)
        assert np.allclose(ph, X.omega(f) * 10.0 * n, rtol=1e-6)


def test_yee_prediction_second_order_in_dx():
    f = np.array([0.3])
    errs = []
    hs = np.array([1 / 20, 1 / 40, 1 / 80, 1 / 160])
    for dx in hs:
        _, tau = X.yee_1d_delay_phase(f, 10.0, 1.0, dx, 0.5 * dx)
        errs.append(abs(tau[0] - 10.0))
    order, _ = X.observed_order(hs, errs)
    assert order == pytest.approx(2.0, abs=0.05)


@pytest.mark.parametrize("n,S", [(1.0, 0.5), (2.0, 0.5), (3.4641, 0.5), (1.0, 0.25)])
def test_yee_group_delay_second_order_closed_form(n, S):
    # derived in docs/phase1_methods.md: k = n w [1 + (w dx)^2 (n^2 - S^2)/24], so tau = n L [1 + (w dx)^2 (n^2 - S^2)/8],
    # with S = dt/dx. The relative remainder must vanish like (w dx)^4.
    f = np.array([0.2, 0.3])
    w = X.omega(f)
    L = 10.0
    rem = []
    for dx in (1 / 400, 1 / 800):
        _, tau = X.yee_1d_delay_phase(f, L, n, dx, S * dx)
        approx = n * L * (1 + (w * dx) ** 2 * (n * n - S * S) / 8.0)
        rem.append(np.max(np.abs(tau - approx) / (tau - n * L)))
    assert rem[0] < 0.02 * (n * n) and rem[1] < rem[0] / 3.5        # relative error of the leading term shrinks ~ (w dx)^2
