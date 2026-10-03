import numpy as np
import pytest

from pcttd import tmm as T
from pcttd import extract as X

F = np.linspace(0.05, 0.7, 700)


def test_single_interface_fresnel():
    for n1, n2 in [(1.0, 2.0), (1.5, 1.0), (1.0, 3.4641)]:
        t, r, Tp, Rp = T.transmission_reflection(np.array([0.3]), [], n_in=n1, n_out=n2)
        assert r[0] == pytest.approx((n1 - n2) / (n1 + n2))
        assert t[0] == pytest.approx(2 * n1 / (n1 + n2))
        assert Rp[0] + Tp[0] == pytest.approx(1.0)


def test_vacuum_layer_is_pure_delay_with_positive_sign():
    d = 7.3
    f = F
    t, r, Tp, Rp = T.transmission_reflection(f, [(1.0, d)])
    assert np.allclose(np.abs(t), 1.0) and np.allclose(np.abs(r), 0.0, atol=1e-12)
    assert np.allclose(t, np.exp(1j * 2 * np.pi * f * d))        # arg t = +k d (exp(-i w t) convention)
    tau, ok = X.group_delay_phase(f, t)
    assert np.allclose(tau[2:-2], d, atol=1e-8)
    assert np.allclose(T.tmm_group_delay(f, [(1.0, d)], referenced=False), d, atol=1e-6)


def test_dielectric_layer_in_matching_medium_delay_is_nL():
    n, d = 2.0, 3.0
    # same index everywhere -> no interfaces; delay n d
    t, *_ = T.transmission_reflection(F, [(n, d)], n_in=n, n_out=n)
    tau, _ = X.group_delay_phase(F, t)
    assert np.allclose(tau[2:-2], n * d, atol=1e-8)


def test_referenced_excess_delay_of_vacuum_layer_is_zero():
    assert np.allclose(T.tmm_group_delay(F, [(1.0, 5.0)], referenced=True), 0.0, atol=1e-6)


@pytest.mark.parametrize("seed", range(5))
def test_energy_conservation_random_lossless_multilayer(seed):
    rng = np.random.default_rng(seed)
    layers = [(float(rng.uniform(1, 4)), float(rng.uniform(0.05, 1.5))) for _ in range(rng.integers(1, 12))]
    n_in, n_out = float(rng.uniform(1, 2)), float(rng.uniform(1, 3))
    t, r, Tp, Rp = T.transmission_reflection(F, layers, n_in, n_out)
    assert np.max(np.abs(Tp + Rp - 1.0)) < 1e-12


def test_reciprocity_of_transmission_under_layer_reversal():
    layers = [(2.0, 0.3), (3.0, 0.5), (1.5, 0.2)]
    t1, *_ = T.transmission_reflection(F, layers, 1.0, 1.0)
    t2, *_ = T.transmission_reflection(F, layers[::-1], 1.0, 1.0)
    assert np.allclose(t1, t2, atol=1e-12)


def test_half_wave_layer_is_transparent():
    n, f0 = 2.5, 0.2
    d = 0.5 / (n * f0)                      # half-wave optical thickness at f0
    _, _, Tp, _ = T.transmission_reflection(np.array([f0]), [(n, d)])
    assert Tp[0] == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("N", [1, 2, 4, 8])
def test_quarter_wave_peak_reflectance_closed_form(N):
    n_lo, n_hi = 1.0, np.sqrt(12.0)
    d_lo, d_hi, f0 = T.quarter_wave_bilayer(n_lo, n_hi)
    layers = T.bilayer(n_lo, d_lo, n_hi, d_hi, N)
    _, _, Tp, Rp = T.transmission_reflection(np.array([f0]), layers, 1.0, 1.0)
    assert Rp[0] == pytest.approx(T.quarter_wave_peak_reflectance(n_lo, n_hi, N), rel=1e-10)
    assert Rp[0] + Tp[0] == pytest.approx(1.0, abs=1e-12)


def test_quarter_wave_gap_edges_match_bloch_condition():
    n_lo, n_hi = 1.0, np.sqrt(12.0)
    d_lo, d_hi, f0 = T.quarter_wave_bilayer(n_lo, n_hi)
    lo, hi = T.quarter_wave_gap_edges(n_lo, n_hi)
    assert lo < f0 < hi
    for edge in (lo, hi):
        assert abs(T.bloch_cos(edge, n_lo, d_lo, n_hi, d_hi)) == pytest.approx(1.0, abs=1e-9)
    # |cos| < 1 (pass band) just outside, > 1 (stop band) just inside
    assert abs(T.bloch_cos(lo * 0.98, n_lo, d_lo, n_hi, d_hi)) < 1
    assert abs(T.bloch_cos(lo * 1.02, n_lo, d_lo, n_hi, d_hi)) > 1
    assert abs(T.bloch_cos(hi * 1.02, n_lo, d_lo, n_hi, d_hi)) < 1
    assert abs(T.bloch_cos(hi * 0.98, n_lo, d_lo, n_hi, d_hi)) > 1


def test_finite_stack_transmission_follows_bloch_decay_in_stop_band():
    n_lo, n_hi = 1.0, np.sqrt(12.0)
    d_lo, d_hi, f0 = T.quarter_wave_bilayer(n_lo, n_hi)
    f = np.array([f0 * 0.9, f0, f0 * 1.1])
    rate = T.decay_rate_per_period(f, n_lo, d_lo, n_hi, d_hi, 1.0)
    assert np.all(np.isfinite(rate))
    lnT = []
    for N in (6, 7, 8, 9, 10):
        _, _, Tp, _ = T.transmission_reflection(f, T.bilayer(n_lo, d_lo, n_hi, d_hi, N))
        lnT.append(np.log(Tp))
    slope = np.polyfit([6, 7, 8, 9, 10], np.array(lnT), 1)[0]
    assert np.allclose(-slope, rate, rtol=1e-3)


def test_delay_in_stop_band_saturates_with_length():
    # Evanescent tunneling (Hartman effect): mid-gap group delay converges to a length-independent value, while the
    # free-propagation delay N*(optical period) keeps growing. Values from the exact TMM: tau(8) and tau(10) agree to < 1e-6,
    # tau(10) is ~1/11 of the free-propagation delay.
    n_lo, n_hi = 1.0, np.sqrt(12.0)
    d_lo, d_hi, f0 = T.quarter_wave_bilayer(n_lo, n_hi)
    tau = {N: T.tmm_group_delay(np.array([f0]), T.bilayer(n_lo, d_lo, n_hi, d_hi, N), referenced=False)[0]
           for N in (4, 6, 8, 10)}
    assert abs(tau[10] - tau[8]) < 1e-6
    free10 = 10 * (n_lo * d_lo + n_hi * d_hi)
    assert tau[10] < 0.1 * free10
