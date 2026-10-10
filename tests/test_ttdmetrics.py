import numpy as np
import pytest

from pcttd import ttdmetrics as TM


def test_fractional_bandwidth_and_errors():
    assert TM.fractional_bandwidth(0.9, 1.1) == pytest.approx(0.2)
    assert TM.fractional_bandwidth(1.0, 1.0) == 0.0
    with pytest.raises(ValueError):
        TM.fractional_bandwidth(1.1, 0.9)


def test_gdd_of_quadratic_delay_and_units():
    f = np.linspace(0.2, 0.3, 201)
    w = 2 * np.pi * f
    tau = 4.0 + 3.0 * (w - w[100])            # tau linear in omega -> GDD = 3 (a/c)^2 everywhere
    assert np.allclose(TM.gdd(f, tau), 3.0)
    tau2 = 1.0 + 0.5 * w ** 2                # d tau / d omega = w (exact for central differences interior)
    assert np.allclose(TM.gdd(f, tau2)[1:-1], w[1:-1])
    assert np.allclose(TM.gdd(f, np.full_like(f, 5.0)), 0.0)


def test_median_ripple_and_deviation():
    tau = np.array([3.9, 4.0, 4.1, 4.0, 4.0])
    assert TM.median_delay(tau) == 4.0
    assert TM.delay_ripple(tau) == pytest.approx(0.2 / 4.0)
    assert TM.max_relative_deviation(tau) == pytest.approx(0.1 / 4.0)
    assert TM.delay_ripple(tau, mask=[0, 1, 0, 1, 1]) == 0.0
    assert TM.delay_ripple(np.array([2.0, np.nan, 2.0])) == 0.0


def test_contiguous_runs():
    assert TM.contiguous_runs([0, 1, 1, 0, 1, 0, 1, 1]) == [(1, 2), (4, 4), (6, 7)]
    assert TM.contiguous_runs([0, 0]) == []
    assert TM.contiguous_runs([1, 1]) == [(0, 1)]


def test_usable_band_picks_widest_flat_window():
    f = np.linspace(0.20, 0.30, 101)
    tau = np.where(f < 0.24, 4.0 + 40 * (0.24 - f), 4.0)      # steep rise below 0.24, flat above
    ub = TM.usable_band(f, tau, np.ones_like(f, bool), tol=0.10)
    assert ub["f_high"] == pytest.approx(0.30)
    assert ub["f_low"] == pytest.approx(0.23, abs=1e-9)   # ramp reaches +10% of the (flat) median 4.0 at f = 0.23
    assert ub["max_dev"] <= 0.10 and ub["ripple"] <= 0.2
    ub5 = TM.usable_band(f, tau, np.ones_like(f, bool), tol=0.05)
    assert ub5["frac_bw"] <= ub["frac_bw"] + 1e-12
    ub25 = TM.usable_band(f, tau, np.ones_like(f, bool), tol=0.25)
    assert ub25["frac_bw"] >= ub["frac_bw"] - 1e-12


def test_usable_band_respects_mask_and_contiguity():
    f = np.linspace(0.20, 0.30, 101)
    tau = np.full_like(f, 4.0)
    ok = np.ones_like(f, bool); ok[40:45] = False           # a transmission dip splits the flat region
    ub = TM.usable_band(f, tau, ok)
    assert (ub["i_low"], ub["i_high"]) == (45, 100)          # the longer side; never spans the gap
    ok2 = np.zeros_like(f, bool); ok2[10] = True               # single sample: below min_samples
    assert TM.usable_band(f, tau, ok2) is None
    tau_nan = tau.copy(); tau_nan[50:] = np.nan
    ub = TM.usable_band(f, tau_nan, np.ones_like(f, bool))
    assert ub["i_high"] == 49 and ub["i_low"] == 0


def test_usable_band_none_when_nothing_flat():
    f = np.linspace(0.2, 0.3, 11)
    tau = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024.0])
    assert TM.usable_band(f, tau, np.ones(11, bool), tol=0.10, min_samples=2) is None
    assert TM.usable_band(f, tau, np.zeros(11, bool)) is None


def test_usable_band_median_is_interval_specific():
    f = np.linspace(0.2, 0.3, 1001)
    tau = np.linspace(3.5, 4.5, 1001)
    ub = TM.usable_band(f, tau, np.ones_like(f, bool), tol=0.10, min_samples=2)
    assert ub["max_dev"] <= 0.10 + 1e-12
    span = tau[ub["i_high"]] - tau[ub["i_low"]]
    assert span == pytest.approx(0.2 * ub["median_tau"], rel=0.02)     # a ramp centered on its median: total span = 2 tol m


def test_gbp_length_target_and_common_band():
    assert TM.delay_bandwidth_product(0.15, 4.0) == pytest.approx(0.6)
    assert TM.required_length_periods(100.0, 4.0) == pytest.approx(25.0)
    with pytest.raises(ValueError):
        TM.required_length_periods(100.0, 0.0)
    assert TM.target_delay_from_cycles(25, 0.25) == pytest.approx(100.0)
    assert TM.common_band([(0.23, 0.26), (0.24, 0.27)]) == (0.24, 0.26)
    assert TM.common_band([(0.23, 0.24), (0.25, 0.27)]) is None
    assert TM.common_band([(0.23, 0.24), None]) is None


def test_interp_on_grid_has_no_extrapolation():
    y = TM.interp_on_grid([0.5, 1.5, 3.0], [0, 1, 2], [0, 10, 20])
    assert y[0] == 5 and y[1] == 15 and np.isnan(y[2])
