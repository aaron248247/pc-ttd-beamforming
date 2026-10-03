"""MEEP-based verification with analytically known answers (Phase 1A-1C). Marked ``meep``.

Tolerances come from the measured numerical uncertainty recorded in results/phase1 (documented next to each assertion),
not from convenience: where a quantity is known to carry second-order discretization error the test checks that scaling.
"""
from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("meep")
from pcttd import extract as X  # noqa: E402
from pcttd import fdtd1d as F1  # noqa: E402
from pcttd import tmm as T  # noqa: E402

pytestmark = pytest.mark.meep
BAND = (0.16, 0.44)
BASE = F1.Cfg1D(res=20)


def H_of(r, i=1, j=0):
    return r["E"][i] / r["E"][j]


def test_dft_frequency_axis_is_the_requested_one():
    r = F1.run_1d(BASE)
    assert np.array_equal(r["f"], np.array(BASE.freqs))      # regression: (fcen, df, nfreq) vs (fmin, fmax, nfreq) mix-up


def test_vacuum_delay_has_positive_sign_and_equals_L_over_c_up_to_yee_dispersion():
    r = F1.run_1d(BASE)
    H = H_of(r)
    m = X.band_mask(r["f"], *BAND)
    tp, ok = X.group_delay_phase(r["f"], H, sign=+1)
    tn, _ = X.group_delay_phase(r["f"], H, sign=-1)
    assert ok
    assert np.all(tp[m] > 0) and np.all(tn[m] < 0)               # sign established numerically: +1 is right
    assert X.MEEP_SIGN == +1
    assert np.allclose(np.abs(H[m]), 1.0, atol=1e-4)
    assert abs(np.mean(tp[m]) - 10.0) < 0.015                    # res 20 measured error 0.0090 (phase-velocity Yee dispersion)
    _, pred = X.yee_1d_delay_phase(r["f"], 10.0, 1.0, r["dx"], r["dt"])
    assert np.max(np.abs(tp[m] - pred[m])) < 1e-3               # measured 3.1e-4: extraction == exact Yee prediction
    assert np.max(np.abs(tp[m] - X.group_delay_ratio(r["f"], H)[m])) < 1e-6   # two independent estimators agree


def test_vacuum_error_is_second_order_in_resolution():
    errs = []
    for res in (20, 40):
        r = F1.run_1d(replace(BASE, res=res))
        m = X.band_mask(r["f"], *BAND)
        errs.append(abs(np.mean(X.group_delay_phase(r["f"], H_of(r))[0][m]) - 10.0))
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.1)       # measured order 2.00 over res 10..160


@pytest.mark.parametrize("n", [1.5, 2.0])
def test_homogeneous_dielectric_delay_is_nL_and_linear_in_length(n):
    Ls = [5.0, 10.0, 20.0]
    taus, preds = [], []
    for L in Ls:
        r = F1.run_1d(replace(BASE, eps_bg=n * n, monitors=(0.0, L), z_hi=L + 8.0))
        m = X.band_mask(r["f"], *BAND)
        taus.append(np.mean(X.group_delay_phase(r["f"], H_of(r))[0][m]))
        preds.append(np.mean(X.yee_1d_delay_phase(r["f"], L, n, r["dx"], r["dt"])[1][m]))
    slope, icpt, resid = X.fit_delay_vs_length(Ls, taus)
    pslope, _, _ = X.fit_delay_vs_length(Ls, preds)
    assert slope == pytest.approx(pslope, rel=2e-4)               # extraction reproduces the exact Yee-grid prediction
    assert slope == pytest.approx(n, rel=6e-3)                    # vs ideal n: measured slope/n-1 = 2.4e-3 (n=1.5), 4.5e-3 (n=2) at res 20
    assert abs(icpt) < 1e-3 and resid < 1e-3                      # delay is proportional to length (no offset)


def test_far_monitor_is_not_cut_short_by_decay_criterion():
    # regression: stop_when_dft_decayed fired before the pulse reached a monitor 45 units away (|H| came out 0 or noise)
    L = 40.0
    r = F1.run_1d(replace(BASE, monitors=(0.0, L), z_hi=L + 8.0))
    m = X.band_mask(r["f"], *BAND)
    assert r["t_end"] > r["t_transit"]
    assert np.allclose(np.abs(H_of(r))[m], 1.0, atol=1e-4)
    assert np.mean(X.group_delay_phase(r["f"], H_of(r))[0][m]) == pytest.approx(L, abs=0.05)


def test_dielectric_slab_matches_tmm_and_conserves_energy():
    n, Ls, res = 2.0, 5.0, 40
    zl = 5.0 - Ls / 2
    cfg = replace(BASE, res=res, decay_tol=1e-12, blocks=((zl, zl + Ls, n * n),))
    rs, rr = F1.run_1d(cfg), F1.run_1d(replace(cfg, blocks=()))
    f = rs["f"]
    t = rs["E"][1] / rr["E"][1]
    r = (rs["E"][0] - rr["E"][0]) / rr["E"][0]
    m = X.band_mask(f, *BAND)
    t_tmm = T.vacuum_referenced_t(f, [(n, Ls)])
    _, _, Tt, Rt = T.transmission_reflection(f, [(n, Ls)])
    assert np.max(np.abs(np.abs(r) ** 2 + np.abs(t) ** 2 - 1)[m]) < 1e-4       # measured 2.6e-5
    assert np.max(np.abs(np.abs(t) ** 2 - Tt)[m]) < 2e-2                       # measured 8.6e-3 at res 40, second order
    assert np.max(np.abs(np.abs(r) ** 2 - Rt)[m]) < 2e-2
    tau_m, _ = X.group_delay_phase(f, t)
    tau_t, _ = X.group_delay_phase(f, t_tmm)                                   # same estimator on both
    assert np.max(np.abs(tau_m - tau_t)[m]) < 0.3                              # measured 0.117 (of 7.5 max) at res 40


def test_bragg_stack_stopband_and_bloch_decay_match_tmm():
    n_lo, n_hi = 1.0, float(np.sqrt(12.0))
    d_lo, d_hi, f0 = T.quarter_wave_bilayer(n_lo, n_hi)
    N = 4
    freqs = tuple(np.round(np.arange(0.10, 0.60 + 1e-9, 0.002), 6))
    blocks = tuple((p + d_lo, p + 1.0, n_hi ** 2) for p in range(N))
    cfg = F1.Cfg1D(res=60, z_hi=N + 8.0, monitors=(-2.0, N + 2.0), freqs=freqs, src_fcen=0.35, src_fwidth=0.45,
                   decay_tol=1e-12, blocks=blocks)
    rs, rr = F1.run_1d(cfg), F1.run_1d(replace(cfg, blocks=()))
    f = rs["f"]
    Tm = np.abs(rs["E"][1] / rr["E"][1]) ** 2
    Rm = np.abs((rs["E"][0] - rr["E"][0]) / rr["E"][0]) ** 2
    layers = T.bilayer(n_lo, d_lo, n_hi, d_hi, N)
    _, _, Tt, Rt = T.transmission_reflection(f, layers)
    assert np.max(np.abs(Tm + Rm - 1)[X.band_mask(f, 0.12, 0.58)]) < 1e-3
    # stop band: T < 1e-2 over the analytic gap interior, and rms |T_meep - T_tmm| small
    lo, hi = T.quarter_wave_gap_edges(n_lo, n_hi)
    inside = (f > lo * 1.1) & (f < hi * 0.9)
    assert Tm[inside].max() < 1e-2 and Tt[inside].max() < 1e-2
    assert np.sqrt(np.mean((Tm - Tt)[X.band_mask(f, 0.12, 0.58)] ** 2)) < 2e-2    # measured ~5e-3 for N=4/N=8 at res 40-80
    # mid-gap attenuation vs the infinite-lattice Bloch rate: |T_meep(N=4)| ~ exp(-rate*N) within the finite-size prefactor
    j = int(np.argmin(np.abs(f - 0.25)))
    rate = float(T.decay_rate_per_period(np.array([0.25]), n_lo, d_lo, n_hi, d_hi, 1.0)[0])
    assert Tt[j] > 0 and abs(np.log(Tm[j] / Tt[j])) < 0.2           # same finite-N value to ~20% in T at the test resolution
    assert rate == pytest.approx(1.9857, abs=1e-3)
