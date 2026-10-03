import numpy as np
import pytest

from pcttd import edges as E
from pcttd import tmm as T

N_LO, N_HI = 1.0, float(np.sqrt(12.0))
D_LO, D_HI, F0 = T.quarter_wave_bilayer(N_LO, N_HI)
LO, HI = T.quarter_wave_gap_edges(N_LO, N_HI)
F = np.round(np.arange(0.12, 0.52, 0.0005), 6)


def lnT_for(Ns):
    return np.array([np.log(T.transmission_reflection(F, T.bilayer(N_LO, D_LO, N_HI, D_HI, N))[2]) for N in Ns])


def test_attenuation_per_period_matches_bloch_rate_in_gap():
    Ns = [6, 8, 10, 12]
    kappa, rms = E.attenuation_per_period(Ns, lnT_for(Ns))
    f_test = np.array([0.23, 0.30, 0.38])
    rate = T.decay_rate_per_period(f_test, N_LO, D_LO, N_HI, D_HI, 1.0)
    got = np.interp(f_test, F, kappa)
    assert np.allclose(got, rate, rtol=2e-3)


def test_pass_band_has_no_systematic_decay():
    Ns = [6, 8, 10, 12]
    kappa, rms = E.attenuation_per_period(Ns, lnT_for(Ns))
    m = (F > 0.46) & (F < 0.52)
    assert np.median(np.abs(kappa[m])) < 0.06           # measured 0.039 (Fabry-Perot fringes make single-frequency slopes scatter)
    gap = (F > 0.25) & (F < 0.40)
    assert np.min(kappa[gap]) > 0.5                      # versus > 0.5 per period everywhere inside the gap


@pytest.mark.parametrize("Ns,tol", [([6, 8, 10, 12], 6e-4), ([8, 12, 16, 20], 3e-4), ([4, 6, 8, 10], 1.5e-3)])
def test_kappa2_extrapolation_recovers_known_bragg_edges(Ns, tol):
    kappa, rms = E.attenuation_per_period(Ns, lnT_for(Ns))
    lo, ilo = E.edge_from_kappa2(F, kappa, rms, F0, "lo")
    hi, ihi = E.edge_from_kappa2(F, kappa, rms, F0, "hi")
    # bias of the estimator vs the analytic infinite-lattice edges: measured 3.8e-4/3.9e-4 (N 6-12), 1.4e-4/1.3e-4 (N 8-20),
    # 9.3e-4/8.8e-4 (N 4-10); the tolerances are ~1.5x those values. Shorter sections -> larger bias.
    assert abs(lo - LO) < tol and abs(hi - HI) < tol
    assert ilo["n_points"] >= 6 and ihi["n_points"] >= 6


def test_estimator_is_insensitive_to_a_length_independent_prefactor():
    # multiplying T by any frequency-dependent, N-independent factor (interface transmission) must not move the edges
    Ns = [6, 8, 10, 12]
    L = lnT_for(Ns) + np.log(0.4 + 0.3 * np.sin(40 * F))[None, :]
    k0, r0 = E.attenuation_per_period(Ns, lnT_for(Ns))
    k1, r1 = E.attenuation_per_period(Ns, L)
    assert np.allclose(k0, k1, atol=1e-9) and np.allclose(r0, r1, atol=1e-9)


def test_too_few_points_returns_nan():
    kappa = np.full_like(F, 5.0)
    e, info = E.edge_from_kappa2(F, kappa, np.zeros_like(F), F0, "lo")
    assert np.isnan(e) and info["n_points"] == 0
