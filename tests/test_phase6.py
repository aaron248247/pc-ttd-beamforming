"""Phase 6: hierarchical subarray model, perturbation primitives, Monte Carlo reproducibility, analytic error link, handoff consistency, Phase 4/5 data integrity."""
import hashlib
from pathlib import Path

import numpy as np
import pytest

from pcttd import array as A, beam as B, robustness as R, subarray as S

ROOT = Path(__file__).resolve().parents[1]
P4 = ROOT / "results" / "phase4" / "phase4_handoff.npz"
F0, TAU = 0.2495, 3.9323
LAM0 = 1 / F0; DE = LAM0 / 2; TH0 = 30.0
FR = 0.231 + 5e-4 * np.arange(75)                         # the Phase 4 grid: contains f0 = 0.2495 exactly

# sha256 of Phase 4/5 result files as committed at the Phase 5 gate (959e4c9); Phase 6 must not modify them
HASH = {"results/phase5/phase5_handoff.npz": "7f0c4bc2cae6150227ced2de75dc515d29f4f3978449bc0697aa516d24c0d23d",
        "results/phase5/phase5_handoff.json": "63810cf0e5a1d54f9e8592e93a5c4e6ef094ecc16069bf6ae8b7aa7d008e8713",
        "results/phase5/5a_array.npz": "996fd5c6ba7f7979162cf046ff796337aedf60682dcbd8f945b257d359dfce8b",
        "results/phase4/phase4_handoff.npz": "20bdd9371a57a0627bc0a85ebd5334772ca16089a8cef06669c9132d68fcbce3",
        "results/phase4/phase4_handoff.json": "2e843854a8f4e3217cd92911a26b49c56be21047c524dee29505d4811e0942df"}


def rand_G(nq=4, seed=1):
    rng = np.random.default_rng(seed)
    return (rng.uniform(0.7, 1.3, (nq, len(FR))) * np.exp(1j * rng.uniform(-3, 3, (nq, len(FR)))))


def test_msub1_reduces_to_ordinary_array():
    G = rand_G()
    h = S.Hier(FR, 4, 1, DE, n=2001)
    U = S.intra_fixed(1, DE, TH0, F0, len(FR))
    p = h.power(G, h.element_factor(U))
    ref = np.abs(A.array_factor(G, FR, DE, h.s)) ** 2 / 16
    assert np.allclose(p, ref, atol=1e-12)


@pytest.mark.parametrize("msub", [2, 4])
def test_hierarchy_equals_element_level_array_factor(msub):
    G = rand_G(seed=msub)
    U = S.intra_ttd(msub, DE, TH0, FR)
    h = S.Hier(FR, 4, msub, DE, n=1501)
    p = h.power(G, h.element_factor(U))
    W = S.hier_weights(G, U)
    assert np.allclose(p, np.abs(A.array_factor(W, FR, DE, h.s)) ** 2 / (4 * msub) ** 2, atol=1e-12)


def test_hier_analyze_matches_beam_analyze_case():
    msub = 2; G = S.fixed_channels(FR, 4, S.tau_required(msub, DE, TH0), F0) * np.exp(1j * 0.1 * np.arange(4))[:, None]
    U = S.intra_fixed(msub, DE, TH0, F0, len(FR)); win = S.window_halfwidth(FR.max(), msub, DE)
    ref = B.analyze_case(S.hier_weights(G, U), FR, DE, TH0, win)
    got = S.Hier(FR, 4, msub, DE).analyze(G, U, TH0, win)
    for k in ("theta_peak", "gain", "hpbw", "strongest_db"):
        assert np.allclose(ref[k], got[k], atol=1e-5)


def test_ideal_ttd_points_to_theta0():
    for msub in (1, 2, 4, 8):
        tr = S.tau_required(msub, DE, TH0); win = S.window_halfwidth(FR.max(), msub, DE)
        h = S.Hier(FR, 4, msub, DE)
        el = h.analyze(S.ideal_channels(FR, 4, tr), S.intra_ttd(msub, DE, TH0, FR), TH0, win)
        assert np.abs(el["theta_peak"] - TH0).max() < 1e-3 and np.allclose(el["gain"], 1.0, atol=1e-6)
        sub = h.analyze(S.ideal_channels(FR, 4, tr), S.intra_fixed(msub, DE, TH0, F0, len(FR)), TH0, win)
        i0 = int(np.argmin(abs(FR - F0)))
        assert abs(sub["theta_peak"][i0] - TH0) < 1e-3                         # exact at f0; intra-subarray fixed phase pulls the lobe slightly elsewhere
        assert np.abs(sub["theta_peak"] - TH0).max() < 0.25                    # ... by less than 0.25 deg for these sizes


def test_half_wave_spacing_has_no_spatial_alias_family():
    gl = B.grating_lobe_angles(np.full(len(FR), TH0), FR, DE)
    assert max(len(g) for g in gl) == 0
    msub = 1; h = S.Hier(FR, 4, msub, DE, n=20001)
    r = h.analyze(S.ideal_channels(FR, 4, S.tau_required(1, DE, TH0)), S.intra_ttd(1, DE, TH0, FR), TH0, 0.9)
    assert np.all(np.isinf(r["grating_db"]) | (r["grating_db"] < -300)) and np.all(r["n_replicas"] == 0)
    assert r["strongest_db"].max() < -10.0                                      # ordinary sidelobes only (no equal-strength alias)
    # the Phase 5 pitch (about 1.96 lambda0) does alias
    assert len(B.grating_lobe_angles([TH0], [F0], TAU / np.sin(np.deg2rad(TH0)))[0]) >= 2


def test_msub4_required_delay_is_one_carrier_cycle():
    assert F0 * S.tau_required(4, DE, TH0) == pytest.approx(1.0, abs=1e-12)
    assert F0 * S.tau_required(2, DE, TH0) == pytest.approx(0.5) and F0 * S.tau_required(8, DE, TH0) == pytest.approx(2.0)
    rho, mm, dc = S.mismatch(TAU, S.tau_required(4, DE, TH0), F0)
    assert rho == pytest.approx(0.98111, abs=1e-4) and mm == pytest.approx(-0.01889, abs=1e-4) and dc == pytest.approx(-0.01889, abs=1e-4)
    assert TAU / (DE * np.sin(np.deg2rad(TH0))) == pytest.approx(3.92444, abs=1e-4)


def test_retrim_changes_only_constant_phase():
    G = rand_G(); tr = S.tau_required(4, DE, TH0)
    Gt, delta = S.retrim(G, FR, tr, F0)
    r = Gt / G
    assert np.allclose(np.abs(r), 1) and np.allclose(r, r[:, [0]])
    i0 = int(np.argmin(abs(FR - F0)))
    rel = np.angle(Gt[:, i0] / Gt[0, i0]); want = np.angle(np.exp(1j * 2 * np.pi * F0 * np.arange(4) * tr))
    assert np.allclose(np.exp(1j * rel), np.exp(1j * want))


def test_phase_perturbation_changes_only_phase():
    G = rand_G(); eps = np.array([0.1, -0.2, 0.3, 0.05])
    Gp = R.perturb(G, FR, eps=eps)
    assert np.allclose(np.abs(Gp), np.abs(G)) and np.allclose(Gp / G, np.exp(1j * eps)[:, None] * np.ones((1, len(FR))))


def test_amplitude_perturbation_changes_only_amplitude():
    G = rand_G(); adb = np.array([0.5, -0.5, 1.0, 0.0])
    Gp = R.perturb(G, FR, amp_db=adb)
    assert np.allclose(np.angle(Gp / G), 0) and np.allclose(20 * np.log10(np.abs(Gp) / np.abs(G)), adb[:, None])


def test_delay_perturbation_phase_slope():
    G = rand_G(); dt = np.array([0.0, 0.02, -0.03, 0.05])
    for fref in (0.0, F0):
        r = np.unwrap(np.angle(R.perturb(G, FR, dtau=dt, f_ref=fref) / G), axis=1)
        assert np.allclose(np.gradient(r, FR, axis=1), 2 * np.pi * dt[:, None])
        assert np.allclose(np.abs(R.perturb(G, FR, dtau=dt, f_ref=fref)), np.abs(G))
    assert np.allclose(np.angle(R.perturb(G, FR, dtau=dt, f_ref=F0)[:, np.argmin(abs(FR - F0))] / G[:, np.argmin(abs(FR - F0))]), 0, atol=1e-12)   # calibrated: no change at f0


def test_zero_perturbation_is_exactly_nominal():
    G = rand_G()
    assert np.array_equal(R.perturb(G, FR), G)
    assert np.array_equal(R.perturb(G, FR, eps=np.zeros(4), dtau=np.zeros(4), amp_db=np.zeros(4)), G)


def _mc(seed, n=6):
    f = FR; msub = 4; tr = S.tau_required(msub, DE, TH0)
    G = S.ideal_channels(f, 4, tr); U = S.intra_fixed(msub, DE, TH0, F0, len(f)); win = S.window_halfwidth(f.max(), msub, DE)
    h = S.Hier(f, 4, msub, DE, n=2001); E = h.element_factor(U); th, g, _ = h.peak(h.power(G, E), TH0, win)
    return R.monte_carlo(h, E, G, f, TH0, win, g, 2.7, np.random.SeedSequence(seed), n, dict(eps_deg=5, delay_frac=0.01, amp_db=0.5), np.arange(4) * TAU)


def test_monte_carlo_reproducible_with_fixed_seed():
    a, b, c = _mc(5), _mc(5), _mc(6)
    for k in a:
        assert np.array_equal(a[k], b[k])
    assert not np.array_equal(a["max_err"], c["max_err"])


def test_first_order_phase_error_prediction_small_error_limit():
    msub = 4; tr = S.tau_required(msub, DE, TH0); win = S.window_halfwidth(FR.max(), msub, DE)
    h = S.Hier(FR, 4, msub, DE); U = S.intra_fixed(msub, DE, TH0, F0, len(FR)); E = h.element_factor(U)
    G = S.ideal_channels(FR, 4, tr); th0 = h.peak(h.power(G, E), TH0, win)[0]
    prev = None
    for a in (0.02, 0.2):
        eps = np.deg2rad(a) * (2 * np.arange(4) - 3) / 3
        sh = h.peak(h.power(R.perturb(G, FR, eps=eps), E), TH0, win)[0] - th0
        pred, _ = B.phase_slope_pointing(np.repeat(eps[:, None], len(FR), 1), FR, msub * DE, TH0)
        assert np.abs(sh - pred).max() < 0.10 * np.abs(pred).max()          # first order; the element factor stiffens the peak by a few percent
        assert np.abs(sh).max() > 0.9 * np.abs(pred).max()


def test_delay_error_prediction_uncalibrated_is_frequency_independent():
    dt = 0.01 * np.arange(4) * TAU
    pred, _ = B.phase_slope_pointing(2 * np.pi * FR[None, :] * dt[:, None], FR, 4 * DE, TH0)
    assert np.ptp(pred) < 1e-9


def test_budget_interpolation():
    assert R.budget([1, 2, 3], [0.2, 0.6, 1.2], 0.4) == pytest.approx(1.5)
    assert R.budget([1, 2], [0.6, 0.7], 0.5) is None and R.budget([1, 2], [0.1, 0.2], 0.5) == float("inf")


@pytest.mark.skipif(not P4.exists(), reason="Phase 4 handoff absent")
def test_res30_res40_handled_consistently():
    z = np.load(P4); mk = z["mask_common"]; f = z["f"][mk]; tr = S.tau_required(4, DE, TH0)
    for sfx in ("", "_res40"):
        H = z["H_raw" + sfx][:, mk]
        Gt, delta = S.retrim(H, f, tr, F0)
        assert Gt.shape == (4, len(f)) and np.allclose(np.abs(Gt), np.abs(H))
        assert np.allclose(np.gradient(np.unwrap(np.angle(Gt), axis=1), f, axis=1), np.gradient(np.unwrap(np.angle(H), axis=1), f, axis=1), atol=1e-8)
    assert z["H_raw"].shape == z["H_raw_res40"].shape and float(z["f0"]) == F0 and float(z["dtau_target"]) == TAU


@pytest.mark.parametrize("rel", sorted(HASH))
def test_phase4_phase5_results_unmodified(rel):
    assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == HASH[rel]
