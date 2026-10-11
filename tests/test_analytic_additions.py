"""Analytic additions (steering-angle sweep, match table, conventional propagation reference): consistency with Phase 4-6 results and internal arithmetic."""
import json
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
AS = json.load(open(R / "final" / "angle_sweep.json"))
CR = json.load(open(R / "final" / "conventional_reference.json"))
H6 = json.load(open(R / "phase6" / "phase6_handoff.json"))
SW = {k: np.array(v) for k, v in AS["sweep"].items()}
I30 = int(np.argmin(abs(SW["theta0"] - 30.0)))


def test_sweep_at_30_reproduces_phase6():
    m = H6["metrics_res30"]["4"]
    assert SW["theta0"][I30] == pytest.approx(30.0)
    assert SW["max_err_trim"][I30] == pytest.approx(m["max_abs_err"]["trim"], rel=1e-9)
    assert SW["max_err_fixed"][I30] == pytest.approx(m["max_abs_err"]["fixed"], rel=1e-9)
    assert SW["max_err_ideal_sub"][I30] == pytest.approx(m["max_abs_err"]["ideal_sub"], rel=1e-9)   # ideal-subarray convention
    assert SW["R_squint"][I30] == pytest.approx(m["R_squint"]["trim"], rel=1e-9)
    assert SW["vs_ideal_sub_max"][I30] == pytest.approx(m["vs_ideal_sub_max"]["trim"], rel=1e-9)
    assert SW["strongest_db_max"][I30] == pytest.approx(m["strongest_db_max"]["trim"], abs=1e-9)
    assert SW["grating_db_max"][I30] == pytest.approx(m["grating_db_max"]["trim"], abs=1e-9)
    assert SW["gain_min"][I30] == pytest.approx(m["gain_min"]["trim"], rel=1e-9)
    assert SW["rho"][I30] == pytest.approx(H6["architectures"]["4"]["rho"], rel=1e-12)


def test_exact_match_angle_is_numerical_minimum_mismatch():
    s = AS["summary"]
    q = s["q_phc"]
    assert s["theta_match_exact_deg"] == pytest.approx(np.rad2deg(np.arcsin(2 * q / 4)), abs=1e-9)
    i = int(np.argmin(abs(SW["theta0"] - s["theta_match_exact_deg"])))
    assert abs(SW["mismatch_pct"][i]) < 1e-6 and SW["rho"][i] == pytest.approx(1.0, abs=1e-6)
    assert s["theta_min_abs_mismatch_grid"] == pytest.approx(s["theta_match_exact_deg"], abs=1e-4)
    assert s["theta_min_vs_ideal_sub_grid"] == pytest.approx(s["theta_match_exact_deg"], abs=1e-4)        # deviation from ideal subarray TTD is smallest at the match
    assert abs(s["theta_min_max_err_grid"] - s["theta_match_exact_deg"]) < 3.0                           # absolute error minimum within ~1 intra-subarray offset of it


def test_match_table_arithmetic():
    q = AS["summary"]["q_phc"]
    for t in AS["match_table"]:
        for e in t["matches"]:
            assert t["msub"] * e["sin_theta"] / 2 == pytest.approx(e["k"] * q, rel=1e-12)
            assert e["sin_theta"] <= 1 + 1e-12
        assert 2 * (t["k_max"] + 1) * q / t["msub"] > 1 or t["msub"] == 0          # the next k would not be visible
    assert [t["k_max"] for t in AS["match_table"]] == [0, 1, 2, 4, 8]


def test_conventional_reference_consistency():
    r = CR["references"]
    f0, tau = CR["f0"], CR["tau_phc_a_over_c"]
    for k in r:
        assert r[k]["L_over_lambda0"] * r[k]["n_g"] == pytest.approx(1.0, rel=1e-12)             # exactly one carrier cycle: L n_g / c = 1/f0
        assert r[k]["L_over_a"] == pytest.approx(r[k]["L_over_lambda0"] / f0, rel=1e-12)
    assert r["air"]["L_over_lambda0"] == pytest.approx(1.0)
    assert r["uniform_eps12"]["n_g"] == pytest.approx(np.sqrt(12))
    assert r["W1_PhC"]["L_over_a"] == pytest.approx(1 / (f0 * tau), rel=1e-12)                    # period route
    assert r["W1_PhC"]["delay_per_period_cycles"] == pytest.approx(AS["summary"]["q_phc"], rel=1e-12)
    assert r["W1_PhC"]["relative_length_vs_W1"] == pytest.approx(1.0)
    assert r["uniform_eps12"]["relative_length_vs_W1"] == pytest.approx(tau / np.sqrt(12), rel=1e-12)


def test_phase4_5_6_results_unchanged():
    try:
        out = subprocess.check_output(["git", "status", "--porcelain", "--", "results/phase4", "results/phase5", "results/phase6"], cwd=ROOT, text=True)
    except Exception:
        pytest.skip("git unavailable")
    assert out.strip() == ""
