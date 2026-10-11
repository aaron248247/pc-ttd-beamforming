"""Phase 7: final results/handoff consistency, figure/table/document presence, Phase 4-6 data integrity, README script references."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
FR = json.load(open(R / "final" / "final_results.json"))
HO = json.load(open(R / "final" / "final_handoff.json"))
V = lambda k: FR[k]["value"]


def test_final_results_match_phase6_handoff():
    h6 = json.load(open(R / "phase6" / "phase6_handoff.json")); m = h6["metrics_res30"]["4"]
    assert V("I_msub4_trim_max_pointing_error") == pytest.approx(m["max_abs_err"]["trim"], rel=1e-12)
    assert V("J_msub4_trim_R_squint") == pytest.approx(m["R_squint"]["trim"], rel=1e-12)
    assert V("L_msub4_trim_strongest_replica_lobe_db") == pytest.approx(m["grating_db_max"]["trim"], rel=1e-12)
    assert V("H_msub4_mismatch_percent") == pytest.approx(100 * h6["architectures"]["4"]["mismatch_rel"], rel=1e-12)
    assert V("P_combined_L2_p95_max_error") == pytest.approx(json.load(open(R / "phase6" / "6b_robustness.json"))["combined_L2_p95_max_err"], rel=1e-12)


def test_headline_values_match_reports():
    """Values quoted in the abstract / report (rounded) against the machine-readable results."""
    assert V("A_fractional_bandwidth") * 100 == pytest.approx(14.8, abs=0.05)
    assert V("B_D0_delay_quantum_cycles") == pytest.approx(0.981, abs=5e-4)
    assert V("H_msub4_mismatch_percent") == pytest.approx(-1.89, abs=0.005)
    assert V("I_msub4_trim_max_pointing_error") == pytest.approx(0.189, abs=5e-4)
    assert V("I_msub4_fixed_max_pointing_error") == pytest.approx(2.69, abs=0.005)
    assert V("J_msub4_trim_R_squint") == pytest.approx(0.070, abs=5e-4)
    assert V("K_msub4_trim_strongest_other_lobe_db") == pytest.approx(-12.4, abs=0.05)
    assert V("L_msub4_trim_strongest_replica_lobe_db") == pytest.approx(-16.6, abs=0.05)
    assert V("P_combined_L2_p95_max_error") == pytest.approx(0.35, abs=0.005)
    assert V("C_increment_max_abs_error_percent") < 0.04
    assert V("D_amp_imbalance_max_db") < 0.023
    assert V("E_trimmed_phase_residual_max_rad") == pytest.approx(0.026, abs=5e-4)
    assert abs(V("G_phase5_strongest_other_lobe_db")) < 1e-4          # 0 dB grating lobes


def test_phase4_recomputation_matches_handoff_trim_error_and_phase5_link():
    p5 = json.load(open(R / "phase5" / "phase5_handoff.json"))
    assert V("E_trimmed_phase_residual_max_rad") == pytest.approx(p5["30"]["phase_link_trim"]["max_abs_eps_rad"], rel=1e-9)


def test_every_value_has_source():
    for k, v in FR.items():
        assert v["source"] and "value" in v, k


def test_handoff_references_exist():
    for f in list(HO["figures"]) + HO["figures_pdf"] + HO["documents"] + list(HO["upstream_handoffs"].values()) + [HO["final_results_file"]]:
        assert (ROOT / f).exists(), f
    for f in ("docs/final_tables.md",):
        assert (ROOT / f).exists()
    assert HO["authoritative_commits"]["phase6_final_head"] == "3c10780"
    assert len(HO["figures"]) >= 6
    assert HO["wording_decisions"]["dispersion_engineered"].startswith("dropped")


def test_phase4_5_6_results_unchanged():
    for rel, h in HO["protected_result_sha256"].items():
        assert hashlib.sha256(open(R / rel, "rb").read()).hexdigest() == h, rel


def test_protected_files_equal_phase6_commit():
    """The recorded hashes must also equal the Phase 6 final commit (skipped silently if git history is unavailable)."""
    for rel, h in list(HO["protected_result_sha256"].items())[:6]:
        try:
            blob = subprocess.run(["git", "show", f"3c10780:results/{rel}"], cwd=ROOT, capture_output=True, check=True).stdout
        except Exception:
            pytest.skip("git history unavailable")
        assert hashlib.sha256(blob).hexdigest() == h, rel


def test_final_tables_numbers_come_from_results():
    t = (ROOT / "docs" / "final_tables.md").read_text()
    for s in (f"{V('I_msub4_trim_max_pointing_error'):.3f}", f"{V('J_msub4_trim_R_squint'):.3f}", f"{V('L_msub4_trim_strongest_replica_lobe_db'):.1f}", f"{V('H_msub4_mismatch_percent'):+.2f}",
              f"{V('P_combined_L2_p95_max_error'):.3f}", f"{V('M_msub4_trim_min_peak_power'):.4f}"):
        assert s in t, s


def test_manuscript_key_numbers_present():
    t = (ROOT / "docs" / "manuscript_draft.md").read_text()
    for s in ("14.8 %", "0.981", "1.9 %", "0.19°", "0.070", "−12.4 dB", "−16.6 dB", "0.35°", "2.69°"):
        assert s in t, s
    # revised manuscript: bare [REF] placeholders replaced by verified citations
    assert "[REF]" not in t and "not dispersion-engineered" in t
    low = t.lower()
    assert "first demonstration" not in low and "novel method" not in low
    assert "0.2887" in t and "29.377" in t and "26–34°" in t


def test_readme_script_references_exist():
    t = (ROOT / "README.md").read_text()
    missing = []
    for m in re.findall(r"`experiments/(phase\w+)/(\w+\.py)`", t):
        if not (ROOT / "experiments" / m[0] / m[1]).exists():
            missing.append(m)
    for m in re.findall(r"\$PY (\d\w+\.py)", t):
        if not list((ROOT / "experiments").glob(f"phase*/{m}")):
            missing.append(m)
    assert not missing, missing


def test_doc_references_exist():
    miss = []
    for f in list((ROOT / "docs").glob("*.md")) + [ROOT / "README.md"]:
        for m in set(re.findall(r"\]\(([\w./\-]+\.(?:md|py|png|json))\)", f.read_text())):
            if not (f.parent / m).exists() and not (ROOT / m).exists():
                miss.append((f.name, m))
    assert not miss, miss
