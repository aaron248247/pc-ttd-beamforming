"""Phase 7D: results/final/final_handoff.json (claim, metrics, commits, files, provenance). No large arrays; pointers to the Phase 4/5/6 handoffs."""
from _common import RES, RESULTS, ROOT, FIG, COMMITS, PROTECTED, jload, json
import hashlib
import subprocess

FR = jload("final/final_results.json")
V = lambda k: FR[k]["value"]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def git(*a):
    try:
        return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


FIGS = {"fig01_concept": "W1 channel and subarray architecture (2-D FDTD domain vs analytic beamforming model)", "fig02_verification": "FDTD verification (convergence, Bragg vs TMM, W1 n_g vs MPB)",
        "fig03_tradeoff": "D0/D1/D2 delay-bandwidth-dispersion trade", "fig04_channels": "four-channel transfer functions (amplitude, increments, raw and trimmed phase error)",
        "fig05_phase5_limit": "element-level array: squint and 0 dB grating lobes", "fig06_architecture": "delay match, pointing and sidelobe level versus Msub",
        "fig07_patterns": "Msub = 4 patterns at band low/center/high", "fig08_pointing": "pointing error versus frequency", "fig09_robustness": "system-level perturbation sensitivity"}
claim = ("In a 2-D lossless FDTD/MPB model, a plain W1 photonic-crystal waveguide gives delay increments of 3.93 a/c per period (0.981 carrier cycle, 0.04 % band-fit accuracy) over a validated 14.8 % fractional "
         "bandwidth; with one channel and one constant phase trim per four-element lambda0/2 subarray at 30 deg (required 1.000 cycle, 1.9 % mismatch) the maximum pointing error over the band is 0.19 deg against 2.69 deg "
         "for fixed phase, strongest sidelobe -12.4 dB, replica lobes <= -16.6 dB. This supports the feasibility of subarray-level PhC delay lines for beam-squint reduction; it does not demonstrate a device.")
h = dict(
    final_title="Subarray-Level Delay Lines from W1 Photonic-Crystal Waveguides for Beam-Squint Reduction: A Computational Study",
    final_central_claim=claim,
    wording_decisions=dict(dispersion_engineered="dropped (plain W1 selected; engineering tried and rejected in Phase 3)", broadband="qualified as 14.8 % validated, domain-limited fractional bandwidth",
                           true_time_delay="delay slope only; trimmed family = TTD slope plus one constant phase trim per channel, never ideal TTD", suitable="replaced by 'supports the feasibility of'",
                           coincidence="Msub = 4 match is n_g f0 about 0.98 at 30 deg, not a design principle"),
    authoritative_commits=COMMITS, phase7_starting_head=COMMITS["phase6_final_head"],
    git_head_at_generation=git("rev-parse", "--short", "HEAD"),
    validated_band=dict(low=V("A_band_low"), high=V("A_band_high"), center=V("A_band_center_f0"), fractional=V("A_fractional_bandwidth"), note="domain-limited; upper edge post hoc"),
    selected_architecture="Msub = 4 subarray-level delay: 4 PhC channels (N = 8-11), 16 lambda0/2 elements, one constant phase trim per channel re-chosen for the required step, theta0 = 30 deg",
    key_metrics={k: V(k) for k in ("A_fractional_bandwidth", "B_D0_delay_per_period_a_over_c", "B_D0_delay_quantum_cycles", "C_increments_bandfit", "C_increment_max_abs_error_percent", "D_amp_imbalance_max_db",
                                   "E_trimmed_phase_residual_max_rad", "E_raw_intercept_error_max_abs_rad", "F_phase5_max_pointing_error_trim", "F_phase5_max_pointing_error_fixed", "G_phase5_strongest_other_lobe_db",
                                   "H_msub4_mismatch_percent", "I_msub4_trim_max_pointing_error", "I_msub4_fixed_max_pointing_error", "J_msub4_trim_R_squint", "K_msub4_trim_strongest_other_lobe_db",
                                   "L_msub4_trim_strongest_replica_lobe_db", "M_msub4_trim_min_peak_power", "N_msub4_trim_hpbw_low_f0_high", "O_res30_vs_res40_max_dtheta", "P_combined_L2_p95_max_error")},
    final_results_file="results/final/final_results.json",
    restrictions=["2-D, lossless", "no fabrication disorder", "no channel coupling", "no bends", "no thermal effects", "no realistic couplers", "no propagation loss", "finite validated band only (0.231-0.268)",
                  "post hoc validated upper edge 0.268", "post-data Phase 4 convergence correction (common minimum run time)", "post-data Phase 4 estimator change", "constant phase trim required",
                  "residual subarray replica lobes (<= -16.6 dB)", "robustness is system-level, not fabrication-level", "only four PhC channels", "only one steering angle (30 deg) studied for the architectural match",
                  "Msub = 8 with two-period spacing not simulated", "no conventional-delay-line comparison beyond MPB-only reference", "no literature survey; no novelty claim", "no experiment"],
    historical_corrections=["Phase 0.5 geometry sign bug and branch-tracking bug (rerun)", "Phase 2 validated-band boundary chosen post hoc", "Phase 3 hole-shift sign bug fixed, D1/D2 rerun",
                            "Phase 4 first FDTD attempt discarded (under-convergence)", "Phase 4 estimator changed post hoc", "Phase 4 cascade step-count erratum",
                            "Phase 6 first Monte Carlo sign convention discarded and rerun", "Phase 6 post hoc diagnostic additions"],
    arxiv_readiness="B: nearly ready; manuscript gap = references and literature positioning; scientific gap = no conventional-line comparison beyond MPB-only (restricts claim to feasibility); no numerical validation gap blocking; no further simulation recommended",
    figures={f"figures/final/{k}.png": v for k, v in FIGS.items()}, figures_pdf=[f"figures/final/{k}.pdf" for k in FIGS],
    tables="docs/final_tables.md (Tables 1-7, 5b)",
    documents=["docs/final_report.md", "docs/manuscript_draft.md", "docs/reproducibility.md", "docs/final_tables.md", "docs/decisions.md", "docs/research_plan.md"],
    upstream_handoffs=dict(phase4="results/phase4/phase4_handoff.npz", phase5="results/phase5/phase5_handoff.npz", phase6="results/phase6/phase6_handoff.npz"),
    protected_result_sha256={p: sha(RESULTS / p) for p in PROTECTED},
    provenance=dict(generated_by="experiments/phase7/7d_handoff.py", inputs="results/final/final_results.json plus the Phase 4-6 result files listed in protected_result_sha256", simulations_run_in_phase7=False),
    reproduce="cd experiments/phase7 && python 7a_final_results.py && python 7b_figures.py && python 7c_tables.py && python 7d_handoff.py (see docs/reproducibility.md)")
for k in h["figures"]:
    assert (ROOT / k).exists(), k
for d in h["documents"]:
    assert (ROOT / d).exists(), d
json.dump(h, open(RES / "final_handoff.json", "w"), indent=1)
print("final_handoff.json written;", len(h["protected_result_sha256"]), "protected hashes")
