"""Phase 7A: one authoritative table of final numbers, every value read from (or computed from) machine-readable result files.
Output: results/final/final_results.json  (key -> {value, unit, source, note})."""
from _common import RES, jload, nload, json, np

out = {}


def put(key, value, unit, source, note=""):
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    out[key] = dict(value=value, unit=unit, source=source, note=note)


p4 = nload("phase4/phase4_handoff.npz"); p4j = jload("phase4/phase4_handoff.json")
p5j = jload("phase5/phase5_handoff.json"); p6j = jload("phase6/phase6_handoff.json"); p6r = jload("phase6/6b_robustness.json")

# ---- A. validated band ------------------------------------------------------------------------------------------------------------------
lo, hi = p4j["common_band"]; fc = 0.5 * (lo + hi)
put("A_band_low", lo, "c/a", "phase4_handoff.json:common_band", "post hoc Phase 2 validation limits, not physical band edges")
put("A_band_high", hi, "c/a", "phase4_handoff.json:common_band")
put("A_band_center_f0", p4j["f0"], "c/a", "phase4_handoff.json:f0")
put("A_fractional_bandwidth", (hi - lo) / fc, "fraction", "computed (f_high - f_low)/f_center", "domain-limited lower bound; unrestricted ~18 % in Phase 3 (not validated)")

# ---- B. D0 delay per period (Phase 3) ----------------------------------------------------------------------------------------------------
for d in ("D0", "D1", "D2"):
    for r in (30, 40):
        b = jload(f"phase3/3e_{d}_res{r}.json")["bands"]["tau_fit|validated|tol10"]
        for k, u in (("median_tau", "a/c"), ("frac_bw", "fraction"), ("f_low", "c/a"), ("f_high", "c/a"), ("gbp", "1"), ("ripple", "fraction"), ("L_periods", "periods"),
                     ("L_wavelengths", "lambda"), ("gdd_band_avg_linfit", "(a/c)^2 per period"), ("R_max", "1"), ("T_min", "1")):
            put(f"B_{d}_res{r}_{k}", b[k], u, f"phase3/3e_{d}_res{r}.json:bands[tau_fit|validated|tol10].{k}")
tau = out["B_D0_res30_median_tau"]["value"]; f0 = p4j["f0"]
put("B_D0_delay_per_period_a_over_c", tau, "a/c", "phase3/3e_D0_res30.json (equals phase4_handoff.dtau_target)")
put("B_D0_delay_quantum_cycles", f0 * tau, "carrier cycles", "computed f0 * median_tau")
assert abs(tau - p4j["dtau_target"]) < 1e-4   # handoff stores the value rounded to 4 decimals

# ---- C-E. Phase 4 channel family (recomputed from the handoff) ---------------------------------------------------------------------------
f = p4["f"]; m = p4["mask_common"].astype(bool); H = p4["H_raw"]; Ht = p4["H_trim"]; Nn = p4["N"]
ph = p4["phase_abs_unwrapped"]
slope = np.array([np.polyfit(f[m], ph[i, m], 1)[0] / (2 * np.pi) for i in range(4)])   # delay = dphi/(2 pi df)
inc = np.diff(slope); tgt = float(p4["dtau_target"])
put("C_channel_delays_bandfit", slope.tolist(), "a/c", "computed: LS slope of phase_abs_unwrapped over mask_common / (2 pi)")
put("C_increments_bandfit", inc.tolist(), "a/c", "computed: adjacent differences")
put("C_increment_error_vs_target_percent", (100 * (inc - tgt) / tgt).tolist(), "%", "computed vs dtau_target")
put("C_increment_max_abs_error_percent", float(np.max(np.abs(100 * (inc - tgt) / tgt))), "%", "computed")
put("D_amp_imbalance_max_db", float(np.max(p4["amp_imbalance_db"][m])), "dB", "phase4_handoff.npz:amp_imbalance_db over mask_common")
put("D_T_min_common", float(np.min(p4["T"][:, m])), "1", "phase4_handoff.npz:T")
put("D_T_max_common", float(np.max(p4["T"][:, m])), "1", "phase4_handoff.npz:T", "excess over 1 is numerical")
rel = np.angle(Ht / Ht[0]); relu = np.unwrap(rel, axis=1)
i0 = int(np.argmin(np.abs(f - f0)))
err = relu - 2 * np.pi * f[None, :] * (np.arange(4)[:, None] * tgt)
err = err - 2 * np.pi * np.round(err[:, [i0]] / (2 * np.pi))          # branch: error at f0 within (-pi, pi]
put("E_trimmed_phase_residual_max_by_channel_rad", np.max(np.abs(err[:, m] - err[:, [i0]]), axis=1).tolist(), "rad", "computed from phase4_handoff.npz:H_trim")
put("E_trimmed_phase_residual_max_rad", float(np.max(np.abs(err[:, m] - err[:, [i0]]))), "rad", "computed")
put("E_trimmed_phase_residual_max_deg", float(np.degrees(np.max(np.abs(err[:, m] - err[:, [i0]])))), "deg", "computed")
put("E_raw_intercept_error_at_f0_rad", p4["trim_delta"].tolist(), "rad", "phase4_handoff.npz:trim_delta", "wrapped; up to 2.7 rad")
put("E_raw_intercept_error_max_abs_rad", float(np.max(np.abs(p4["trim_delta"]))), "rad", "phase4_handoff.npz:trim_delta")

# ---- F-G. Phase 5 ----------------------------------------------------------------------------------------------------------------------
s5 = p5j["setup"]; c5 = p5j["30"]
put("F_phase5_pitch_over_lambda0", s5["d_over_lambda0"], "lambda0", "phase5_handoff.json:setup.d_over_lambda0")
for k in ("fixed", "ideal", "raw", "trim"):
    put(f"F_phase5_max_pointing_error_{k}", c5["max_abs_err"][k], "deg", f"phase5_handoff.json:30.max_abs_err.{k}")
    put(f"F_phase5_R_squint_{k}", c5["R_squint"][k], "1", f"phase5_handoff.json:30.R_squint.{k}")
put("F_phase5_trim_vs_ideal_max", c5["vs_ideal_max"]["trim"], "deg", "phase5_handoff.json:30.vs_ideal_max.trim")
put("G_phase5_strongest_other_lobe_db", max(c5["strongest_db_max"][k] for k in ("fixed", "ideal", "raw", "trim")), "dB re design lobe", "phase5_handoff.json:30.strongest_db_max",
    "0 dB: equal-strength spatial-alias grating lobes (d = 1.96 lambda0)")
put("G_phase5_fixed_squint_low_f0_high", c5["err_low_f0_high"]["fixed"], "deg", "phase5_handoff.json:30.err_low_f0_high.fixed")

# ---- H-P. Phase 6 ----------------------------------------------------------------------------------------------------------------------
a4 = p6j["architectures"]["4"]; m4 = p6j["metrics_res30"]["4"]
for M in ("1", "2", "4", "8"):
    A = p6j["architectures"][M]
    put(f"H_msub{M}_required_cycles", A["req_cycles_f0"], "carrier cycles", f"phase6_handoff.json:architectures.{M}.req_cycles_f0")
    put(f"H_msub{M}_rho", A["rho"], "1", f"phase6_handoff.json:architectures.{M}.rho")
    put(f"H_msub{M}_mismatch_percent", 100 * A["mismatch_rel"], "%", f"phase6_handoff.json:architectures.{M}.mismatch_rel")
    mm = p6j["metrics_res30"][M]
    put(f"I_msub{M}_trim_max_pointing_error", mm["max_abs_err"]["trim"], "deg", f"phase6_handoff.json:metrics_res30.{M}.max_abs_err.trim")
    put(f"J_msub{M}_trim_R_squint", mm["R_squint"]["trim"], "1", f"phase6_handoff.json:metrics_res30.{M}.R_squint.trim")
    put(f"K_msub{M}_trim_strongest_other_lobe_db", mm["strongest_db_max"]["trim"], "dB", f"phase6_handoff.json:metrics_res30.{M}.strongest_db_max.trim")
put("H_phc_delay_quantum_cycles", a4["phc_cycles_f0"], "carrier cycles", "phase6_handoff.json:architectures.4.phc_cycles_f0")
put("H_msub4_required_delay_cycles", a4["req_cycles_f0"], "carrier cycles", "phase6_handoff.json:architectures.4.req_cycles_f0")
put("H_msub4_mismatch_percent", 100 * a4["mismatch_rel"], "%", "phase6_handoff.json:architectures.4.mismatch_rel")
put("I_msub4_fixed_max_pointing_error", m4["max_abs_err"]["fixed"], "deg", "phase6_handoff.json:metrics_res30.4.max_abs_err.fixed")
put("I_msub4_ideal_sub_max_pointing_error", m4["max_abs_err"]["ideal_sub"], "deg", "phase6_handoff.json:metrics_res30.4.max_abs_err.ideal_sub", "fixed intra-subarray phases squint even for ideal TTD")
put("I_msub4_trim_rms_pointing_error", m4["rms_err"]["trim"], "deg", "phase6_handoff.json:metrics_res30.4.rms_err.trim")
put("I_msub4_trim_vs_ideal_sub_max", m4["vs_ideal_sub_max"]["trim"], "deg", "phase6_handoff.json:metrics_res30.4.vs_ideal_sub_max.trim")
put("I_msub4_trim_err_low_f0_high", m4["err_low_f0_high"]["trim"], "deg", "phase6_handoff.json:metrics_res30.4.err_low_f0_high.trim")
put("I_msub4_tracking_tolerance", m4["tol_pt"], "deg", "phase6_handoff.json:metrics_res30.4.tol_pt")
put("I_msub4_raw_max_pointing_error", m4["max_abs_err"]["raw"], "deg", "phase6_handoff.json:metrics_res30.4.max_abs_err.raw")
put("I_msub4_trim_p4_max_pointing_error", p6j["metrics_res30"]["4"]["max_abs_err"].get("trim_p4"), "deg", "phase6_handoff.json:metrics_res30.4.max_abs_err.trim_p4", "Phase 4 trim used unchanged")
put("J_msub4_trim_R_squint", m4["R_squint"]["trim"], "1", "phase6_handoff.json:metrics_res30.4.R_squint.trim")
put("K_msub4_trim_strongest_other_lobe_db", m4["strongest_db_max"]["trim"], "dB", "phase6_handoff.json:metrics_res30.4.strongest_db_max.trim", "ordinary sidelobe")
put("L_msub4_trim_strongest_replica_lobe_db", m4["grating_db_max"]["trim"], "dB", "phase6_handoff.json:metrics_res30.4.grating_db_max.trim", "subarray replica (quantization) lobe at analytic positions")
put("M_msub4_trim_min_peak_power", m4["gain_min"]["trim"], "1 (re ideal 1)", "phase6_handoff.json:metrics_res30.4.gain_min.trim")
put("M_msub4_trim_min_peak_power_db", 10 * np.log10(m4["gain_min"]["trim"]), "dB", "computed")
put("N_msub4_trim_hpbw_low_f0_high", m4["hpbw_low_f0_high"]["trim"], "deg", "phase6_handoff.json:metrics_res30.4.hpbw_low_f0_high.trim")
rc = p6j["res_check"]["4"]["trim"]
put("O_res30_vs_res40_max_dtheta", rc["max_dtheta_deg"], "deg", "phase6_handoff.json:res_check.4.trim.max_dtheta_deg")
put("O_res30_vs_res40_limit", rc["tol_third"], "deg", "phase6_handoff.json:res_check.4.trim.tol_third")
put("O_res30_vs_res40_dR_squint", rc["d_R_squint"], "1", "phase6_handoff.json:res_check.4.trim.d_R_squint")
put("O_res30_vs_res40_dstrongest_db", rc["d_strongest_db"], "dB", "phase6_handoff.json:res_check.4.trim.d_strongest_db")
put("P_combined_L2_p95_max_error", p6r["combined_L2_p95_max_err"], "deg", "phase6/6b_robustness.json:combined_L2_p95_max_err")
put("P_nominal_max_error", p6r["nominal_max_err_deg"], "deg", "phase6/6b_robustness.json:nominal_max_err_deg")
bud = p6r["budget"]
put("P_budget_phase_p95_limit_0p5deg", bud["phase"]["p95"]["0.5"], "deg", "phase6/6b_robustness.json:budget.phase.p95")
put("P_budget_phase_p95_limit_1p0deg", bud["phase"]["p95"]["1.0"], "deg", "phase6/6b_robustness.json:budget.phase.p95")
put("P_budget_delay_uncal_p95_limit_0p5deg_percent", bud["delay_uncal"]["p95"]["0.5"], "% of q*tau_PhC", "phase6/6b_robustness.json:budget.delay_uncal.p95")
put("P_budget_delay_uncal_p95_limit_1p0deg_percent", bud["delay_uncal"]["p95"]["1.0"], "% of q*tau_PhC", "phase6/6b_robustness.json:budget.delay_uncal.p95")
put("P_mc_seed", p6r["seed"], "1", "phase6/6b_robustness.json:seed"); put("P_mc_samples", p6r["n_samples"], "1", "phase6/6b_robustness.json:n_samples")
put("H_subarray_pitch_over_lambda0", a4["subarray_pitch_over_lambda0"], "lambda0", "phase6_handoff.json:architectures.4")
put("H_n_elements", a4["n_el"], "1", "phase6_handoff.json:architectures.4")
put("H_element_pitch_over_lambda0", a4["d_e_over_lambda0"], "lambda0", "phase6_handoff.json:architectures.4")

json.dump(out, open(RES / "final_results.json", "w"), indent=1)
print(len(out), "values written to results/final/final_results.json")
for k in ("A_fractional_bandwidth", "B_D0_delay_quantum_cycles", "C_increments_bandfit", "C_increment_max_abs_error_percent", "D_amp_imbalance_max_db", "E_trimmed_phase_residual_max_rad", "E_raw_intercept_error_max_abs_rad",
          "F_phase5_max_pointing_error_trim", "G_phase5_strongest_other_lobe_db", "H_msub4_mismatch_percent", "I_msub4_trim_max_pointing_error", "J_msub4_trim_R_squint", "K_msub4_trim_strongest_other_lobe_db",
          "L_msub4_trim_strongest_replica_lobe_db", "M_msub4_trim_min_peak_power", "N_msub4_trim_hpbw_low_f0_high", "O_res30_vs_res40_max_dtheta", "P_combined_L2_p95_max_error", "I_msub4_trim_p4_max_pointing_error"):
    print(k, out[k]["value"])
