"""Phase 6E: results/phase6/phase6_handoff.{npz,json} from 6a_subarrays and 6b_robustness (selected arrays, definitions, restrictions, provenance)."""
from _common import RES, np, CASES, MSUBS, THETA0
import json
from pcttd.io import save_json

z = np.load(RES / "6a_subarrays.npz"); j = json.load(open(RES / "6a_subarrays.json")); zb = np.load(RES / "6b_robustness.npz"); jb = json.load(open(RES / "6b_robustness.json"))
best = j["best"]; BM = best["msub"]
f = z["f"]
keep = dict(f=f, mask_common=np.ones(len(f), bool), f0=z["f0"], theta0=THETA0, d_elem=z["d_e"], lambda0=z["lam0"], tau_phc=z["tau_phc"], theta_pat=z["theta_pat"], pattern_freqs=f[z["i_lpf"]],
            best_msub=BM, best_case=best["case"], best_id=best["id"], best_n_elements=4 * BM, best_subarray_pitch=BM * float(z["d_e"]),
            best_tau_required=j["arch"][str(BM)]["tau_req"], best_delay_mismatch_rel=j["arch"][str(BM)]["mismatch_rel"], best_delay_mismatch_cycles=j["arch"][str(BM)]["mismatch_cycles"],
            msub_list=np.array(MSUBS), tau_required=np.array([j["arch"][str(M)]["tau_req"] for M in MSUBS]), rho=np.array([j["arch"][str(M)]["rho"] for M in MSUBS]),
            delay_mismatch_rel=np.array([j["arch"][str(M)]["mismatch_rel"] for M in MSUBS]), delay_mismatch_cycles=np.array([j["arch"][str(M)]["mismatch_cycles"] for M in MSUBS]),
            n_elements=np.array([j["arch"][str(M)]["n_el"] for M in MSUBS]), subarray_pitch=np.array([j["arch"][str(M)]["subarray_pitch"] for M in MSUBS]),
            rho_grid=z["rho_grid"], mismatch_max_err_retrim=z["mm4_max_err_retrim"], mismatch_err_f0_noretrim=z["mm4_err_f0_noretrim"], pattern_cases=np.array(CASES), mc_nominal_theta_peak=zb["theta_nominal"])
for res in (30, 40):
    for M in MSUBS:
        for k in CASES + ["raw_unit", "trim_unit"]:
            for q in ("theta_peak", "gain", "hpbw", "strongest_db", "strongest_theta", "grating_db", "n_3db"):
                keep[f"res{res}_M{M}_{k}_{q}"] = z[f"res{res}_M{M}_{k}_{q}"]
            keep[f"res{res}_M{M}_{k}_pointing_error"] = z[f"res{res}_M{M}_{k}_theta_peak"] - THETA0
for M in MSUBS:
    for k in CASES:
        keep[f"patterns_res30_M{M}_{k}"] = z[f"pat_M{M}_{k}"]
np.savez(RES / "phase6_handoff.npz", **keep)
meta = dict(source="results/phase6/6a_subarrays.{npz,json}, 6b_robustness.{npz,json}; Phase 4/5 inputs unchanged", frequency="mask_common 0.231-0.268 (domain-limited); all arrays are over this band",
            best=j["best"], candidates=j["candidates"], architectures=j["arch"], case_B=j["arch_B"], sign_conventions=j["definitions"], restrictions=None,
            metrics_res30={k: j["30"][k] for k in j["30"]}, metrics_res40={k: j["40"][k] for k in j["40"]}, res_check=j["res_check"], quantization_analytic=j["quantization"], mismatch_sweep_msub4_note="ideal-delay hypothetical, arrays mismatch_* in the npz",
            robustness=dict(seed=jb["seed"], n_samples=jb["n_samples"], deterministic=jb["deterministic"], monte_carlo=jb["monte_carlo"], budget=jb["budget"], combined_L2_p95_max_err=jb["combined_L2_p95_max_err"],
                            combined_L2_ok=jb["combined_L2_ok"], engine_validation=jb["engine_validation"], quantization_scan=jb["quantization_scan"], definitions=jb["definitions"],
                            p4_trim_intercept_pred_deg=jb["p4_trim_intercept_pred_deg"], retrim_slope_pred_low_high_deg=jb["retrim_slope_pred_low_high_deg"]),
            subarray_replica_geometry=dict(count_low_f0_high=j["subarray_level_replica_count_low_f0_high"], angles_f0=j["subarray_level_replica_angles_f0"], half_wave_visible_replicas_max=j["half_wave_visible_replicas_max"]),
            restrictions_list=["quantitative use only inside mask_common 0.231-0.268 (post hoc Phase 2 validation edge; domain-limited)",
                               "2-D, lossless, uniform, no disorder/coupling/thermal/bends/length-dependent loss; no fabrication model (tolerances are system-level estimates)",
                               "four channels N = 8-11 of one flat-n_g D0 mode (one-period increments only); Msub = 8 with 2-period steps and fractional periods not simulated",
                               "raw PhC has constant phase-intercept errors up to ~2.7 rad; trimmed = constant phase per channel at f0 (never ideal TTD); delay slope untouched",
                               "delay quantum 0.9811 carrier cycle per period; Msub = 4 matches to -1.89 %, other Msub do not",
                               "design lobe tracked in a window of half the subarray-level grating spacing; subarray-level replica (quantization) lobes are not removed, only suppressed by the element factor",
                               "intra-subarray phase is fixed (set at f0): the subarray element factor squints and pulls the lobe by up to ~0.16 deg even with ideal inter-subarray TTD",
                               "physical mapping is by normalized scaling only; no claim that the 2-D dielectric PhC is a microwave feed network; RF-photonic use is a different (envelope-delay) regime",
                               "common minimum run time 3500 and local-quadratic/band-fit estimators inherited from Phase 4"],
            provenance_note="see provenance field; Phase 4/5 inputs verified by sha256 in tests/test_phase6.py")
meta.pop("restrictions")
save_json(RES / "phase6_handoff.json", meta)
chk = np.load(RES / "phase6_handoff.npz"); print(len(chk.files), "arrays; best", str(chk["best_id"]))
