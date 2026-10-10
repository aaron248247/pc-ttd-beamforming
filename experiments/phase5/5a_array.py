"""Phase 5A: four-element array with the actual Phase 4 transfer functions (res 30 production, res 40 confirmation).
Writes results/phase5/5a_array.{npz,json}."""
from _common import P4, RES, np, THETA0, CASES
from pcttd import array as A, beam as B
from pcttd.io import save_json, provenance

z = np.load(P4 / "phase4_handoff.npz", allow_pickle=True)
fa = z["f"]; mk = z["mask_common"]; f = fa[mk]
f0 = float(z["f0"]); dtau = float(z["dtau_target"]); N = z["N"]
d = dtau / np.sin(np.deg2rad(THETA0)); lam0 = 1 / f0
fmax = f.max(); win = 1 / (2 * fmax * d)
tau = np.arange(4) * dtau
i0 = int(np.argmin(abs(f - f0))); I = [0, i0, len(f) - 1]
assert abs(f[i0] - f0) < 1e-9


def case_weights(res):
    Hr = z["H_raw" + ("_res40" if res == 40 else "")][:, mk]; Ht = z["H_trim" + ("_res40" if res == 40 else "")][:, mk]
    return dict(fixed=A.fixed_phase_weights(f, tau, f0), ideal=A.ideal_ttd_weights(f, tau), raw=Hr, trim=Ht,
                raw_unit=Hr / np.abs(Hr), trim_unit=Ht / np.abs(Ht))


def eps_phase(H):
    """Residual phase error vs ideal TTD, continuous in f and on the linear-in-m branch at f0 (rad)."""
    rel = np.unwrap(np.angle(H / H[0]), axis=1)
    e = rel - 2 * np.pi * f[None, :] * (np.arange(4) * dtau)[:, None]
    e = e + 2 * np.pi * np.round((np.angle(np.exp(1j * e[:, i0])) - e[:, i0]) / (2 * np.pi))[:, None]   # wrap f0 value to (-pi, pi]
    ref = np.unwrap(e[:, i0]); e = e + (ref - e[:, i0])[:, None]                                      # unwrap along m so the m-dependence is linear
    return e


out = {}; js = {}
for res in (30, 40):
    W = case_weights(res)
    R = {k: B.analyze_case(w, f, d, THETA0, win) for k, w in W.items()}
    out[res] = R
    err = {k: R[k]["theta_peak"] - THETA0 for k in R}
    fixed_max = np.abs(err["fixed"]).max()
    js[res] = dict(
        max_abs_err={k: float(np.abs(err[k]).max()) for k in R}, rms_err={k: float(np.sqrt(np.mean(err[k] ** 2))) for k in R},
        err_low_f0_high={k: [float(err[k][i]) for i in I] for k in R},
        R_squint={k: float(np.abs(err[k]).max() / fixed_max) for k in R},
        R_squint_spread={k: float(np.abs(R[k]["theta_peak"] - R[k]["theta_peak"][i0]).max() / fixed_max) for k in R},
        vs_ideal_max={k: float(np.abs(R[k]["theta_peak"] - R["ideal"]["theta_peak"]).max()) for k in R},
        vs_ideal_rms={k: float(np.sqrt(np.mean((R[k]["theta_peak"] - R["ideal"]["theta_peak"]) ** 2))) for k in R},
        gain_low_f0_high={k: [float(R[k]["gain"][i]) for i in I] for k in R}, gain_min={k: float(R[k]["gain"].min()) for k in R},
        hpbw_low_f0_high={k: [float(R[k]["hpbw"][i]) for i in I] for k in R},
        strongest_db_low_f0_high={k: [float(R[k]["strongest_db"][i]) for i in I] for k in R},
        strongest_theta_low_f0_high={k: [float(R[k]["strongest_theta"][i]) for i in I] for k in R},
        strongest_db_max={k: float(R[k]["strongest_db"].max()) for k in R}, n3db_low_f0_high={k: [int(R[k]["n_3db"][i]) for i in I] for k in R},
        edge_db_low_f0_high={k: [[float(x) for x in R[k]["edge_db"][i]] for i in I] for k in R})
    W_ = W
    for key in ("raw", "trim"):
        eps = eps_phase(W_[key])
        pred, g = B.phase_slope_pointing(eps - eps[:, [0]] * 0 if False else eps, f, d, THETA0)
        js[res][f"phase_link_{key}"] = dict(max_abs_eps_rad=float(np.abs(eps - eps[0]).max()), g_slope_rad_per_element_f0=float(g[i0]),
                                            pred_shift_max_deg=float(np.abs(pred).max()),
                                            numeric_shift_vs_ideal_max_deg=float(np.abs(R[key]["theta_peak"] - R["ideal"]["theta_peak"]).max()),
                                            pred_minus_numeric_max_deg=float(np.abs(pred - (R[key]["theta_peak"] - R["ideal"]["theta_peak"])).max()))
        out[res][key + "_eps"] = eps; out[res][key + "_pred"] = pred

# window membership: is the raw lobe near the window edge?
W30 = case_weights(30)
s_peak = {k: np.sin(np.deg2rad(out[30][k]["theta_peak"])) for k in CASES}
js["window_use_max_frac"] = {k: float(np.abs(s_peak[k] - np.sin(np.deg2rad(THETA0))).max() / win) for k in CASES}

# patterns (theta grid 0.01 deg) at low / f0 / high, res 30
th = np.linspace(-90, 90, 18001)
pat = np.array([[A.array_factor(W30[k][:, [i]], [f[i]], d, np.sin(np.deg2rad(th)))[0] for i in I] for k in CASES])
pat = np.abs(pat) ** 2 / 16

# analytic fixed-phase check + grating lobes
analytic = np.rad2deg(np.arcsin(f0 / f * np.sin(np.deg2rad(THETA0))))
js["fixed_analytic_check_max_dev_deg"] = float(np.abs(out[30]["fixed"]["theta_peak"] - analytic).max())
grat = B.grating_lobe_angles(out[30]["ideal"]["theta_peak"], f, d)
js["grating_ideal_visible_count"] = [int(len(g)) for g in grat]
js["grating_ideal_angles_low_f0_high"] = [[float(x) for x in grat[i]] for i in I]

# unit-amplitude vs actual amplitude
R30 = out[30]
js["amp_diag"] = {k: dict(max_dtheta_deg=float(np.abs(R30[k]["theta_peak"] - R30[k + "_unit"]["theta_peak"]).max()),
                          max_dgain=float(np.abs(R30[k]["gain"] - R30[k + "_unit"]["gain"]).max()),
                          max_dgain_db=float(np.abs(10 * np.log10(R30[k]["gain"] / R30[k + "_unit"]["gain"])).max()),
                          max_dstrongest_db=float(np.abs(R30[k]["strongest_db"] - R30[k + "_unit"]["strongest_db"]).max())) for k in ("raw", "trim")}
js["amp_imbalance_db_max"] = float(z["amp_imbalance_db"][mk].max())

# res 30 vs 40
a, b = out[30], out[40]
js["res_check"] = {}
for k in CASES:
    js["res_check"][k] = dict(max_dtheta_deg=float(np.abs(a[k]["theta_peak"] - b[k]["theta_peak"]).max()),
                              d_max_abs_err=float(abs(js[30]["max_abs_err"][k] - js[40]["max_abs_err"][k])),
                              d_R_squint=float(abs(js[30]["R_squint"][k] - js[40]["R_squint"][k])),
                              d_gain=float(np.abs(a[k]["gain"] - b[k]["gain"]).max()),
                              d_vs_ideal_max=float(abs(js[30]["vs_ideal_max"][k] - js[40]["vs_ideal_max"][k])),
                              d_strongest_db=float(np.abs(a[k]["strongest_db"] - b[k]["strongest_db"]).max()))

hp_ideal0 = float(R30["ideal"]["hpbw"][i0])
js["setup"] = dict(theta0_deg=THETA0, f0=f0, dtau_target=dtau, d=d, d_over_lambda0=d / lam0, lambda0=lam0, window_halfwidth_sin=win,
                   f_low=float(f[0]), f_high=float(fmax), n_freq=int(len(f)), N=[int(x) for x in N], hpbw_ideal_f0_deg=hp_ideal0,
                   tol_pt_deg=0.05 * hp_ideal0, d_over_lambda_band=[d * f[0], d * fmax],
                   delay_quantum_cycles=f0 * dtau, window_theta_deg=[float(np.rad2deg(np.arcsin(0.5 - win))), float(np.rad2deg(np.arcsin(0.5 + win)))])
js["definitions"] = dict(power="|AF|^2/16 (normalized power; 1 = ideal equal-amplitude 4-element array at its peak); gain_db = 10 log10",
                         sign="MEEP exp(-i omega t): delay tau -> exp(+i omega tau); AF(theta,f) = sum_m H_m exp(-i 2 pi f m d sin theta)",
                         R_squint="max_f|theta_PhC - theta0| / max_f|theta_fixed - theta0| (includes absolute offset); R_squint_spread removes theta_PhC(f0)",
                         strongest="strongest local maximum of the full -90..90 power pattern other than the design lobe, dB relative to the design-lobe peak",
                         design_lobe="maximum inside |sin(theta)-sin(theta0)| <= 1/(2 f_max d), parabolic refinement, 40001-point grid",
                         unit_cases="|H_m| set to 1, phases kept; diagnostic only")
js["restrictions"] = ["quantitative use only inside mask_common 0.231-0.268 (post hoc Phase 2 validation edge; domain-limited)", "raw PhC has constant phase-intercept errors up to ~2.7 rad",
                      "trimmed = phase-intercept-trimmed PhC, one constant per channel at f0, amplitude and delay unchanged; never ideal TTD",
                      "delay quantum ~0.98 carrier cycle per lattice period", "2-D, lossless, uniform, no disorder/coupling/thermal/bends/length-dependent loss",
                      "four-channel agreement reflects one flat-n_g D0 mode at four lengths", "common minimum run time 3500 (post-data method correction)",
                      "pointwise central-difference increment noise-limited; local-quadratic and band-fit used"]
save_json(RES / "5a_array.json", js)
arrs = dict(f=f, mask_common_phase4=mk, f_all=fa, theta_grid=th, d=d, theta0=THETA0, f0=f0, dtau_target=dtau, window=win, pat_freqs=f[I])
for res in (30, 40):
    for k, r in out[res].items():
        if isinstance(r, dict):
            for q, v in r.items():
                arrs[f"res{res}_{k}_{q}"] = v
        else:
            arrs[f"res{res}_{k}"] = r
arrs["patterns_res30"] = pat; arrs["pattern_cases"] = np.array(CASES)
np.savez(RES / "5a_array.npz", **arrs)

import json
print(json.dumps({k: js[k] for k in ("setup", "fixed_analytic_check_max_dev_deg", "window_use_max_frac")}, indent=1))
for res in (30, 40):
    print("RES", res)
    for q in ("max_abs_err", "rms_err", "R_squint", "R_squint_spread", "vs_ideal_max", "gain_min", "hpbw_low_f0_high", "strongest_db_low_f0_high", "err_low_f0_high", "n3db_low_f0_high", "edge_db_low_f0_high"):
        print(q, json.dumps(js[res][q], default=float))
    for key in ("raw", "trim"):
        print("link", key, js[res][f"phase_link_{key}"])
print("amp", js["amp_diag"], js["amp_imbalance_db_max"]); print("res", json.dumps(js["res_check"], indent=0)); print(js["grating_ideal_visible_count"][0::20], js["grating_ideal_angles_low_f0_high"])
