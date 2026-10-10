"""Phase 6A: delay-granularity mitigation.  Hierarchical subarray arrays (Msub = 1, 2, 4, 8; four PhC channels = four subarrays) with the actual
Phase 4 transfer functions (res 30 production, res 40 check).  Writes results/phase6/6a_subarrays.{npz,json}."""
from _common import P4, P5, RES, np, THETA0, MSUBS, NQ, CASES, DIAG
import json
from pcttd import array as A, beam as B, subarray as S
from pcttd.io import save_json

z = np.load(P4 / "phase4_handoff.npz", allow_pickle=True)
fa = z["f"]; mk = z["mask_common"]; f = fa[mk]
f0 = float(z["f0"]); tau_phc = float(z["dtau_target"]); lam0 = 1 / f0; d_e = lam0 / 2
i0 = int(np.argmin(abs(f - f0))); I = [0, i0, len(f) - 1]; fmax = f.max()
assert abs(f[i0] - f0) < 1e-9
s0 = np.sin(np.deg2rad(THETA0))

# ---------------------------------------------------------------- architecture table
arch = {}
for M in MSUBS:
    tr = S.tau_required(M, d_e, THETA0)
    rho, mm, dcyc = S.mismatch(tau_phc, tr, f0)
    arch[M] = dict(msub=M, n_el=NQ * M, d_e=d_e, d_e_over_lambda0=0.5, subarray_pitch=M * d_e, subarray_pitch_over_lambda0=M / 2, tau_req=tr,
                   req_cycles_f0=f0 * tr, tau_phc=tau_phc, phc_cycles_f0=f0 * tau_phc, rho=rho, mismatch_rel=mm, mismatch_cycles=dcyc,
                   aperture_lambda0=NQ * M * d_e / lam0, aperture_span_lambda0=(NQ * M - 1) * d_e / lam0,
                   window_sin=S.window_halfwidth(fmax, M, d_e))
arch_B = dict(tau_half=d_e * s0, tau_half_cycles=f0 * d_e * s0, ratio_phc_over_half=tau_phc / (d_e * s0), d_e=d_e, n_el=4)

# ---------------------------------------------------------------- case builders
def channels(res, M):
    """Subarray channel transfer functions G[case][q, f] for one resolution and subarray size."""
    sfx = "_res40" if res == 40 else ""
    Hr = z["H_raw" + sfx][:, mk]; Hp4 = z["H_trim" + sfx][:, mk]
    tr = S.tau_required(M, d_e, THETA0)
    Gtrim, delta = S.retrim(Hr, f, tr, f0)
    G = dict(fixed=S.fixed_channels(f, NQ, tr, f0), ideal_elem=S.ideal_channels(f, NQ, tr), ideal_sub=S.ideal_channels(f, NQ, tr),
             raw=Hr, trim=Gtrim, trim_p4=Hp4, raw_unit=Hr / np.abs(Hr), trim_unit=Gtrim / np.abs(Gtrim))
    U = dict(fixed=S.intra_fixed(M, d_e, THETA0, f0, len(f)), ideal_elem=S.intra_ttd(M, d_e, THETA0, f))
    for k in CASES + DIAG:
        U.setdefault(k, U["fixed"])
    return G, U, delta


out = {}; js = {"arch": {str(M): arch[M] for M in MSUBS}, "arch_B": arch_B}
for res in (30, 40):
    out[res] = {}; js[str(res)] = {}
    for M in MSUBS:
        G, U, delta = channels(res, M)
        eng = S.Hier(f, NQ, M, d_e)
        win = arch[M]["window_sin"]
        R = {k: eng.analyze(G[k], U[k], THETA0, win) for k in CASES + DIAG}
        out[res][M] = R
        err = {k: R[k]["theta_peak"] - THETA0 for k in R}
        fixed_max = np.abs(err["fixed"]).max()
        hp_ideal0 = float(R["ideal_elem"]["hpbw"][i0]); tol = 0.05 * hp_ideal0
        r = dict(tol_pt=tol, hpbw_ideal_f0=hp_ideal0, retrim_delta_rad=delta,
                 max_abs_err={k: float(np.abs(err[k]).max()) for k in R}, rms_err={k: float(np.sqrt(np.mean(err[k] ** 2))) for k in R},
                 err_low_f0_high={k: [float(err[k][i]) for i in I] for k in R},
                 spread_ptp={k: float(np.ptp(R[k]["theta_peak"])) for k in R},
                 R_squint={k: float(np.abs(err[k]).max() / fixed_max) for k in R},
                 R_squint_spread={k: float(np.abs(R[k]["theta_peak"] - R[k]["theta_peak"][i0]).max() / fixed_max) for k in R},
                 vs_ideal_elem_max={k: float(np.abs(R[k]["theta_peak"] - R["ideal_elem"]["theta_peak"]).max()) for k in R},
                 vs_ideal_sub_max={k: float(np.abs(R[k]["theta_peak"] - R["ideal_sub"]["theta_peak"]).max()) for k in R},
                 gain_low_f0_high={k: [float(R[k]["gain"][i]) for i in I] for k in R}, gain_min={k: float(R[k]["gain"].min()) for k in R},
                 hpbw_low_f0_high={k: [float(R[k]["hpbw"][i]) for i in I] for k in R},
                 strongest_db_low_f0_high={k: [float(R[k]["strongest_db"][i]) for i in I] for k in R},
                 strongest_theta_low_f0_high={k: [float(R[k]["strongest_theta"][i]) for i in I] for k in R},
                 strongest_db_max={k: float(R[k]["strongest_db"].max()) for k in R},
                 n3db_low_f0_high={k: [int(R[k]["n_3db"][i]) for i in I] for k in R}, n3db_max={k: int(R[k]["n_3db"].max()) for k in R},
                 grating_db_max={k: float(R[k]["grating_db"].max()) for k in R}, grating_db_low_f0_high={k: [float(R[k]["grating_db"][i]) for i in I] for k in R},
                 n_replicas_low_f0_high={k: [int(R[k]["n_replicas"][i]) for i in I] for k in R},
                 ambiguity={k: bool(R[k]["strongest_db"].max() > -10.0) for k in R},
                 window_use_max_frac={k: float(np.abs(np.sin(np.deg2rad(R[k]["theta_peak"])) - s0).max() / win) for k in R})
        # criterion 6 flags
        r["tracking_ok"] = {k: bool(np.abs(err[k]).max() <= tol) for k in R}
        r["squint_suppressed"] = {k: bool(r["R_squint_spread"][k] < 0.5) for k in R}
        r["lobe_mitigated"] = {k: bool(R[k]["strongest_db"].max() <= -10.0) for k in R}
        # analytic subarray-level squint law for a retrimmed channel family: s - s0 = (f0/f - 1) s0 (1 - rho)
        rho = arch[M]["rho"]
        pred = np.rad2deg(np.arcsin(s0 + (f0 / f - 1) * s0 * (1 - rho))) if abs((f0 / f.min() - 1) * s0 * (1 - rho)) < 1 - s0 else np.full(len(f), np.nan)
        r["analytic_law_trim"] = dict(max_abs_dev_deg=float(np.nanmax(np.abs(pred - R["trim"]["theta_peak"]))), pred_err_low_f0_high=[float(pred[i] - THETA0) for i in I])
        out[res][M]["trim_analytic"] = dict(pred=pred)
        js[str(res)][str(M)] = r

        # phase-link: residual phase vs ideal subarray TTD (subarray-level), slope over q -> pointing shift (first order)
        for key in ("raw", "trim"):
            rel = np.unwrap(np.angle(G[key] / G[key][0]), axis=1)
            e = rel - 2 * np.pi * f[None, :] * (np.arange(NQ) * S.tau_required(M, d_e, THETA0))[:, None]
            e = e + 2 * np.pi * np.round((np.angle(np.exp(1j * e[:, i0])) - e[:, i0]) / (2 * np.pi))[:, None]
            ref = np.unwrap(e[:, i0]); e = e + (ref - e[:, i0])[:, None]
            pred2, g = B.phase_slope_pointing(e, f, M * d_e, THETA0)
            num = R[key]["theta_peak"] - R["ideal_sub"]["theta_peak"]
            r[f"phase_link_{key}"] = dict(max_abs_eps_rad=float(np.abs(e - e[0]).max()), pred_shift_max_deg=float(np.abs(pred2).max()),
                                          numeric_vs_ideal_sub_max_deg=float(np.abs(num).max()), pred_minus_numeric_max_deg=float(np.abs(pred2 - num).max()))
            out[res][M][key + "_eps"] = e; out[res][M][key + "_pred"] = pred2

# ---------------------------------------------------------------- best architecture (criterion 6)
cand = []
for M in MSUBS:
    r = js["30"][str(M)]
    for k in ("raw", "trim", "trim_p4"):
        ok = r["tracking_ok"][k] and r["squint_suppressed"][k] and r["lobe_mitigated"][k]
        cand.append(dict(id=f"M{M}_{k}", msub=M, case=k, qualifies=bool(ok), max_abs_err=r["max_abs_err"][k], strongest_db_max=r["strongest_db_max"][k],
                         tracking_ok=r["tracking_ok"][k], squint_suppressed=r["squint_suppressed"][k], lobe_mitigated=r["lobe_mitigated"][k]))
q = [c for c in cand if c["qualifies"]]
pool = q if q else cand
best = min(pool, key=lambda c: (c["max_abs_err"], c["strongest_db_max"]))
js["candidates"] = cand; js["best"] = dict(best, any_qualifies=bool(q))

# ---------------------------------------------------------------- Msub=4 specific: delay mismatch -> pointing (hypothetical delay scale rho)
rho_grid = np.linspace(0.5, 2.0, 61)
M = 4; tr = S.tau_required(M, d_e, THETA0); eng4 = S.Hier(f, NQ, M, d_e); U4 = S.intra_fixed(M, d_e, THETA0, f0, len(f)); E4 = eng4.element_factor(U4); win4 = arch[M]["window_sin"]
mm = {"rho": rho_grid, "max_err_retrim": [], "err_f0_retrim": [], "max_err_noretrim": [], "err_f0_noretrim": []}
for rho in rho_grid:
    Gi = S.ideal_channels(f, NQ, rho * tr)               # ideal delay of scale rho x required, intercept = 0 (exact at f0 only if rho = 1 mod 1/(f0 ...))
    Gr, _ = S.retrim(Gi, f, tr, f0)
    for key, GG in (("retrim", Gr), ("noretrim", Gi)):
        th = eng4.peak(eng4.power(GG, E4), THETA0, win4)[0] - THETA0
        mm[f"max_err_{key}"].append(float(np.abs(th).max())); mm[f"err_f0_{key}"].append(float(th[i0]))
for k in ("max_err_retrim", "err_f0_retrim", "max_err_noretrim", "err_f0_noretrim"):
    mm[k] = np.array(mm[k])
js["mismatch_sweep_msub4"] = {k: v for k, v in mm.items()}

# ---------------------------------------------------------------- period quantization table (analytic): rho(Msub, k)
js["quantization"] = {f"M{Mq}_k{k}": dict(msub=Mq, k=k, tau_step=k * tau_phc, tau_req=Mq * lam0 / 4, rho=k * tau_phc / (Mq * lam0 / 4),
                                         squint_fraction=abs(1 - k * tau_phc / (Mq * lam0 / 4)))
                      for Mq in range(1, 13) for k in (1, 2, 3)}

# ---------------------------------------------------------------- lambda0/2 spatial aliasing check, ideal element-level TTD
gl = B.grating_lobe_angles(np.full(len(f), THETA0), f, d_e)
js["half_wave_visible_replicas_max"] = int(max(len(g) for g in gl))
gl_sub = {str(M): B.grating_lobe_angles(np.full(len(f), THETA0), f, M * d_e) for M in MSUBS}
js["subarray_level_replica_count_low_f0_high"] = {str(M): [int(len(gl_sub[str(M)][i])) for i in I] for M in MSUBS}
js["subarray_level_replica_angles_f0"] = {str(M): [float(x) for x in gl_sub[str(M)][i0]] for M in MSUBS}

# ---------------------------------------------------------------- patterns (res 30), low/f0/high
eng_by_M = {M: S.Hier(f, NQ, M, d_e) for M in MSUBS}
pat = {}
for M in MSUBS:
    G, U, _ = channels(30, M)
    for k in CASES:
        pat[f"M{M}_{k}"] = np.array([eng_by_M[M].pattern(G[k][:, i], U[k][:, i], i)[1][::2] for i in I]).astype(np.float32)
theta_pat = np.rad2deg(np.arcsin(eng_by_M[1].s))[::2]

# ---------------------------------------------------------------- Phase 5 reference (read only)
p5 = np.load(P5 / "phase5_handoff.npz"); p5j = json.load(open(P5 / "phase5_handoff.json"))
assert np.allclose(p5["f"], f)

# ---------------------------------------------------------------- res 30 vs 40 for every M and case (best is selected above from res 30)
js["res_check"] = {}
for M in MSUBS:
    tol = js["30"][str(M)]["tol_pt"]
    js["res_check"][str(M)] = {k: dict(max_dtheta_deg=float(np.abs(out[30][M][k]["theta_peak"] - out[40][M][k]["theta_peak"]).max()),
                                       d_max_abs_err=float(abs(js["30"][str(M)]["max_abs_err"][k] - js["40"][str(M)]["max_abs_err"][k])),
                                       d_R_squint=float(abs(js["30"][str(M)]["R_squint"][k] - js["40"][str(M)]["R_squint"][k])),
                                       d_gain=float(np.abs(out[30][M][k]["gain"] - out[40][M][k]["gain"]).max()),
                                       d_strongest_db=float(np.abs(out[30][M][k]["strongest_db"] - out[40][M][k]["strongest_db"]).max()),
                                       d_vs_ideal_elem=float(abs(js["30"][str(M)]["vs_ideal_elem_max"][k] - js["40"][str(M)]["vs_ideal_elem_max"][k])),
                                       tol_third=tol / 3) for k in CASES}

js["setup"] = dict(theta0_deg=THETA0, f0=f0, lambda0=lam0, d_e=d_e, tau_phc=tau_phc, nq=NQ, N_channels=[int(x) for x in z["N"]], f_low=float(f[0]), f_high=float(fmax), n_freq=len(f),
                   grid_points=int(len(eng_by_M[1].s)), seed_none="deterministic part", delay_quantum_cycles=f0 * tau_phc)
js["definitions"] = dict(power="|AF|^2 / N_el^2 (1 = ideal equal-amplitude array at its peak)", sign="exp(-i omega t); delay tau -> exp(+i omega tau)",
                         hierarchy="W_(q,p) = G_q(f) U_p(f); element (q,p) at (q Msub + p) lambda0/2; exact factorized full far field, not point-source subarrays",
                         trim="trim = constant phase per channel at f0 chosen for the REQUIRED step q tau_req(Msub); trim_p4 = Phase 4 H_trim as is",
                         window="|sin - sin(theta0)| <= 1/(2 f_max Msub d_e); 40001-point sin grid, parabolic refinement (production grid finer than the 20001 lobe grid of Revision 6, same method)",
                         strongest="strongest local maximum other than the design lobe over -90..90 deg, dB re design peak; ambiguity = strongest > -10 dB")
save_json(RES / "6a_subarrays.json", js)
arrs = dict(f=f, f0=f0, i_lpf=np.array(I), theta_pat=theta_pat, rho_grid=rho_grid, tau_phc=tau_phc, d_e=d_e, lam0=lam0)
for k, v in pat.items():
    arrs[f"pat_{k}"] = v
for res in (30, 40):
    for M in MSUBS:
        for k in CASES + DIAG:
            for q_, v in out[res][M][k].items():
                arrs[f"res{res}_M{M}_{k}_{q_}"] = v
        for key in ("raw", "trim"):
            arrs[f"res{res}_M{M}_{key}_eps"] = out[res][M][key + "_eps"]; arrs[f"res{res}_M{M}_{key}_pred"] = out[res][M][key + "_pred"]
        arrs[f"res{res}_M{M}_trim_analytic"] = out[res][M]["trim_analytic"]["pred"]
for k, v in mm.items():
    arrs[f"mm4_{k}"] = v
np.savez(RES / "6a_subarrays.npz", **arrs)

print(json.dumps(arch_B, indent=0))
for M in MSUBS:
    a = arch[M]; print("M", M, {k: round(v, 5) if isinstance(v, float) else v for k, v in a.items()})
for M in MSUBS:
    r = js["30"][str(M)]
    print("== M", M, "tol", round(r["tol_pt"], 3), "hpbw_ideal", round(r["hpbw_ideal_f0"], 3))
    for k in CASES + DIAG:
        print(f"  {k:10s} maxerr {r['max_abs_err'][k]:8.4f} rms {r['rms_err'][k]:7.4f} Rsq {r['R_squint'][k]:7.3f} Rspr {r['R_squint_spread'][k]:7.3f} vsIdealSub {r['vs_ideal_sub_max'][k]:7.4f} "
              f"gmin {r['gain_min'][k]:6.3f} strong {r['strongest_db_max'][k]:7.2f} n3 {r['n3db_max'][k]} winfrac {r['window_use_max_frac'][k]:.2f} err[l,0,h] {np.round(r['err_low_f0_high'][k],3)}")
    print("  grating_db_max", {k: round(r["grating_db_max"][k], 2) for k in CASES}, "n_rep", r["n_replicas_low_f0_high"]["trim"])
    print("  analytic law dev", r["analytic_law_trim"], "link", r["phase_link_trim"])
print("best", js["best"]); print("half-wave replicas", js["half_wave_visible_replicas_max"], js["subarray_level_replica_count_low_f0_high"], js["subarray_level_replica_angles_f0"])
for M in MSUBS:
    print("res", M, {k: (round(v["max_dtheta_deg"], 4), round(v["d_R_squint"], 4), round(v["d_gain"], 4), round(v["d_strongest_db"], 3), round(v["tol_third"], 3)) for k, v in js["res_check"][str(M)].items() if k in ("trim", "raw", "trim_p4")})
