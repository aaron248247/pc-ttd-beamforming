"""Phase 4B: complex channel transfer functions H_m = S21_m of the four directly simulated D0 channels, relative phase/delay, errors, trims,
common usable band, cascade check.  usage: python 4b_analyze.py <res>  -> results/phase4/4b_res<res>.{npz,json}
Definitions: docs/phase2_5_criteria.md, Revision 4 (2026-10-10).  Conventions: a = c = 1, MEEP_SIGN = +1, tau_g = +d arg H / d omega.
"""
from _common import RES, RAW, ROOT, np, LENGTHS, DN, F0, DTAU_TARGET, VALID, N0
import json, sys
from pcttd import guided as GD, w1analysis as WA, ttdmetrics as TM, channels as C
from pcttd.io import save_json, provenance

P2 = ROOT / "results" / "phase2"
P3 = ROOT / "results" / "phase3"
D_MON = 3
CRIT = dict(T_min=0.5, defect_max=0.02, res_max=0.10, inc_tol=0.10, imbalance_db=1.0)


def lead_profiles(f):
    z = np.load(P2 / "2c_mode_profiles.npz")
    if not np.array_equal(z["f"], f):
        raise RuntimeError("frequency grid differs from the cached Phase 2 lead-mode profiles")
    return z["y"], z["Ey"], z["Hz"], float(z["H"])


def build(res):
    Ns = LENGTHS
    raws = {N: dict(np.load(RAW / f"D0_res{res}_N{N}.npz")) for N in Ns}
    f = raws[Ns[0]]["f"]
    mp3 = np.load(P3 / "3d_mpb_analysis.npz")
    assert np.allclose(f, mp3["f"])
    y, Ey, Hz, H = lead_profiles(f)
    S = {N: WA.modal_sparams(raws[N], y, Ey, Hz, H) for N in Ns}
    inband = mp3["D0_inband"].astype(bool)
    ref_pp = -2 * np.pi * np.where(inband, mp3["D0_kappa"], 0.0)          # approximate unwrapped phase advance per period (MPB)
    out = dict(f=f, ng_mpb=mp3["D0_ng"], ref_pp=ref_pp, inband=inband, N=np.array(Ns))
    H21 = np.array([S[N]["S21"] for N in Ns]); H11 = np.array([S[N]["S11"] for N in Ns])
    out.update(H=H21, S11=H11, T=np.abs(H21) ** 2, R=np.abs(H11) ** 2, res_in=np.array([S[N]["res_in"] for N in Ns]),
               res_out=np.array([S[N]["res_out"] for N in Ns]), about=np.array([S[N]["ab_out"] for N in Ns]),
               defect=np.array([WA.unitarity_defect(S[N]["S11"], S[N]["S21"]) for N in Ns]),
               t_end=np.array([raws[N]["t_end"] for N in Ns]), wall_s=np.array([raws[N]["wall_s"] for N in Ns]))
    # absolute unwrapped phase of H_m (branch from MPB over the full length N + 2 d_mon) and a branch-consistency check against plain unwrapping
    phi_abs = np.array([GD.unwrap_theta(np.angle(H21[i]), ref=ref_pp * (N + 2 * D_MON)) for i, N in enumerate(Ns)])
    out["phi_abs"] = phi_abs
    out["tau_abs"] = C.delays(f, phi_abs)
    out["phi_rel"] = C.relative_phase(H21, ref_per_period=ref_pp, dN=DN)
    for meth, key in (("central", ""), ("poly", "_poly")):
        out["Dtau" + key] = C.delays(f, out["phi_rel"], method=meth)
        out["inc" + key] = C.incremental_delays(f, H21, ref_per_period=ref_pp, dN=DN, method=meth)
    out["raw_cfg"] = np.array([raws[N]["cfg_json"] for N in Ns])
    out["raw_prov"] = np.array([raws[N]["prov_json"] for N in Ns])
    out["hit_max"] = np.array([raws[N]["hit_max_time"] for N in Ns])
    return out, S, raws


def masks(out, f, inc_key="inc", inc_tol=CRIT["inc_tol"]):
    dom = (f >= VALID[0] - 1e-9) & (f <= VALID[1] + 1e-9)
    m = dict(
        domain=dom, inband=out["inband"],
        T=np.all(out["T"] >= CRIT["T_min"], axis=0),
        unitarity=np.all(np.abs(out["defect"]) <= CRIT["defect_max"], axis=0),
        residual=np.all(np.maximum(out["res_in"], out["res_out"]) <= CRIT["res_max"], axis=0),
        increment=np.all(np.abs(out[inc_key] / DTAU_TARGET - 1.0) <= inc_tol, axis=0) & np.all(np.isfinite(out[inc_key]), axis=0),
        amplitude=C.amplitude_imbalance_db(out["H"]) <= CRIT["imbalance_db"])
    allm = C.common_usable_mask(list(m.values()))
    return m, allm


def band_stats(out, f, i0, i1, delta, eps_raw, eps_trim, inc_key="inc"):
    sl = slice(i0, i1 + 1)
    w = TM.omega(f)
    d = dict(i_low=int(i0), i_high=int(i1), f_low=float(f[i0]), f_high=float(f[i1]), f_center=float(0.5 * (f[i0] + f[i1])),
             frac_bw=TM.fractional_bandwidth(f[i0], f[i1]), n_samples=int(i1 - i0 + 1))
    d["T_min"] = float(out["T"][:, sl].min()); d["T_max"] = float(out["T"][:, sl].max()); d["R_max"] = float(out["R"][:, sl].max())
    d["res_max"] = float(max(out["res_in"][:, sl].max(), out["res_out"][:, sl].max()))
    d["unitarity_defect_absmax"] = float(np.abs(out["defect"][:, sl]).max())
    imb = C.amplitude_imbalance_db(out["H"])
    d["amp_imbalance_db_max"] = float(imb[sl].max()); d["amp_imbalance_db_f0"] = float(imb[int(np.argmin(abs(f - F0)))])
    d["amp_imbalance_note"] = "20 log10(max_m|H_m| / min_m|H_m|), amplitude ratio"
    d["amp_rel_to_ch0_range"] = [float((np.abs(out["H"][1:, sl]) / np.abs(out["H"][0, sl])).min()), float((np.abs(out["H"][1:, sl]) / np.abs(out["H"][0, sl])).max())]
    inc = out[inc_key][:, sl]
    d["inc_estimator"] = inc_key
    d["inc_rel_err_max"] = float(np.max(np.abs(inc / DTAU_TARGET - 1)))
    d["inc_band_mean"] = [float(np.mean(r)) for r in inc]
    d["inc_band_mean_rel_err"] = [float(np.mean(r) / DTAU_TARGET - 1) for r in inc]
    d["inc_fit_slope"] = [float(np.polyfit(w[sl], np.unwrap(np.angle(out["H"][m + 1, sl] / out["H"][m, sl])), 1)[0]) for m in range(3)]
    _, rel = C.delay_error(out["Dtau"], DTAU_TARGET)
    d["Dtau_rel_err_max"] = [None] + [float(np.max(np.abs(rel[m, sl]))) for m in range(1, 4)]
    d["Dtau_f0"] = [float(out["Dtau"][m, int(np.argmin(abs(f - F0)))]) for m in range(4)]
    d["eps_raw_max"] = [float(np.max(np.abs(eps_raw[m, sl]))) for m in range(4)]
    d["eps_trim_max"] = [float(np.max(np.abs(eps_trim[m, sl]))) for m in range(4)]
    d["eps_raw_max_all"] = float(max(d["eps_raw_max"])); d["eps_trim_max_all"] = float(max(d["eps_trim_max"]))
    # decomposition of the trimmed error: best-fit line in f over the band (mean delay error) and the curvature that remains
    curv = []
    for m in range(1, 4):
        p = np.polyfit(f[sl], eps_trim[m, sl], 1)
        curv.append(float(np.max(np.abs(eps_trim[m, sl] - np.polyval(p, f[sl])))))
    d["eps_trim_after_linear_removal_max"] = [0.0] + curv
    d["trim_delta"] = [float(x) for x in delta]
    return d


if __name__ == "__main__":
    res = int(sys.argv[1])
    out, S, raws = build(res)
    f = out["f"]
    cfg = json.loads(str(out["raw_cfg"][0]))
    assert cfg["decay_tol"] <= 1e-7 + 1e-15, "production requires decay_tol <= 1e-7"
    assert not out["hit_max"].any(), "a run hit max_time"
    i0f = int(np.argmin(abs(f - F0)))
    assert abs(f[i0f] - F0) < 1e-9
    eps_raw = C.raw_phase_error(f, out["phi_rel"], DTAU_TARGET, F0)
    delta = C.trim_constants(eps_raw, f, F0)
    eps_trim = eps_raw - delta[:, None]
    Htrim = C.apply_trim(out["H"], delta)
    out.update(eps_phi_raw=eps_raw, eps_phi_trim=eps_trim, H_trim=Htrim, delta=delta, imbalance_db=C.amplitude_imbalance_db(out["H"]))
    _, out["eps_tau_rel"] = C.delay_error(out["Dtau"], DTAU_TARGET)
    out["eps_tau"], _ = C.delay_error(out["Dtau"], DTAU_TARGET)
    summ = dict(res=res, Ns=LENGTHS, N0=N0, dN=DN, f0=F0, dtau_target=DTAU_TARGET, valid_domain=list(VALID), criteria=CRIT, cfg=cfg,
                t_end={N: float(out["t_end"][i]) for i, N in enumerate(LENGTHS)}, wall_s={N: float(out["wall_s"][i]) for i, N in enumerate(LENGTHS)},
                run_provenance={N: json.loads(str(out["raw_prov"][i])) for i, N in enumerate(LENGTHS)})
    mk, allm = masks(out, f)
    out["mask_common"] = allm
    for k, v in mk.items():
        out["mask_" + k] = v
    summ["criteria_pass_counts_in_domain"] = {k: int((v & mk["domain"]).sum()) for k, v in mk.items()}
    summ["n_domain_samples"] = int(mk["domain"].sum())
    summ["bands"] = {}
    # Revision 5: the local-quadratic estimator is the primary one for the band mask; the pre-declared central difference is reported alongside
    for lab, key, tol in (("primary_poly_tol10", "inc_poly", 0.10), ("predeclared_central_tol10", "inc", 0.10), ("poly_tol5", "inc_poly", 0.05),
                          ("poly_tol25", "inc_poly", 0.25), ("central_tol5", "inc", 0.05), ("central_tol25", "inc", 0.25)):
        m2, a2 = masks(out, f, inc_key=key, inc_tol=tol)
        r = C.widest_run(a2, f)
        summ["bands"][lab] = band_stats(out, f, r[0], r[1], delta, eps_raw, eps_trim, inc_key=key) if r else None
        summ.setdefault("runs", {})[lab] = [[float(f[a]), float(f[b])] for a, b in TM.contiguous_runs(a2)]
        if lab == "primary_poly_tol10":
            out["mask_common"] = a2
    dom = mk["domain"]
    summ["inc_noise_rms_central_minus_poly_over_target_valid"] = [float(np.sqrt(np.mean(((out["inc"][m] - out["inc_poly"][m])[dom] / DTAU_TARGET) ** 2))) for m in range(3)]
    # whole-grid diagnostics (not only the common band)
    summ["diag_valid_domain"] = {f"N{N}": dict(T=[float(out["T"][i][mk["domain"]].min()), float(out["T"][i][mk["domain"]].max())], R_max=float(out["R"][i][mk["domain"]].max()),
                                              res_max=float(max(out["res_in"][i][mk["domain"]].max(), out["res_out"][i][mk["domain"]].max())),
                                              defect_absmax=float(np.abs(out["defect"][i][mk["domain"]]).max())) for i, N in enumerate(LENGTHS)}
    # monotonicity of the channel delay with N (absolute tau, full phase) and of the relative delay, inside the primary band / validated domain
    prim = summ["bands"]["primary_poly_tol10"]
    if prim:
        sl = slice(prim["i_low"], prim["i_high"] + 1)
        summ["monotonic_all_samples_in_band"] = bool(np.all(np.diff(out["tau_abs"][:, sl], axis=0) > 0))
    summ["monotonic_all_samples_in_valid_domain"] = bool(np.all(np.diff(out["tau_abs"][:, mk["domain"]], axis=0) > 0))
    # MPB absolute-phase branch consistency: plain unwrap along f vs MPB-referenced, differences must be a constant multiple of 2 pi
    summ["branch_check_max_dev_rad"] = [float(np.max(np.abs(np.diff(np.unwrap(np.angle(out["H"][i])) - out["phi_abs"][i])))) for i in range(4)]
    # cascade check: predict N = 11 from (8, 9) with 3 steps and from (9, 10) with 2 steps
    casc = {}
    for (a, b, c) in ((0, 1, 3), (1, 2, 3), (0, 1, 2)):          # predict channel c from the pair (a, b): T_c = [T_b T_a^-1]^((N_c - N_b)/(N_b - N_a)) T_b
        steps = (LENGTHS[c] - LENGTHS[b]) // (LENGTHS[b] - LENGTHS[a])
        S11p, S21p = C.cascade_predict(out["S11"][a], out["H"][a], out["S11"][b], out["H"][b], steps)
        key = f"{LENGTHS[a]}_{LENGTHS[b]}_to_{LENGTHS[c]}"
        casc[key] = (S11p, S21p, c)
        out[f"casc_S21_{key}"] = S21p; out[f"casc_S11_{key}"] = S11p
    prim_i = (prim["i_low"], prim["i_high"]) if prim else (int(np.where(mk["domain"])[0][0]), int(np.where(mk["domain"])[0][-1]))
    sl = slice(prim_i[0], prim_i[1] + 1)
    summ["cascade"] = {}
    for key, (S11p, S21p, c) in casc.items():
        dphi = np.angle(S21p / out["H"][c]); ratio = np.abs(S21p) / np.abs(out["H"][c]) - 1
        # delay of the prediction vs direct: derivative of the unwrapped phase of predicted/direct
        dtau = GD.delay_from_phase(f, np.unwrap(dphi))
        summ["cascade"][key] = dict(max_abs_phase_diff=float(np.max(np.abs(dphi[sl]))), max_abs_amp_rel_diff=float(np.max(np.abs(ratio[sl]))),
                                    max_abs_delay_diff_ac=float(np.max(np.abs(dtau[sl]))), max_abs_S11_diff=float(np.max(np.abs(S11p - out["S11"][c])[sl])),
                                    max_abs_phase_diff_valid_domain=float(np.max(np.abs(dphi[mk["domain"]]))),
                                    max_abs_amp_rel_diff_valid_domain=float(np.max(np.abs(ratio[mk["domain"]]))),
                                    accepted=bool(np.max(np.abs(dphi[sl])) <= 0.1 and np.max(np.abs(ratio[sl])) <= 0.03), target_N=LENGTHS[c])
    summ["prov_analysis"] = provenance()
    np.savez_compressed(RES / f"4b_res{res}.npz", **out)
    save_json(RES / f"4b_res{res}.json", summ, with_provenance=False)
    for k, v in summ["bands"].items():
        if v:
            print(k, f"[{v['f_low']:.4f},{v['f_high']:.4f}] BW {v['frac_bw']*100:.1f}% inc_err_max {v['inc_rel_err_max']*100:.1f}% imb {v['amp_imbalance_db_max']:.3f} dB "
                  f"eps_raw {v['eps_raw_max_all']:.3f} eps_trim {v['eps_trim_max_all']:.3f}")
        else:
            print(k, None)
    print("cascade", {k: (round(v['max_abs_phase_diff'], 4), round(v['max_abs_amp_rel_diff'], 4)) for k, v in summ["cascade"].items()})
    print("noise", summ["inc_noise_rms_central_minus_poly_over_target_valid"])
