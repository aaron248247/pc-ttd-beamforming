"""Phase 3E: modal S-parameters of the engineered devices (lead = plain W1), multi-length delay extraction, usable bands, GDD.

usage: python 3e_analyze.py <design D0|D1|D2> <res> [N ...]  (reads results/phase3/raw/<design>_res<res>_N<N>.npz; requires 3d_mpb_analysis.npz)
-> results/phase3/3e_<design>_res<res>.{npz,json}

Definitions (a = c = 1; MEEP_SIGN = +1):
  theta_pp(f)      section phase advance per period.  PRIMARY: slope of the unwrapped S21 phase versus N (3-length linear fit of arg[S21(N)/S21(N0)],
                   2 pi branch fixed by the MPB device-section Bloch phase).  CROSS-CHECK: exact wave-cascade (T-matrix) phases of the pairs.
  tau_pp = d theta_pp / d omega (central difference; local-quadratic variant as sensitivity); n_g = tau_pp in these units.
  GDD_pp = d tau_pp / d omega.  Because the central tau is differentiated twice from FDTD phase, GDD uses the local-quadratic (+-3 sample) tau.
"""
from _common import RES, RAW, ROOT, np
import json, sys
from pcttd import guided as GD, w1analysis as WA, ttdmetrics as TM
from pcttd.io import save_json, provenance

VALID = (0.231, 0.268)
BAND = 8
CYC = 25
P2 = ROOT / "results" / "phase2"


def lead_profiles(f):
    z = np.load(P2 / "2c_mode_profiles.npz")
    if not np.array_equal(z["f"], f):
        raise RuntimeError("frequency grid differs from the cached Phase 2 lead-mode profiles")
    return z["y"], z["Ey"], z["Hz"], float(z["H"])


def analyze(design, res, Ns):
    mp3 = np.load(RES / "3d_mpb_analysis.npz")
    raws = {N: dict(np.load(RAW / f"{design}_res{res}_N{N}.npz")) for N in Ns}
    f = raws[Ns[0]]["f"]
    assert np.allclose(f, mp3["f"])
    y, Ey, Hz, H = lead_profiles(f)
    S = {N: WA.modal_sparams(raws[N], y, Ey, Hz, H) for N in Ns}
    S21 = {N: S[N]["S21"] for N in Ns}; S11 = {N: S[N]["S11"] for N in Ns}
    kap = mp3[design + "_kappa"]; ng_mpb = mp3[design + "_ng"]; inband = mp3[design + "_inband"]
    ref_pp = -2 * np.pi * np.where(inband, kap, 0.0)
    out = dict(f=f, ng_mpb=ng_mpb, kappa_mpb=kap, inband=inband, gdd_mpb=mp3[design + "_gdd"])
    for N in Ns:
        out[f"S21_N{N}"] = S21[N]; out[f"S11_N{N}"] = S11[N]; out[f"T_N{N}"] = np.abs(S21[N]) ** 2; out[f"R_N{N}"] = np.abs(S11[N]) ** 2
        out[f"res_in_N{N}"] = S[N]["res_in"]; out[f"res_out_N{N}"] = S[N]["res_out"]; out[f"about_N{N}"] = S[N]["ab_out"]
        out[f"defect_N{N}"] = WA.unitarity_defect(S11[N], S21[N])
    N0 = Ns[0]
    P = [np.zeros_like(f)]
    for N in Ns[1:]:
        P.append(WA.unwrap_to_reference(np.angle(S21[N] / S21[N0]), ref_pp * (N - N0)))
    sl, ic, rms = WA.linear_length_fit(Ns, np.array(P))
    out["theta_fit"], out["fit_intercept"], out["fit_rms"] = sl, ic, rms
    out["tau_fit"] = WA.ng_from_theta(f, sl, "central"); out["tau_fit_poly"] = WA.ng_from_theta(f, sl, "poly", half=3)
    for a in range(len(Ns)):
        for b in range(a + 1, len(Ns)):
            th, mod = WA.cascade_theta(S11, S21, Ns[a], Ns[b])
            thu = WA.unwrap_to_reference(th, ref_pp * (Ns[b] - Ns[a])) / (Ns[b] - Ns[a])
            out[f"theta_cas_{Ns[a]}_{Ns[b]}"] = thu; out[f"eigmod_{Ns[a]}_{Ns[b]}"] = mod
            out[f"tau_cas_{Ns[a]}_{Ns[b]}"] = WA.ng_from_theta(f, thu, "central")
            out[f"tau_cas_poly_{Ns[a]}_{Ns[b]}"] = WA.ng_from_theta(f, thu, "poly", half=3)
    # pairwise S21-phase slopes (no de-embedding of interface multiple reflections), for the stability-with-length comparison
    for a in range(len(Ns)):
        for b in range(a + 1, len(Ns)):
            ph = WA.unwrap_to_reference(np.angle(S21[Ns[b]] / S21[Ns[a]]), ref_pp * (Ns[b] - Ns[a])) / (Ns[b] - Ns[a])
            out[f"tau_s21_{Ns[a]}_{Ns[b]}"] = WA.ng_from_theta(f, ph, "central")
    out["Tmin"] = np.min([out[f"T_N{N}"] for N in Ns], axis=0)
    out["Rmax"] = np.max([out[f"R_N{N}"] for N in Ns], axis=0)
    out["gdd_fit"] = GD.local_poly_derivative(TM.omega(f), out["tau_fit_poly"], half=3)
    return out, raws, Ns


def band_metrics(out, Ns, ub, f):
    sel = slice(ub["i_low"], ub["i_high"] + 1)
    d = dict(ub)
    tau = out["tau_fit"]; w = TM.omega(f)
    d["T_min"] = float(out["Tmin"][sel].min()); d["R_max"] = float(out["Rmax"][sel].max())
    d["res_max"] = float(max(np.max(out[f"res_in_N{N}"][sel]) for N in Ns) if Ns else 0.0)
    d["res_max_out"] = float(max(np.max(out[f"res_out_N{N}"][sel]) for N in Ns))
    d["res_median"] = float(np.median([np.median(out[f"res_in_N{N}"][sel]) for N in Ns] + [np.median(out[f"res_out_N{N}"][sel]) for N in Ns]))
    d["unitarity_defect_absmax"] = float(max(np.max(np.abs(out[f"defect_N{N}"][sel])) for N in Ns))
    d["fit_rms_rad_max"] = float(out["fit_rms"][sel].max())
    d["gbp"] = TM.delay_bandwidth_product(d["frac_bw"], d["median_tau"])
    g = out["gdd_fit"]
    d["gdd_center"] = float(np.interp(ub["f_center"], f, g)); d["gdd_absmax"] = float(np.abs(g[sel]).max())
    d["gdd_band_avg_linfit"] = float(np.polyfit(w[sel], out["tau_fit_poly"][sel], 1)[0])
    d["gdd_mpb_center"] = float(np.nan_to_num(np.interp(ub["f_center"], f, np.nan_to_num(out["gdd_mpb"], nan=np.nan))))
    d["target_delay_ac"] = TM.target_delay_from_cycles(CYC, ub["f_center"])
    d["L_periods"] = TM.required_length_periods(d["target_delay_ac"], d["median_tau"]); d["L_wavelengths"] = d["L_periods"] * ub["f_center"]
    # agreement with MPB and across lengths inside the band
    d["tau_vs_mpb_rel_median"] = float(np.median(tau[sel] / out["ng_mpb"][sel] - 1)); d["tau_vs_mpb_rel_maxabs"] = float(np.max(np.abs(tau[sel] / out["ng_mpb"][sel] - 1)))
    keys = [k for k in out if k.startswith("tau_cas_") and not k.startswith("tau_cas_poly")]
    d["cascade_vs_fit_rel_maxabs"] = {k[8:]: float(np.max(np.abs(out[k][sel] / tau[sel] - 1))) for k in keys}
    d["cascade_vs_fit_rel_median"] = {k[8:]: float(np.median(out[k][sel] / tau[sel] - 1)) for k in keys}
    d["tau_ripple_cascade_pairs"] = {k[8:]: float(TM.delay_ripple(out[k][sel])) for k in keys}
    return d


if __name__ == "__main__":
    design, res = sys.argv[1], int(sys.argv[2])
    Ns = [int(x) for x in sys.argv[3:]] or [4, 8, 12]
    out, raws, Ns = analyze(design, res, Ns)
    f = out["f"]
    dom = (f >= VALID[0] - 1e-9) & (f <= VALID[1] + 1e-9)
    ok_base = (out["Tmin"] >= 0.5) & out["inband"]
    summ = dict(design=design, res=res, Ns=Ns, cfg=json.loads(str(raws[Ns[0]]["cfg_json"])), t_end={N: float(raws[N]["t_end"]) for N in Ns},
                wall_s={N: float(raws[N]["wall_s"]) for N in Ns}, hit_max_time={N: bool(raws[N]["hit_max_time"]) for N in Ns},
                run_provenance={N: json.loads(str(raws[N]["prov_json"])) for N in Ns}, valid_domain=list(VALID), target_cycles=CYC)
    assert summ["cfg"]["decay_tol"] <= 1e-7 + 1e-15, "production requires decay_tol <= 1e-7"
    summ["bands"] = {}
    for variant in ("tau_fit", "tau_fit_poly"):
        for lab, dm in (("validated", dom), ("full_0p22_0p275", np.ones_like(dom))):
            for tol in (0.05, 0.10, 0.25):
                ub = TM.usable_band(f, out[variant], ok_base & dm, tol=tol)
                summ["bands"][f"{variant}|{lab}|tol{int(round(tol*100))}"] = band_metrics(out, Ns, ub, f) if ub else None
    # per-N global diagnostics inside the validated domain and in the whole device band
    summ["diag"] = {}
    for N in Ns:
        m = dom & out["inband"]
        summ["diag"][f"N{N}"] = dict(T_range_valid_inband=[float(out[f"T_N{N}"][m].min()), float(out[f"T_N{N}"][m].max())] if m.any() else None,
                                    T_max_all=float(out[f"T_N{N}"].max()), R_max_valid_inband=float(out[f"R_N{N}"][m].max()) if m.any() else None,
                                    defect_absmax_valid_inband=float(np.abs(out[f"defect_N{N}"][m]).max()) if m.any() else None,
                                    res_in_max_valid_inband=float(out[f"res_in_N{N}"][m].max()) if m.any() else None,
                                    res_out_max_valid_inband=float(out[f"res_out_N{N}"][m].max()) if m.any() else None)
    summ["prov_analysis"] = provenance()
    np.savez_compressed(RES / f"3e_{design}_res{res}.npz", **out)
    save_json(RES / f"3e_{design}_res{res}.json", summ, with_provenance=False)
    for k, v in summ["bands"].items():
        if v and k.startswith("tau_fit|"):
            print(k, f"[{v['f_low']:.4f},{v['f_high']:.4f}] BW {v['frac_bw']*100:.1f}% tau {v['median_tau']:.2f} ripple {v['ripple']*100:.1f}% Tmin {v['T_min']:.2f} Rmax {v['R_max']:.3f} gdd_c {v['gdd_center']:.1f}")
        elif k.startswith("tau_fit|"):
            print(k, None)
