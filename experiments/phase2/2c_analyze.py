"""Phase 2C: modal S-parameters, multi-length de-embedding, group index vs MPB, estimator and sampling sensitivity, length-scaling fit.

usage: python 2c_analyze.py <tag> <res> [N ...]   (reads results/phase2/raw/<tag>_res<res>_N<N>.npz)  -> results/phase2/2c_<tag>_res<res>.{npz,json}
Mode profiles (plain-W1 forward band-8 mode at every analysis frequency) are cached in results/phase2/2c_mode_profiles.npz.
"""
from _common import RES, RAW, np
import json, sys
from scipy.interpolate import CubicSpline
from pcttd import w1modes as WM, guided as GD, w1analysis as WA, extract as X
from pcttd.io import save_json, provenance

BAND, M = 8, 6
USEFUL = None


def mode_profiles(f, J=6):
    p = RES / ("2c_mode_profiles.npz" if J == 6 else f"2c_mode_profiles_M{J}.npz")
    if p.exists():
        z = np.load(p)
        if np.array_equal(z["f"], f):
            return z["y"], z["Ey"], z["Hz"], float(z["H"]), z["kappa"]
    b = dict(np.load(RES / "2a_mpb_plain_fine.npz")) if J == 6 else dict(np.load(RES / f"2a_mpb_M{J}.npz"))
    kap = WM.kappa_of_f(b, J + 2, f)
    P = WM.forward_profiles(kap, J + 2, M=J, res=32)
    # sanity: the solver's frequency at kappa equals the requested frequency
    assert np.max(np.abs(P["freq"] - f)) < 2e-4, np.max(np.abs(P["freq"] - f))
    np.savez(p, f=f, y=P["y"], Ey=P["Ey"], Hz=P["Hz"], H=P["H"], kappa=kap, vg=P["vg"])
    return P["y"], P["Ey"], P["Hz"], float(P["H"]), kap


def mpb_ng(f):
    b = dict(np.load(RES / "2a_mpb_plain_fine.npz"))
    keep = b["vg"][:, BAND - 1] < -2e-3
    fr, ng = b["freq"][keep, BAND - 1], 1.0 / np.abs(b["vg"][keep, BAND - 1])
    return CubicSpline(fr[::-1], ng[::-1])(f)


def analyze(tag, res, Ns):
    raws = {N: dict(np.load(RAW / f"{tag}_res{res}_N{N}.npz")) for N in Ns}
    f = raws[Ns[0]]["f"]
    J = int(json.loads(str(raws[Ns[0]]["cfg_json"]))["J"])
    y, Ey, Hz, H, kap = mode_profiles(f, J)
    S = {N: WA.modal_sparams(raws[N], y, Ey, Hz, H) for N in Ns}
    S21 = {N: S[N]["S21"] for N in Ns}; S11 = {N: S[N]["S11"] for N in Ns}
    ng_ref = mpb_ng(f)
    out = dict(f=f, ng_mpb=ng_ref, kappa_mpb=kap)
    summ = dict(tag=tag, res=res, Ns=list(Ns), t_end={N: float(raws[N]["t_end"]) for N in Ns}, wall_s={N: float(raws[N]["wall_s"]) for N in Ns},
                hit_max_time={N: bool(raws[N]["hit_max_time"]) for N in Ns}, cfg=json.loads(str(raws[Ns[0]]["cfg_json"])))
    for N in Ns:
        out[f"S21_N{N}"] = S21[N]; out[f"S11_N{N}"] = S11[N]; out[f"res_in_N{N}"] = S[N]["res_in"]; out[f"res_out_N{N}"] = S[N]["res_out"]
        out[f"about_N{N}"] = S[N]["ab_out"]
    # pairwise cascade phases and group index
    pairs = [(Ns[i], Ns[j]) for i in range(len(Ns)) for j in range(i + 1, len(Ns))]
    theta_ref_pp = -2 * np.pi * kap                           # forward Bloch phase per period in radians (MPB k is in units of 2 pi / a; k_fwd = -kappa)
    ng_pairs = {}
    for a, b in pairs:
        th, mod = WA.cascade_theta(S11, S21, a, b)
        thu = WA.unwrap_to_reference(th, theta_ref_pp * (b - a)) / (b - a)
        out[f"theta_pp_{a}_{b}"] = thu; out[f"eigmod_{a}_{b}"] = mod
        out[f"ng_{a}_{b}"] = WA.ng_from_theta(f, thu, "central"); out[f"ng_poly_{a}_{b}"] = WA.ng_from_theta(f, thu, "poly", half=3)
        ng_pairs[(a, b)] = out[f"ng_{a}_{b}"]
    # direct length fit of the S21 phase (needs >= 3 lengths): phase relative to the shortest length, branch fixed by MPB
    if len(Ns) >= 3:
        N0 = Ns[0]
        P = [np.zeros_like(f)]
        for N in Ns[1:]:
            ph = np.angle(S21[N] / S21[N0])
            P.append(WA.unwrap_to_reference(ph, theta_ref_pp * (N - N0)))
        sl, ic, rms = WA.linear_length_fit(Ns, np.array(P))
        out["fit_slope"], out["fit_intercept"], out["fit_rms"] = sl, ic + sl * N0, rms   # intercept re-referenced to N0: phase offset at the shortest length (0 if the fit is exact)
        out["ng_fit"] = WA.ng_from_theta(f, sl, "central")
    return out, summ, f


if __name__ == "__main__":
    tag, res = sys.argv[1], int(sys.argv[2])
    Ns = [int(x) for x in sys.argv[3:]] or [6, 12, 24]
    out, summ, f = analyze(tag, res, Ns)
    np.savez_compressed(RES / f"2c_{tag}_res{res}.npz", **out)
    ints = json.loads((RES / "2a_mpb_plain.json").read_text())["useful_interval"]
    m_int = (f >= ints[0]) & (f <= ints[1]); m_all = (f > 0.225) & (f < 0.272)
    pk = [k for k in out if k.startswith("ng_") and k[3:4].isdigit()]
    for k in ["ng_fit"] if "ng_fit" in out else []: pk.append(k)
    for k in pk:
        summ[k] = dict(useful=WA.ng_vs_reference(f, out[k], out["ng_mpb"], m_int), all_0p225_0p272=WA.ng_vs_reference(f, out[k], out["ng_mpb"], m_all))
        for kk in summ[k]: summ[k][kk] = {a: b for a, b in summ[k][kk].items() if a not in ("abs", "rel")}
    summ["useful_interval"] = ints
    for N in Ns:
        summ[f"N{N}"] = dict(S21_abs_useful=[float(np.abs(out[f"S21_N{N}"])[m_int].min()), float(np.abs(out[f"S21_N{N}"])[m_int].max())],
                             S11_abs_max_useful=float(np.abs(out[f"S11_N{N}"])[m_int].max()),
                             unitarity_defect_max_useful=float(np.abs(WA.unitarity_defect(out[f"S11_N{N}"], out[f"S21_N{N}"]))[m_int].max()),
                             res_in_range_useful=[float(out[f"res_in_N{N}"][m_int].min()), float(out[f"res_in_N{N}"][m_int].max())],
                             res_out_range_useful=[float(out[f"res_out_N{N}"][m_int].min()), float(out[f"res_out_N{N}"][m_int].max())],
                             backward_at_output_max_useful=float(np.abs(out[f"about_N{N}"])[m_int].max()))
    if "fit_rms" in out:
        summ["length_fit"] = dict(rms_residual_rad_useful_max=float(out["fit_rms"][m_int].max()), intercept_rad_useful=[float(out["fit_intercept"][m_int].min()), float(out["fit_intercept"][m_int].max())])
    save_json(RES / f"2c_{tag}_res{res}.json", summ)
    print(json.dumps({k: v for k, v in summ.items() if k not in ("cfg",)}, indent=1)[:6000])
