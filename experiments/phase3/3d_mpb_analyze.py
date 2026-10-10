"""Phase 3D: MPB-only analysis of D0/D1/D2 (guided band 8) and the conventional R reference: n_g(f), kappa(f), band identity/isolation,
k-grid and resolution sensitivity, MPB-predicted usable bands (same +-10/5/25% rule on n_g, NO transmission information), GDD, required length.
Writes results/phase3/3d_mpb_analysis.{npz,json}.  MPB-only quantities; the FDTD-based usable bands are in 3e_analyze.py."""
from _common import RES, np
import json
from scipy.interpolate import CubicSpline
from pcttd import ttdmetrics as TM
from pcttd.io import save_json

VALID = (0.231, 0.268)
BAND = 8
FG = np.round(np.arange(0.220, 0.275 + 1e-9, 0.0005), 6)        # same grid as the FDTD runs
z = np.load(RES / "3a_mpb_designs.npz")
z48 = np.load(RES / "3c_mpb_res48.npz") if (RES / "3c_mpb_res48.npz").exists() else None
CYC = 25


def branch_arrays(k, freq, vg, band=BAND, sign=-1):
    keep = sign * vg[:, band - 1] > 2e-3
    f = freq[keep, band - 1]; ng = 1.0 / np.abs(vg[keep, band - 1]); kk = k[keep]
    o = np.argsort(f)
    return f[o], ng[o], kk[o]


def curve(k, freq, vg, fgrid, sign=-1):
    f, ng, kk = branch_arrays(k, freq, vg, sign=sign)
    inside = (fgrid >= f.min()) & (fgrid <= f.max())
    cs = CubicSpline(f, ng)
    n = np.where(inside, cs(np.clip(fgrid, f.min(), f.max())), np.nan)
    return n, inside, cs


def band_report(name, ng, inside, f, dom_mask):
    r = {}
    for tol in (0.05, 0.10, 0.25):
        for lab, m in (("full", inside), ("validated", inside & dom_mask)):
            ub = TM.usable_band(f, ng, m, tol=tol)
            if ub is not None:
                ub["gbp"] = TM.delay_bandwidth_product(ub["frac_bw"], ub["median_tau"])
                g = TM.gdd(f, np.nan_to_num(ng, nan=0.0))
                seg = g[ub["i_low"]:ub["i_high"] + 1]
                ub["gdd_center"] = float(np.interp(ub["f_center"], f, g)); ub["gdd_absmax"] = float(np.abs(seg).max())
                ub["target_cycles"] = CYC; tt = TM.target_delay_from_cycles(CYC, ub["f_center"])
                ub["target_delay_ac"] = tt; ub["L_periods"] = TM.required_length_periods(tt, ub["median_tau"])
                ub["L_wavelengths"] = ub["L_periods"] * ub["f_center"]
            r[f"tol{int(tol*100)}_{lab}"] = ub
    return r


out, arr = {}, {"f": FG}
dom = (FG >= VALID[0]) & (FG <= VALID[1])
for name in ("D0", "D1", "D2"):
    k, fr, vg = z[name + "_k"], z[name + "_freq"], z[name + "_vg"]
    ng, inside, cs = curve(k, fr, vg, FG)
    fb, ngb, kb = branch_arrays(k, fr, vg)
    kap = np.where(inside, CubicSpline(fb[::-1] if False else fb, kb)(np.clip(FG, fb.min(), fb.max())), np.nan)
    arr[name + "_ng"] = ng; arr[name + "_kappa"] = kap; arr[name + "_inband"] = inside
    arr[name + "_gdd"] = np.where(inside, TM.gdd(FG, np.nan_to_num(ng, nan=0.0)), np.nan)
    d = dict(band_edge_low=float(fr[:, BAND - 1].min()), band_top=float(fr[:, BAND - 1].max()), ng_min=float(ngb.min()), f_at_ng_min=float(fb[np.argmin(ngb)]))
    # band identity / isolation: monotone decreasing guided band (no k-turnover except the flat ends), gap to neighbors within the even class
    b8 = fr[:, BAND - 1]; d["monotone_decreasing_vg_lt0_count"] = int((vg[:, BAND - 1] < -2e-3).sum()); d["nk"] = int(len(k))
    d["min_gap_to_band7"] = float((b8 - fr[:, BAND - 2]).min()); d["min_gap_to_band9"] = float((fr[:, BAND] - b8).min())
    # k-grid sensitivity: HF n_g from the nk = 101 grid vs the n_g obtained by differentiating the band on a 51-point subgrid (spline)
    sub = slice(0, None, 2)
    kk2, ff2 = k[sub], fr[sub, BAND - 1]; vg_sp = CubicSpline(kk2, ff2).derivative()(k)
    keep = (vg[:, BAND - 1] < -2e-3)
    ref = 1.0 / np.abs(vg[keep, BAND - 1]); alt = 1.0 / np.abs(vg_sp[keep]); fm = fr[keep, BAND - 1]
    m = (fm >= VALID[0]) & (fm <= VALID[1])
    d["kgrid51_spline_vs_HF_ng_maxrel_validated"] = float(np.max(np.abs(alt[m] / ref[m] - 1)))
    if z48 is not None and name in ("D1", "D2"):
        n48, in48, _ = curve(z48[name + "_k"], z48[name + "_freq"], z48[name + "_vg"], FG)
        both = inside & in48 & dom
        d["res48_vs_res32_ng_rel_validated"] = dict(median=float(np.nanmedian(n48[both] / ng[both] - 1)), maxabs=float(np.nanmax(np.abs(n48[both] / ng[both] - 1))))
        d["res48_band_edge"] = float(z48[name + "_freq"][:, BAND - 1].min())
        arr[name + "_ng_res48"] = n48
        d["res48_worst_f"] = float(FG[both][np.nanargmax(np.abs(n48[both] / ng[both] - 1))])
    d["bands"] = band_report(name, ng, inside, FG, dom)
    ub = d["bands"]["tol10_validated"]
    if z48 is not None and name in ("D1", "D2") and ub:
        m = np.zeros_like(dom); m[ub["i_low"]:ub["i_high"] + 1] = True
        d["res48_vs_res32_ng_rel_in_mpb_tol10_validated_band"] = float(np.nanmax(np.abs(n48[m] / ng[m] - 1)))
    out[name] = d

# conventional reference R (square lattice, TM, r/a = 0.19, M = 6, res 24, nk 61; band selected in 3a)
kR, fR, vR = z["R_k"], z["R_freq"], z["R_vg"]
ordr = np.argsort(fR); fR, vR, kR = fR[ordr], vR[ordr], kR[ordr]
ngR = 1.0 / vR; csR = CubicSpline(fR, ngR)
FR = np.round(np.arange(fR.min(), fR.max(), 0.0005), 6)
ngRg = csR(FR); insideR = np.ones_like(FR, bool)
arr["R_f"] = FR; arr["R_ng"] = ngRg; arr["R_gdd"] = TM.gdd(FR, ngRg); arr["R_strip"] = np.interp(FR, fR, z["R_strip"][ordr])
dR = dict(note="incomplete reference comparison - MPB only (res 24, nk = 61, M = 6, no FDTD, no interface/transmission information)",
          f_range=[float(fR.min()), float(fR.max())], ng_min=float(ngR.min()), edge_energy_max=float(z["R_edge"].max()),
          strip_energy_min=float(z["R_strip"].min()))
dR["bands"] = band_report("R", ngRg, insideR, FR, np.ones_like(FR, bool))
out["R"] = dR
np.savez(RES / "3d_mpb_analysis.npz", **arr)
save_json(RES / "3d_mpb_analysis.json", out)
def line(n, b):
    return f"{n}: [{b['f_low']:.4f}, {b['f_high']:.4f}] BW {b['frac_bw']*100:.1f}% ng {b['median_tau']:.2f} ripple {b['ripple']*100:.1f}% gdd_c {b['gdd_center']:.2f} gbp {b['gbp']:.3f} L {b['L_periods']:.1f}" if b else f"{n}: none"
for n in out:
    print(n, {k: v for k, v in out[n].items() if k != "bands"})
    for lab, b in out[n]["bands"].items():
        print("  ", line(lab, b))
