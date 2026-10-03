"""Phase 1C: 1-D Bragg stack, MEEP vs transfer-matrix method (T, R, complex t, excess delay, stop band, R+T=1, Bloch decay)."""
from dataclasses import replace
import numpy as np
from _common import *  # noqa
from scipy.optimize import minimize_scalar
from pcttd import fdtd1d as F1, extract as X, io as pio, tmm as T

N_LO, N_HI = 1.0, float(np.sqrt(12.0))
D_LO, D_HI, F0 = T.quarter_wave_bilayer(N_LO, N_HI)      # period 1; optical thickness of each layer equal
FREQS = tuple(np.round(np.arange(0.10, 0.60 + 1e-9, 0.0004), 6))
BAND = (0.12, 0.58)
DELAY_MASK_T = 1e-3          # compare delays only where T > this (phase of a vanishing t is noise)
BASE = F1.Cfg1D(res=80, dpml=2.0, z_lo=-8.0, z_hi=16.0, z_src=-5.0, monitors=(-2.0, 10.0), freqs=FREQS,
                src_fcen=0.35, src_fwidth=0.45, decay_tol=1e-12)


def layers_for(N):
    return T.bilayer(N_LO, D_LO, N_HI, D_HI, N)


def stack_cfg(N, **kw):
    D = float(N)
    blocks = tuple((p + D_LO, p + 1.0, N_HI ** 2) for p in range(N))
    base = replace(BASE, z_hi=D + 8.0, monitors=(-2.0, D + 2.0), blocks=blocks)
    return replace(base, **kw)


def ref_cfg(cfg):
    return replace(cfg, blocks=())


def edges(f, Tp, level=0.5, center=F0):
    """Stop band as the interval around ``center`` where T < level; edges by linear interpolation."""
    i0 = int(np.argmin(np.abs(f - center)))
    lo = i0
    while lo > 0 and Tp[lo] < level:
        lo -= 1
    hi = i0
    while hi < len(f) - 1 and Tp[hi] < level:
        hi += 1
    def interp(a, b):
        return float(f[a] + (level - Tp[a]) * (f[b] - f[a]) / (Tp[b] - Tp[a]))
    return interp(lo, lo + 1), interp(hi - 1, hi)


def compare(r_s, r_ref, N):
    f = r_s["f"]
    D = float(N)
    t_meep = r_s["E"][1] / r_ref["E"][1]
    rr = (r_s["E"][0] - r_ref["E"][0]) / r_ref["E"][0]
    layers = layers_for(N)
    t_tmm = T.vacuum_referenced_t(f, layers)
    _, _, Tt, Rt = T.transmission_reflection(f, layers)
    Tm, Rm = np.abs(t_meep) ** 2, np.abs(rr) ** 2
    m = X.band_mask(f, *BAND)
    big = m & (Tt > DELAY_MASK_T)
    tau_m, ok = X.group_delay_phase(f, np.where(big, t_meep, np.nan))
    tau_t, _ = X.group_delay_phase(f, np.where(big, t_tmm, np.nan))
    tau_exact = T.tmm_group_delay(f, layers, referenced=True)
    ok_ = np.isfinite(tau_m) & np.isfinite(tau_t) & big
    lo_t, hi_t = edges(f, Tt)
    lo_m, hi_m = edges(f, Tm)
    # Effective optical-thickness error: best global frequency scale s with T_meep(f) ~ T_tmm(s f)
    def rms(sc):
        return float(np.sqrt(np.mean((Tm[m] - T.transmission_reflection(f[m] * sc, layers)[2]) ** 2)))
    sres = minimize_scalar(rms, bounds=(0.98, 1.02), method="bounded", options=dict(xatol=1e-9))
    s_fit = float(sres.x)
    Tt_s = T.transmission_reflection(f * s_fit, layers)[2]
    # complex-t error where T is not tiny
    rel = np.abs(t_meep - t_tmm)[big] / np.abs(t_tmm)[big]
    return dict(
        N=N, max_abs_dT=float(np.max(np.abs(Tm - Tt)[m])), max_abs_dR=float(np.max(np.abs(Rm - Rt)[m])),
        max_rel_complex_t=float(rel.max()), median_rel_complex_t=float(np.median(rel)),
        max_abs_R_plus_T_minus_1=float(np.max(np.abs(Rm + Tm - 1.0)[m])),
        max_abs_delay_err=float(np.max(np.abs(tau_m - tau_t)[ok_])), max_tmm_delay=float(np.max(np.abs(tau_exact[big]))),
        estimator_truncation=float(np.max(np.abs(tau_t - tau_exact)[ok_])),
        median_abs_delay_err=float(np.median(np.abs(tau_m - tau_t)[ok_])),
        freq_scale_s=s_fit, freq_scale_minus1=s_fit - 1.0, rms_dT_unscaled=rms(1.0), rms_dT_after_scale=float(sres.fun),
        max_abs_dT_after_scale=float(np.max(np.abs(Tm - Tt_s)[m])),
        edge_lo_tmm=lo_t, edge_hi_tmm=hi_t, edge_lo_meep=lo_m, edge_hi_meep=hi_m,
        edge_err_lo=lo_m - lo_t, edge_err_hi=hi_m - hi_t,
        min_T_tmm_in_band=float(Tt[m].min()), min_T_meep_in_band=float(Tm[m].min()),
        min_T_resolved_ratio=float(Tm[X.band_mask(f, 0.3, 0.34)].mean() / Tt[X.band_mask(f, 0.3, 0.34)].mean()),
        unwrap_ok=bool(ok), t_end=r_s["t_end"], hit_max_time=r_s["hit_max_time"],
    )


def run_pair(cfg):
    return pool_map(F1.run_1d, [cfg, ref_cfg(cfg)])


def self_compare(r_s, r_ref, r_s0, r_ref0, N):
    """Run vs. the reference setting of the same study (isolates the studied parameter from the TMM discretization error)."""
    f = r_s["f"]
    m = X.band_mask(f, *BAND)
    t, t0 = r_s["E"][1] / r_ref["E"][1], r_s0["E"][1] / r_ref0["E"][1]
    rr, rr0 = (r_s["E"][0] - r_ref["E"][0]) / r_ref["E"][0], (r_s0["E"][0] - r_ref0["E"][0]) / r_ref0["E"][0]
    big = m & (np.abs(t0) ** 2 > DELAY_MASK_T)
    tau, _ = X.group_delay_phase(f, np.where(big, t, np.nan))
    tau0, _ = X.group_delay_phase(f, np.where(big, t0, np.nan))
    g = np.isfinite(tau) & np.isfinite(tau0) & big
    return dict(self_max_abs_dT=float(np.max(np.abs(np.abs(t) ** 2 - np.abs(t0) ** 2)[m])),
                self_max_abs_dR=float(np.max(np.abs(np.abs(rr) ** 2 - np.abs(rr0) ** 2)[m])),
                self_max_abs_delay=float(np.max(np.abs(tau - tau0)[g])) if g.any() else None,
                self_max_rel_complex_t=float(np.max(np.abs(t - t0)[big] / np.abs(t0)[big])))


def study(param, values, make, N=8, ref=None):
    """``ref``: index of the most-converged setting (default: last value) used for the self-comparison."""
    cfgs = []
    for v in values:
        cfgs += [make(v), ref_cfg(make(v))]
    rs = pool_map(F1.run_1d, cfgs)
    k0 = len(values) - 1 if ref is None else ref
    rows = []
    for i, v in enumerate(values):
        row = dict(param=param, value=v, **compare(rs[2 * i], rs[2 * i + 1], N))
        row.update(self_compare(rs[2 * i], rs[2 * i + 1], rs[2 * k0], rs[2 * k0 + 1], N))
        amp = np.abs(rs[2 * i + 1]["E"][0])
        mk = X.band_mask(rs[2 * i + 1]["f"], *BAND)
        row["source_coverage_min_over_max"] = float(amp[mk].min() / amp.max())
        row["self_reference_value"] = values[k0]
        rows.append(row)
    return rows


if __name__ == "__main__":
    f_ = np.array(FREQS)
    tau_max = float(np.max(np.abs(T.tmm_group_delay(f_, layers_for(12), referenced=True))))
    X.assert_sampling_adequate(f_, 1.5 * tau_max)          # a-priori unwrap guard using the TMM-known maximum delay
    out = dict(config=cfg_dict(BASE), band=BAND, n_lo=N_LO, n_hi=N_HI, d_lo=D_LO, d_hi=D_HI, f0=F0,
               analytic_gap_edges_infinite=T.quarter_wave_gap_edges(N_LO, N_HI),
               tmm_max_delay_N12=tau_max, unwrap_guard_df=float(np.max(np.diff(f_))))
    out["base_N8"] = study("N", [8], lambda v: stack_cfg(v), N=8)[0]
    # period count (finite-structure length)
    Ns = [2, 3, 4, 6, 8, 10, 12]
    rows = []
    for Nv in Ns:
        rows += study("N", [Nv], lambda v: stack_cfg(v), N=Nv)
    out["periods"] = rows
    # resolution (subpixel smoothing of non-grid-aligned layers) and observed order
    resv = [20, 40, 80, 160, 320]
    out["resolution"] = study("res", resv, lambda v: stack_cfg(8, res=v))
    h = [1 / v for v in resv]
    out["resolution_orders"] = {k: X.observed_order(h, [r[k] for r in out["resolution"]])[0]
                                for k in ("max_abs_dT", "max_abs_dR", "max_rel_complex_t", "max_abs_delay_err",
                                                  "rms_dT_unscaled", "rms_dT_after_scale")}
    out["resolution_orders"]["abs_freq_scale_minus1"] = X.observed_order(h, [abs(r["freq_scale_minus1"]) for r in out["resolution"]])[0]
    out["resolution_orders"]["abs_edge_err_hi"] = X.observed_order(h, [abs(r["edge_err_hi"]) for r in out["resolution"]])[0]
    out["resolution_orders"]["abs_edge_err_lo"] = X.observed_order(h, [abs(r["edge_err_lo"]) for r in out["resolution"]])[0]
    out["pml"] = study("dpml", [0.25, 0.5, 1.0, 2.0, 4.0], lambda v: stack_cfg(8, dpml=v))
    out["decay"] = study("decay_tol", [1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14], lambda v: stack_cfg(8, decay_tol=v))
    out["source_width"] = study("src_fwidth", [0.15, 0.25, 0.30, 0.35, 0.45, 0.6, 0.9], lambda v: stack_cfg(8, src_fwidth=v))
    out["placement_left"] = study("monitor A distance to stack", [0.5, 1.0, 2.0, 4.0], lambda v: stack_cfg(8, monitors=(-v, 10.0)))
    out["placement_right"] = study("monitor B distance to stack", [0.5, 1.0, 2.0, 4.0], lambda v: stack_cfg(8, monitors=(-2.0, 8.0 + v), z_hi=8.0 + v + 6.0))
    out["placement_source"] = study("source distance to stack", [1.0, 2.0, 5.0, 8.0], lambda v: stack_cfg(8, z_src=-v, z_lo=-v - 3.0))
    # Bloch decay in the stop band: MEEP slope of ln T vs N against the analytic infinite-lattice rate
    fb = [0.215, 0.23, 0.25]
    Nb = [2, 3, 4, 5, 6]
    cfgs = []
    for Nv in Nb:
        cfgs += [stack_cfg(Nv), ref_cfg(stack_cfg(Nv))]
    rs = pool_map(F1.run_1d, cfgs)
    bl = []
    for fv in fb:
        lnT = []
        for i, Nv in enumerate(Nb):
            f = rs[2 * i]["f"]
            j = int(np.argmin(np.abs(f - fv)))
            lnT.append(np.log(np.abs(rs[2 * i]["E"][1][j] / rs[2 * i + 1]["E"][1][j]) ** 2))
        slope = np.polyfit(Nb, lnT, 1)[0]
        rate = float(T.decay_rate_per_period(np.array([fv]), N_LO, D_LO, N_HI, D_HI, 1.0)[0])
        bl.append(dict(f=fv, meep_decay_per_period=float(-slope), analytic_bloch_decay=rate, lnT=[float(x) for x in lnT]))
    out["bloch_decay"] = dict(N=Nb, rows=bl)
    # curves for figures (saved data only)
    curve_cfgs = [(8, 160), (4, 80), (12, 80)]
    cs = []
    for Nv, rv in curve_cfgs:
        cs += [stack_cfg(Nv, res=rv), ref_cfg(stack_cfg(Nv, res=rv))]
    rs2 = pool_map(F1.run_1d, cs)
    np.savez_compressed(RES / "1c_curves.npz", f=rs2[0]["f"],
                        **{f"E_stack_N{Nv}_res{rv}": rs2[2 * i]["E"] for i, (Nv, rv) in enumerate(curve_cfgs)},
                        **{f"E_ref_N{Nv}_res{rv}": rs2[2 * i + 1]["E"] for i, (Nv, rv) in enumerate(curve_cfgs)})
    out["provenance"] = pio.provenance()
    pio.save_json(RES / "1c_bragg.json", out, with_provenance=False)
    b = out["base_N8"]
    print("1C done. base N=8:", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in b.items()})
    print("orders", out["resolution_orders"])
