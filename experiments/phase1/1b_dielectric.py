"""Phase 1B: homogeneous nondispersive dielectric (tau = n L / c) and a single slab vs TMM (reference-subtracted form)."""
from dataclasses import replace
import numpy as np
from _common import *  # noqa
from pcttd import fdtd1d as F1, extract as X, io as pio, tmm as T

BAND = (0.16, 0.44)
BASE = F1.Cfg1D()


def hom(n, L, **kw):
    return replace(BASE, eps_bg=n * n, monitors=(0.0, L), z_hi=L + 8.0, **kw)


def analyze_hom(r, n, L):
    H = r["E"][1] / r["E"][0]
    _, ty = X.yee_1d_delay_phase(r["f"], L, n, r["dx"], r["dt"])
    m = F1.delay_metrics(r["f"], H, n * L, ty, BAND)
    m.update(dt=r["dt"], t_end=r["t_end"], hit_max_time=r["hit_max_time"],
             abs_H_dev=float(np.max(np.abs(np.abs(H[X.band_mask(r["f"], *BAND)]) - 1.0))))
    return m


def slab_cfg(n, Ls, res=40, **kw):
    # monitors at 0 and 10 in vacuum; slab of thickness Ls centered between them
    zl = 5.0 - Ls / 2
    kw.setdefault("decay_tol", 1e-12)
    return replace(BASE, res=res, blocks=((zl, zl + Ls, n * n),), **kw)


def analyze_slab(r_slab, r_ref, n, Ls):
    f = r_slab["f"]
    t_meep = r_slab["E"][1] / r_ref["E"][1]                 # vacuum-referenced transmission (same monitor planes)
    rr = (r_slab["E"][0] - r_ref["E"][0]) / r_ref["E"][0]    # reflected / incident at monitor 1 (phase referenced to monitor)
    layers = [(n, Ls)]
    t_tmm = T.vacuum_referenced_t(f, layers)
    _, _, Tt, Rt = T.transmission_reflection(f, layers)
    m = X.band_mask(f, *BAND)
    tau_m, ok = X.group_delay_phase(f, t_meep)
    tau_t_exact = T.tmm_group_delay(f, layers, referenced=True)        # analytic-derivative TMM delay
    tau_t, _ = X.group_delay_phase(f, t_tmm)                           # SAME estimator and frequency samples as the FDTD data
    return dict(
        estimator_truncation=float(np.max(np.abs(tau_t[m] - tau_t_exact[m]))),
        max_abs_dt=float(np.max(np.abs(t_meep[m] - t_tmm[m]))),
        max_abs_dT=float(np.max(np.abs(np.abs(t_meep[m]) ** 2 - Tt[m]))),
        max_abs_dR=float(np.max(np.abs(np.abs(rr[m]) ** 2 - Rt[m]))),
        max_abs_R_plus_T_minus_1=float(np.max(np.abs(np.abs(rr[m]) ** 2 + np.abs(t_meep[m]) ** 2 - 1.0))),
        max_abs_excess_delay_err=float(np.max(np.abs(tau_m[m] - tau_t[m]))),
        max_excess_delay_tmm=float(np.max(np.abs(tau_t_exact[m]))), unwrap_ok=bool(ok), t_end=r_slab["t_end"],
        mean_excess_delay_meep=float(np.mean(tau_m[m])), mean_excess_delay_tmm=float(np.mean(tau_t[m])),
    )


if __name__ == "__main__":
    out = dict(config=cfg_dict(BASE), band=BAND)
    ns = [1.5, 2.0, float(np.sqrt(12.0))]
    Ls = [2.5, 5.0, 10.0, 20.0, 40.0]
    # 1. delay vs length for each index
    tasks = [(n, L) for n in ns for L in Ls]
    rs = pool_map(F1.run_1d, [hom(n, L) for n, L in tasks])
    out["delay_vs_length"] = []
    for n in ns:
        rows = [dict(n=n, L=L, **analyze_hom(r, n, L)) for (nn, L), r in zip(tasks, rs) if nn == n]
        slope, icpt, resid = X.fit_delay_vs_length([r["L"] for r in rows], [r["mean_tau"] for r in rows])
        out["delay_vs_length"].append(dict(n=n, rows=rows, slope=slope, intercept=icpt, max_residual=resid,
                                           expected_slope=n, slope_over_n=slope / n))
    # 2. resolution (n = sqrt(12), L = 10) and observed order
    n = ns[2]
    resvals = [10, 20, 40, 80, 160]
    rs = pool_map(F1.run_1d, [replace(hom(n, 10.0), res=v) for v in resvals])
    rows = [dict(res=v, **analyze_hom(r, n, 10.0)) for v, r in zip(resvals, rs)]
    h = [1 / v for v in resvals]
    out["resolution_n_sqrt12"] = dict(rows=rows, order_mean_err=X.observed_order(h, [abs(r["mean_err"]) for r in rows])[0],
                                      order_max_err=X.observed_order(h, [r["max_abs_err"] for r in rows])[0])
    # 3. PML thickness and decay tolerance inside the dielectric (n = 2, L = 10)
    n = 2.0
    pv = [0.25, 0.5, 1.0, 2.0, 4.0]
    rs = pool_map(F1.run_1d, [replace(hom(n, 10.0), dpml=v) for v in pv])
    out["pml_n2"] = [dict(dpml=v, **analyze_hom(r, n, 10.0)) for v, r in zip(pv, rs)]
    dv = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
    rs = pool_map(F1.run_1d, [replace(hom(n, 10.0), decay_tol=v) for v in dv])
    out["decay_n2"] = [dict(decay_tol=v, **analyze_hom(r, n, 10.0)) for v, r in zip(dv, rs)]
    # 4. slab with interfaces vs TMM (n = 2 and sqrt(12), thickness 5), resolution study
    out["slab_vs_tmm"] = []
    for n_s in (2.0, float(np.sqrt(12.0))):
        for res in (20, 40, 80, 160):
            r_s, r_ref = pool_map(F1.run_1d, [slab_cfg(n_s, 5.0, res), replace(BASE, res=res, decay_tol=1e-12)])
            if res in (80, 160):
                np.savez_compressed(RES / f"1b_slab_curves_n{n_s:.3f}_res{res}.npz", f=r_s["f"], E_slab=r_s["E"], E_ref=r_ref["E"])
            out["slab_vs_tmm"].append(dict(n=n_s, L_slab=5.0, res=res, **analyze_slab(r_s, r_ref, n_s, 5.0)))
    sl = out["slab_vs_tmm"]
    for n_s in (2.0, float(np.sqrt(12.0))):
        rows = [r for r in sl if r["n"] == n_s]
        out[f"slab_order_n{n_s:.3f}"] = dict(
            order_dt=X.observed_order([1 / r["res"] for r in rows], [r["max_abs_dt"] for r in rows])[0],
            order_delay=X.observed_order([1 / r["res"] for r in rows], [r["max_abs_excess_delay_err"] for r in rows])[0])
    out["provenance"] = pio.provenance()
    pio.save_json(RES / "1b_dielectric.json", out, with_provenance=False)
    print("1B done")
    for d in out["delay_vs_length"]:
        print(" n=%.4f slope/n=%.6f intercept=%.2e max_resid=%.2e" % (d["n"], d["slope_over_n"], d["intercept"], d["max_residual"]))
    print(" res order (n=sqrt12): mean %.3f max %.3f" % (out["resolution_n_sqrt12"]["order_mean_err"], out["resolution_n_sqrt12"]["order_max_err"]))
    for r in sl: print(" slab n=%.3f res=%d |dt|=%.2e |dT|=%.2e |dR|=%.2e R+T-1=%.1e delayerr=%.2e est.trunc=%.2e (max tmm delay %.2f)" % (r["n"], r["res"], r["max_abs_dt"], r["max_abs_dT"], r["max_abs_dR"], r["max_abs_R_plus_T_minus_1"], r["max_abs_excess_delay_err"], r["estimator_truncation"], r["max_excess_delay_tmm"]))
