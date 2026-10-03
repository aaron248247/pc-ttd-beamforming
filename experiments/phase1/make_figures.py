"""Phase 1 figures. Reads ONLY saved results (results/phase1/*.json, *.npz); runs no simulation."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _common import *  # noqa
from pcttd import extract as X, io as pio, tmm as T

FIG = ROOT / "figures" / "phase1"
FIG.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, AQUA, MAG, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#0b0b0b", "#898781"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e1",
                     "grid.linewidth": 0.6, "figure.dpi": 130, "axes.axisbelow": True})
BAND = (0.16, 0.44)


def loglog_with_guide(ax, x, y, order, label, color, marker="o"):
    ax.loglog(x, y, marker + "-", color=color, lw=1.5, ms=4, label=label)
    x0, y0 = x[len(x) // 2], y[len(y) // 2]
    xs = np.array([x[0], x[-1]])
    ax.loglog(xs, y0 * (xs / x0) ** order, ":", color=MUTED, lw=1)


def fig_vacuum_dielectric():
    A = pio.load_json(RES / "1a_vacuum.json")
    B = pio.load_json(RES / "1b_dielectric.json")
    z = np.load(RES / "1a_base_curves.npz")
    f, E, dt, dx = z["f"], z["E"], float(z["dt"]), float(z["dx"])
    H = E[1] / E[0]
    tau, _ = X.group_delay_phase(f, H)
    _, pred = X.yee_1d_delay_phase(f, 10.0, 1.0, dx, dt)
    m = X.band_mask(f, *BAND)
    fig, axs = plt.subplots(2, 3, figsize=(11, 6))
    ax = axs[0, 0]
    ax.plot(f[m], tau[m] - 10.0, color=BLUE, lw=1.6, label="FDTD (res 40)")
    ax.plot(f[m], pred[m] - 10.0, "--", color=ORANGE, lw=1.2, label="exact Yee-grid prediction")
    ax.set_xlabel("f = a/λ"); ax.set_ylabel("τ$_g$ − L/c  (a/c)"); ax.set_title("1A vacuum, L = 10: residual vs ideal L/c", fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[0, 1]
    rows = A["resolution"]["rows"]
    res = np.array([r["value"] for r in rows], float)
    loglog_with_guide(ax, 1 / res, [abs(r["mean_err"]) for r in rows], 2, "|mean error| vs L/c", BLUE)
    ax.loglog(1 / res, [r["max_abs_resid_vs_yee"] for r in rows], "s-", color=ORANGE, lw=1.2, ms=3, label="|FDTD − Yee prediction|")
    ax.set_xlabel("grid spacing 1/res"); ax.set_ylabel("delay error (a/c)")
    ax.set_title("order %.2f (measured); dotted = slope 2" % A["resolution"]["order_mean_err"], fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[0, 2]
    Ls = [r["value"] for r in A["length"]["rows"]]
    ax.plot(Ls, [r["mean_tau"] for r in A["length"]["rows"]], "o", color=BLUE, label="FDTD")
    ax.plot(Ls, Ls, "-", color=MUTED, lw=1, label="L/c")
    ax.set_xlabel("L"); ax.set_ylabel("mean τ$_g$"); ax.set_title("slope %.5f, intercept %.1e" % (A["length"]["slope"], A["length"]["intercept"]), fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[1, 0]
    for d, c in zip(B["delay_vs_length"], (BLUE, ORANGE, AQUA)):
        L = [r["L"] for r in d["rows"]]
        ax.plot(L, [r["mean_tau"] for r in d["rows"]], "o", color=c, label="n = %.3f" % d["n"])
        ax.plot(L, np.array(L) * d["n"], "-", color=c, lw=0.8)
    ax.set_xlabel("L"); ax.set_ylabel("mean τ$_g$"); ax.set_title("1B homogeneous dielectric: τ = nL/c (lines)", fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[1, 1]
    for n_s, c in ((2.0, BLUE), (float(np.sqrt(12)), ORANGE)):
        zz = np.load(RES / f"1b_slab_curves_n{n_s:.3f}_res160.npz")
        ff = zz["f"]; mm = X.band_mask(ff, *BAND)
        t_m = zz["E_slab"][1] / zz["E_ref"][1]
        ax.plot(ff[mm], np.abs(t_m[mm]) ** 2, color=c, lw=1.4, label="FDTD n = %.2f (res 160)" % n_s)
        _, _, Tt, _ = T.transmission_reflection(ff, [(n_s, 5.0)])
        ax.plot(ff[mm], Tt[mm], "--", color="k", lw=0.8)
    ax.set_xlabel("f"); ax.set_ylabel("T"); ax.set_title("slab (thickness 5): FDTD vs TMM (dashed)", fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[1, 2]
    for n_s, c in ((2.0, BLUE), (float(np.sqrt(12)), ORANGE)):
        rows = [r for r in B["slab_vs_tmm"] if abs(r["n"] - n_s) < 1e-6]
        loglog_with_guide(ax, [1 / r["res"] for r in rows], [r["max_abs_dt"] for r in rows], 2, "max |t − t$_{TMM}$|, n = %.2f" % n_s, c)
    ax.set_xlabel("grid spacing 1/res"); ax.set_ylabel("error"); ax.set_title("slab vs TMM: second order (dotted)", fontsize=9); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "p1_vacuum_dielectric.png"); plt.close(fig)


def fig_bragg():
    C = pio.load_json(RES / "1c_bragg.json")
    z = np.load(RES / "1c_curves.npz")
    f = z["f"]
    fig, axs = plt.subplots(2, 3, figsize=(11, 6))
    layers = T.bilayer(C["n_lo"], C["d_lo"], C["n_hi"], C["d_hi"], 8)
    Es, Er = z["E_stack_N8_res160"], z["E_ref_N8_res160"]
    t = Es[1] / Er[1]
    _, _, Tt, Rt = T.transmission_reflection(f, layers)
    lo, hi = C["analytic_gap_edges_infinite"]
    ax = axs[0, 0]
    ax.semilogy(f, np.abs(t) ** 2, color=BLUE, lw=1.3, label="FDTD res 160")
    ax.semilogy(f, Tt, "--", color=ORANGE, lw=1, label="TMM")
    ax.axvspan(lo, hi, color=MUTED, alpha=0.15, lw=0)
    ax.set_ylim(1e-9, 2); ax.set_xlabel("f"); ax.set_ylabel("T"); ax.set_title("N = 8: T (shaded: infinite-lattice gap)", fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[0, 1]
    big = Tt > 1e-3
    tau_m, _ = X.group_delay_phase(f, np.where(big, t, np.nan))
    tau_t, _ = X.group_delay_phase(f, np.where(big, T.vacuum_referenced_t(f, layers), np.nan))
    ax.plot(f, tau_m, color=BLUE, lw=1.3, label="FDTD"); ax.plot(f, tau_t, "--", color=ORANGE, lw=1, label="TMM (same estimator)")
    ax.set_xlabel("f"); ax.set_ylabel("excess delay τ$_g$ − D/c"); ax.set_title("delay where T > 1e-3", fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[0, 2]
    rows = C["resolution"]; h = [1 / r["value"] for r in rows]
    loglog_with_guide(ax, h, [abs(r["freq_scale_minus1"]) for r in rows], 2, "|s − 1| (effective optical-thickness error)", BLUE)
    ax.loglog(h, [abs(r["edge_err_hi"]) for r in rows], "s-", color=ORANGE, lw=1.2, ms=3, label="upper −3 dB edge error")
    ax.loglog(h, [r["rms_dT_unscaled"] for r in rows], "^-", color=AQUA, lw=1.2, ms=3, label="rms ΔT")
    ax.set_xlabel("grid spacing 1/res"); ax.set_title("convergence (N = 8); dotted = slope 2", fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[1, 0]
    P = C["periods"]
    Ns = [r["N"] for r in P]
    ax.plot(Ns, [r["edge_hi_meep"] for r in P], "o", color=BLUE, label="FDTD upper edge (res 80)")
    ax.plot(Ns, [r["edge_hi_tmm"] for r in P], "-", color=ORANGE, label="TMM upper edge")
    ax.plot(Ns, [r["edge_lo_meep"] for r in P], "o", color=BLUE); ax.plot(Ns, [r["edge_lo_tmm"] for r in P], "-", color=ORANGE)
    ax.axhline(hi, color=MUTED, ls=":"); ax.axhline(lo, color=MUTED, ls=":")
    ax.set_xlabel("periods N"); ax.set_ylabel("−3 dB edge (f)"); ax.set_title("finite-length effect → infinite edges (dotted)", fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[1, 1]
    Nb = C["bloch_decay"]["N"]
    for r, c in zip(C["bloch_decay"]["rows"], (BLUE, ORANGE, AQUA)):
        ax.plot(Nb, r["lnT"], "o", color=c, label="f = %.3f: FDTD %.3f vs Bloch %.3f" % (r["f"], r["meep_decay_per_period"], r["analytic_bloch_decay"]))
        ax.plot(Nb, np.polyval([-r["analytic_bloch_decay"], r["lnT"][0] + r["analytic_bloch_decay"] * Nb[0]], Nb), "-", color=c, lw=0.8)
    ax.set_xlabel("periods N"); ax.set_ylabel("ln T"); ax.set_title("stop-band decay per period (lines: Bloch)", fontsize=9); ax.legend(frameon=False, fontsize=6)
    ax = axs[1, 2]
    dec = C["decay"]
    ax.loglog([r["value"] for r in dec], [max(r["max_abs_R_plus_T_minus_1"], 1e-9) for r in dec], "o-", color=BLUE, label="max |R + T − 1|")
    ax.loglog([r["value"] for r in dec], [max(r["self_max_abs_dT"], 1e-9) for r in dec], "s-", color=ORANGE, label="max |ΔT| vs tightest")
    ax.set_xlabel("decay tolerance"); ax.set_title("termination criterion (N = 8)", fontsize=9); ax.legend(frameon=False, fontsize=7); ax.invert_xaxis()
    fig.tight_layout(); fig.savefig(FIG / "p1_bragg.png"); plt.close(fig)


def fig_pc2d():
    D = pio.load_json(RES / "1d_pc2d.json")
    M = pio.load_json(RES / "1d_mpb_directional.json")
    ev = D["mpb_even_gap_res128"]; ev64 = D["mpb_even_gap_res64"]; an = D["mpb_any_gap_res128"]
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    ax = axs[0, 0]
    cols = [BLUE, ORANGE, AQUA, MAG]
    for r, c in zip([x for x in D["length"] if x["N"] in (2, 4, 8, 16)], cols):
        cv = r["curve"]
        ax.semilogy(cv["f"], np.maximum(cv["T"], 1e-9), color=c, lw=1.2, label="N = %d" % r["N"])
    ax.axvspan(ev[0], ev[1], color=MUTED, alpha=0.18, lw=0, label="MPB even-Φ gap (k$_y$=0)")
    ax.axvline(an[1], color=MUTED, ls=":", lw=1)
    ax.set_ylim(1e-8, 2); ax.set_xlabel("f = a/λ"); ax.set_ylabel("T"); ax.set_title("finite crystal transmission (r/a = 0.30, ny = %d)" % D["config"]["ny"], fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[0, 1]
    L = D["length"]; Ns = [r["N"] for r in L]
    ax.plot(Ns, [r["edge_lo_3dB"] for r in L], "o-", color=BLUE, label="−3 dB lower"); ax.plot(Ns, [r["edge_hi_3dB"] for r in L], "o-", color=ORANGE, label="−3 dB upper")
    ax.plot(Ns, [r["edge_lo_20dB"] for r in L], "s--", color=BLUE, ms=3, label="−20 dB lower"); ax.plot(Ns, [r["edge_hi_20dB"] for r in L], "s--", color=ORANGE, ms=3, label="−20 dB upper")
    ax.axhline(ev[0], color=MUTED, ls=":"); ax.axhline(ev[1], color=MUTED, ls=":")
    ax.set_xlabel("columns N"); ax.set_ylabel("stop-band edge (f)"); ax.set_title("edges → MPB band edges (dotted) as N grows", fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[1, 0]
    R_ = D["resolution"]; nyv = np.array([r["ny"] for r in R_], float)
    ax.plot(nyv, [r["edge_lo_20dB"] for r in R_], "s--", color=BLUE, ms=3, label="−20 dB lower"); ax.plot(nyv, [r["edge_hi_20dB"] for r in R_], "s--", color=ORANGE, ms=3, label="−20 dB upper")
    ax.plot(nyv, [r["edge_lo_3dB"] for r in R_], "o-", color=BLUE, label="−3 dB lower"); ax.plot(nyv, [r["edge_hi_3dB"] for r in R_], "o-", color=ORANGE, label="−3 dB upper")
    ax.axhline(ev[0], color=MUTED, ls=":"); ax.axhline(ev[1], color=MUTED, ls=":")
    ax.set_xlabel("ny (pixels across √3 a)"); ax.set_ylabel("edge (f)"); ax.set_title("resolution (N = %d)" % D["config"]["N"], fontsize=9); ax.legend(frameon=False, fontsize=7)
    ax = axs[1, 1]
    for key, c, lab in (("decay", BLUE, "decay tol"), ):
        rr = D[key]
        ax.semilogx([r["value"] for r in rr], [r["max_abs_RplusT_minus_1_all"] for r in rr], "o-", color=c, label="max |R+T−1| (band)")
        ax.semilogx([r["value"] for r in rr], [r["max_abs_RplusT_minus_1_gap"] for r in rr], "s-", color=ORANGE, label="max |R+T−1| (gap)")
    ax.set_yscale("log"); ax.set_xlabel("decay tolerance"); ax.set_title("energy conservation vs termination", fontsize=9); ax.legend(frameon=False, fontsize=7); ax.invert_xaxis()
    fig.tight_layout(); fig.savefig(FIG / "p1_pc2d.png"); plt.close(fig)


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["vd", "bragg", "pc2d"]
    if "vd" in which: fig_vacuum_dielectric()
    if "bragg" in which: fig_bragg()
    if "pc2d" in which: fig_pc2d()
    print("figures written to", FIG)
