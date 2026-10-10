"""Phase 3 figures. Reads ONLY saved results (results/phase3/*.json, *.npz); runs no simulation."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _common import RES, ROOT

FIG = ROOT / "figures" / "phase3"
FIG.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, AQUA, MAG, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#0b0b0b", "#898781"
COL = {"D0": BLUE, "D1": ORANGE, "D2": AQUA}
LAB = {"D0": "D0 (f1=1.0, s=0)", "D1": "D1 (f1=1.2, s=-0.05)", "D2": "D2 (f1=1.2, s=-0.10)"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e1",
                     "grid.linewidth": 0.6, "figure.dpi": 130, "axes.axisbelow": True})
VALID = (0.231, 0.268)
DES = ("D0", "D1", "D2")
Z = np.load(RES / "3a_mpb_designs.npz"); MP = np.load(RES / "3d_mpb_analysis.npz"); MJ = json.loads((RES / "3d_mpb_analysis.json").read_text())
F = {d: np.load(RES / f"3e_{d}_res30.npz") for d in DES}
FJ = {d: json.loads((RES / f"3e_{d}_res30.json").read_text()) for d in DES}
f = MP["f"]
KEY = "tau_fit|validated|tol10"


def valid(ax):
    ax.axvspan(*VALID, color=MUTED, alpha=0.12, lw=0)


def usable(ax, d, color):
    b = FJ[d]["bands"][KEY]
    if b:
        ax.axvspan(b["f_low"], b["f_high"], ymin=0, ymax=0.04, color=color, alpha=0.9, lw=0)


def fig_bands():
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.3))
    for d in DES:
        k, fr = Z[d + "_k"], Z[d + "_freq"][:, 7]
        ax[0].plot(k, fr, color=COL[d], lw=1.6, label=LAB[d])
    ax[0].set_xlabel(r"$\kappa$ ($2\pi/a$)"); ax[0].set_ylabel(r"$f=a/\lambda$"); ax[0].set_ylim(0.20, 0.30); ax[0].legend(frameon=False, fontsize=7)
    ax[0].axhspan(*VALID, color=MUTED, alpha=0.12, lw=0); ax[0].set_title("MPB guided band (Hz-even #8)", fontsize=9)
    for d in DES:
        ax[1].plot(f, MP[d + "_ng"], color=COL[d], lw=1.6, label="MPB " + d)
        pm = np.where((F[d]["Tmin"] >= 0.5) & F[d]["inband"], 1.0, np.nan)
        ax[1].plot(f[::4], (F[d]["tau_fit"] * pm)[::4], ".", color=COL[d], ms=3, alpha=0.7)
    valid(ax[1]); ax[1].set_ylim(3.4, 8); ax[1].set_xlabel("$f$"); ax[1].set_ylabel(r"$n_g$ (a/c per period)")
    ax[1].set_title("MPB lines, FDTD 3-length fit dots; shaded = Phase 2 validated range", fontsize=7.5); ax[1].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "3a_bands_ng.png"); plt.close(fig)


def fig_T():
    fig, ax = plt.subplots(2, 3, figsize=(10.5, 5.2), sharex=True)
    for i, d in enumerate(DES):
        for N, ls in zip(FJ[d]["Ns"], ("-", "--", ":")):
            ax[0, i].plot(f, F[d][f"T_N{N}"], ls, color=COL[d], lw=1.3, label=f"N={N}")
            ax[1, i].semilogy(f, F[d][f"R_N{N}"], ls, color=COL[d], lw=1.3)
        ax[0, i].axhline(0.5, color=INK, lw=0.6); valid(ax[0, i]); valid(ax[1, i]); ax[0, i].set_title(LAB[d], fontsize=9)
        ax[0, i].set_ylim(-0.02, 1.08); ax[0, i].set_ylabel(r"$T=|S_{21}|^2$"); ax[1, i].set_ylabel(r"$R=|S_{11}|^2$"); ax[1, i].set_xlabel("$f$"); ax[1, i].set_ylim(1e-5, 1)
        ax[0, i].legend(frameon=False, fontsize=7, loc="lower right")
        mb = MP[d + "_inband"]
        if (~mb).any():
            ax[0, i].axvspan(f[0], f[mb][0], color=ORANGE, alpha=0.08, lw=0)
    fig.tight_layout(); fig.savefig(FIG / "3b_transmission_reflection.png"); plt.close(fig)


def fig_diag():
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
    for d in DES:
        for N, ls in zip(FJ[d]["Ns"], ("-", "--", ":")):
            ax[0].plot(f, F[d][f"res_out_N{N}"], ls, color=COL[d], lw=1.0)
            ax[1].plot(f, np.abs(F[d][f"defect_N{N}"]), ls, color=COL[d], lw=1.0)
        ax[2].semilogy(f, F[d]["fit_rms"], color=COL[d], lw=1.2)
    ax[0].set_ylabel("projection residual (output plane)"); ax[1].set_ylabel(r"$|1-|S_{11}|^2-|S_{21}|^2|$"); ax[2].set_ylabel("3-length phase-fit rms [rad]")
    for a in ax:
        valid(a); a.set_xlabel("$f$")
    ax[1].axhline(0.02, color=INK, lw=0.6)
    fig.tight_layout(); fig.savefig(FIG / "3c_diagnostics.png"); plt.close(fig)


def fig_delay():
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
    for i, d in enumerate(DES):
        ax[i].plot(f, MP[d + "_ng"], color=INK, lw=1.0, label="MPB")
        pm = np.where((F[d]["Tmin"] >= 0.5) & F[d]["inband"], 1.0, np.nan)      # evanescent / low-T region: delay undefined, not plotted
        ax[i].plot(f, F[d]["tau_fit"] * pm, color=COL[d], lw=1.0, label="FDTD fit, central diff")
        for k, ls in (("4_12", "-"), ("8_12", ":"), ("4_8", "--")):
            ax[i].plot(f, F[d][f"tau_cas_{k}"] * pm, ls, color=COL[d], lw=0.7, alpha=0.7, label="cascade " + k)
        b = FJ[d]["bands"][KEY]
        valid(ax[i])
        if b:
            ax[i].axvspan(b["f_low"], b["f_high"], color=COL[d], alpha=0.12, lw=0)
            ax[i].axhline(b["median_tau"] * 1.1, color=COL[d], lw=0.5); ax[i].axhline(b["median_tau"] * 0.9, color=COL[d], lw=0.5)
        ax[i].set_ylim(3, 7); ax[i].set_xlabel("$f$"); ax[i].set_title(LAB[d], fontsize=9)
    ax[0].set_ylabel(r"$\tau_{pp}$ (a/c per period)"); ax[0].legend(frameon=False, fontsize=6.5)
    fig.suptitle(r"per-period delay; colored band = usable band ($\pm$10%, T$\geq$0.5, validated range), lines = $\pm$10% of its median", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "3d_group_delay.png"); plt.close(fig)


def fig_gdd():
    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.3))
    for d in DES:
        ax[0].plot(f, MP[d + "_gdd"], color=COL[d], lw=1.5, label="MPB " + d)
        ax[1].plot(f, F[d]["gdd_fit"] * np.where((F[d]["Tmin"] >= 0.5) & F[d]["inband"], 1.0, np.nan), color=COL[d], lw=0.9, label="FDTD " + d)
    for a in ax:
        valid(a); a.set_xlabel("$f$"); a.set_ylabel(r"GDD per period $d\tau/d\omega$ [(a/c)$^2$]"); a.legend(frameon=False, fontsize=7)
    ax[0].set_ylim(-12, 6); ax[1].set_ylim(-12, 12)
    ax[1].set_title("FDTD: double derivative of phase, noise-limited", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "3e_gdd.png"); plt.close(fig)


def fig_bands_bars():
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.0), sharey=True)
    for a, src, ttl in ((ax[0], "FDTD", "FDTD usable bands (validated range)"), (ax[1], "MPB", "MPB-only n_g flatness bands (no transmission info)")):
        for i, d in enumerate(DES):
            for j, tol in enumerate((25, 10, 5)):
                b = FJ[d]["bands"][f"tau_fit|validated|tol{tol}"] if src == "FDTD" else MJ[d]["bands"][f"tol{tol}_validated"]
                if b:
                    a.barh(i * 4 + j, b["f_high"] - b["f_low"], left=b["f_low"], color=COL[d], alpha=(0.3, 0.9, 0.55)[j], height=0.8)
                    a.text(b["f_high"] + 0.0005, i * 4 + j, f"{tol}%: {b['frac_bw']*100:.1f}%", va="center", fontsize=6.5)
        a.axvspan(*VALID, color=MUTED, alpha=0.1, lw=0); a.set_xlim(0.225, 0.285); a.set_xlabel("$f$"); a.set_title(ttl, fontsize=8)
    ax[0].set_yticks([1, 5, 9]); ax[0].set_yticklabels(DES)
    fig.tight_layout(); fig.savefig(FIG / "3f_usable_bands.png"); plt.close(fig)


def fig_tradeoff():
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
    for d in DES:
        b = FJ[d]["bands"][KEY]; m = MJ[d]["bands"]["tol10_validated"]
        if b:
            ax[0].plot(b["median_tau"], b["frac_bw"] * 100, "o", color=COL[d], label=d + " FDTD")
            ax[1].plot(b["gbp"], b["L_periods"], "o", color=COL[d])
            ax[2].plot(abs(b["gdd_band_avg_linfit"]), b["ripple"] * 100, "o", color=COL[d])
        if m:
            ax[0].plot(m["median_tau"], m["frac_bw"] * 100, "s", mfc="none", color=COL[d], label=d + " MPB")
            ax[1].plot(m["gbp"], m["L_periods"], "s", mfc="none", color=COL[d]); ax[2].plot(abs(m["gdd_center"]), m["ripple"] * 100, "s", mfc="none", color=COL[d])
    ax[0].set_xlabel(r"median $\tau_{pp}$ (a/c)"); ax[0].set_ylabel("usable fractional BW [%]"); ax[0].legend(frameon=False, fontsize=6.5)
    ax[1].set_xlabel("group-delay-bandwidth product"); ax[1].set_ylabel("periods for 25 carrier cycles")
    ax[2].set_xlabel(r"|GDD| band average (FDTD) / center (MPB) [(a/c)$^2$]"); ax[2].set_ylabel("delay ripple [%]")
    fig.suptitle("tradeoff illustration (three designs; not a Pareto frontier)", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "3g_tradeoff.png"); plt.close(fig)


def fig_R():
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
    fR = MP["R_f"]
    ax[0].plot(fR, MP["R_ng"], color=MAG, lw=1.6, label="R (MPB only)")
    ax[0].plot(f, MP["D0_ng"], color=BLUE, lw=1.4, label="D0 (MPB)")
    ax[0].set_ylim(0, 8); ax[0].set_xlabel("$f$"); ax[0].set_ylabel("$n_g$"); ax[0].legend(frameon=False, fontsize=7)
    x = lambda fr_, fc: (fr_ - fc) / fc
    bR = MJ["R"]["bands"]["tol10_full"]; bD = MJ["D0"]["bands"]["tol10_validated"]
    ax[1].plot(x(fR, bR["f_center"]), MP["R_ng"], color=MAG, lw=1.6); ax[1].plot(x(f, bD["f_center"]), MP["D0_ng"], color=BLUE, lw=1.4)
    ax[1].set_xlim(-0.3, 0.3); ax[1].set_ylim(0, 8); ax[1].set_xlabel(r"$(f-f_c)/f_c$ about each +-10% band center"); ax[1].set_ylabel("$n_g$")
    ax[2].plot(x(fR, bR["f_center"]), MP["R_gdd"], color=MAG, lw=1.6); ax[2].plot(x(f, bD["f_center"]), MP["D0_gdd"], color=BLUE, lw=1.4)
    ax[2].set_xlim(-0.3, 0.3); ax[2].set_ylim(-6, 3); ax[2].set_xlabel(r"$(f-f_c)/f_c$"); ax[2].set_ylabel(r"GDD per period [(a/c)$^2$]")
    fig.suptitle("INCOMPLETE REFERENCE COMPARISON: R is MPB only (res 24, no FDTD); D0 shown as MPB for a like-for-like comparison", fontsize=8, color="#b00020")
    fig.tight_layout(); fig.savefig(FIG / "3h_R_reference_mpb_only.png"); plt.close(fig)


for fn in (fig_bands, fig_T, fig_diag, fig_delay, fig_gdd, fig_bands_bars, fig_tradeoff, fig_R):
    fn()
print("figures written to", FIG)
