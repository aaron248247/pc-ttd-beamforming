"""Phase 2 figures. Reads ONLY saved results (results/phase2/*.json, *.npz); runs no simulation."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _common import RES, ROOT

FIG = ROOT / "figures" / "phase2"
FIG.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, AQUA, MAG, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#0b0b0b", "#898781"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e1",
                     "grid.linewidth": 0.6, "figure.dpi": 130, "axes.axisbelow": True})
J = json.loads((RES / "2a_mpb_plain.json").read_text())
U = J["useful_interval"]
FINAL = np.load(RES / "2c_tol7_res40.npz")        # production run: res 40, decay 1e-7, N = 6, 12, 24
f = FINAL["f"]


def shade(ax):
    ax.axvspan(*U, color=BLUE, alpha=0.08, lw=0)


def fig_mpb():
    b = np.load(RES / "2a_mpb_plain_fine.npz")
    k, fr, vg = b["k"], b["freq"][:, 7], b["vg"][:, 7]
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))
    for j in range(5, 9):
        ax[0].plot(k, b["freq"][:, j], color=MUTED, lw=0.8)
    ax[0].plot(k, fr, color=BLUE, lw=1.8, label="guided band (Hz-even #8)")
    ax[0].set_xlabel(r"$\kappa$ ($2\pi/a$)"); ax[0].set_ylabel(r"$f = a/\lambda$"); ax[0].set_ylim(0.18, 0.34); ax[0].legend(frameon=False, loc="upper right")
    keep = vg < -2e-3
    ng = 1 / np.abs(vg[keep])
    ax[1].semilogy(fr[keep], ng, color=BLUE, lw=1.5)
    ax[1].axvspan(*U, color=BLUE, alpha=0.08, lw=0, label="useful interval (n_g <= 1.1 n_g,min)")
    ax[1].set_xlabel(r"$f$"); ax[1].set_ylabel(r"$n_g$ (MPB, Hellmann-Feynman)"); ax[1].legend(frameon=False, fontsize=7)
    ax[2].plot(fr[keep], ng, color=BLUE, lw=1.5); shade(ax[2]); ax[2].set_ylim(3.5, 6); ax[2].set_xlim(0.215, 0.285)
    ax[2].set_xlabel(r"$f$"); ax[2].set_ylabel(r"$n_g$ (zoom)")
    fig.tight_layout(); fig.savefig(FIG / "2a_mpb_band_ng.png"); plt.close(fig)


def fig_sparams():
    fig, ax = plt.subplots(1, 4, figsize=(12.5, 3.0))
    for N, c in ((6, BLUE), (12, ORANGE), (24, AQUA)):
        ax[0].plot(f, np.abs(FINAL[f"S21_N{N}"]), color=c, lw=1, label=f"N={N}")
        ax[1].semilogy(f, np.abs(FINAL[f"S11_N{N}"]), color=c, lw=1)
        ax[2].semilogy(f, FINAL[f"res_in_N{N}"], color=c, lw=1)
        ax[3].plot(f, 1 - np.abs(FINAL[f"S11_N{N}"]) ** 2 - np.abs(FINAL[f"S21_N{N}"]) ** 2, color=c, lw=1)
    for a, t in zip(ax, [r"$|S_{21}|$", r"$|S_{11}|$", "projection residual (input plane)", r"$1-|S_{11}|^2-|S_{21}|^2$"]):
        shade(a); a.set_xlabel(r"$f$"); a.set_ylabel(t)
    ax[0].legend(frameon=False, fontsize=7); ax[0].set_ylim(0.95, 1.05)
    fig.tight_layout(); fig.savefig(FIG / "2b_sparams.png"); plt.close(fig)


def fig_ng():
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    ax[0].plot(f, FINAL["ng_mpb"], color=INK, lw=1.8, label="MPB (Hellmann-Feynman)")
    ax[0].plot(f, FINAL["ng_fit"], color=BLUE, lw=1, label="FDTD: 3-length phase slope, central difference (res 40)")
    ax[0].plot(f, FINAL["ng_12_24"], color=ORANGE, lw=0.8, alpha=0.8, label="FDTD: cascade pair N = 12, 24")
    for fn, mk in (("2d_pulse_f0.2350_res30.json", "s"), ("2d_pulse_f0.2564_res30.json", "s")):
        p = RES / fn
        if p.exists():
            d = json.loads(p.read_text()); ax[0].plot([d["fcen"]], [d["ng_obs"]], mk, color=AQUA, ms=6, label="time-domain pulse centroid" if "2350" in fn else None)
    shade(ax[0]); ax[0].set_ylim(3, 7); ax[0].set_ylabel(r"$n_g$"); ax[0].legend(frameon=False, fontsize=7)
    ax[1].plot(f, 100 * (FINAL["ng_fit"] / FINAL["ng_mpb"] - 1), color=BLUE, lw=1)
    ax[1].plot(f, 100 * (FINAL["ng_12_24"] / FINAL["ng_mpb"] - 1), color=ORANGE, lw=0.8, alpha=0.8)
    ax[1].axhline(0, color=MUTED, lw=0.8); shade(ax[1]); ax[1].set_ylim(-12, 12)
    ax[1].set_ylabel(r"$n_g^{FDTD}/n_g^{MPB}-1$ (%)"); ax[1].set_xlabel(r"$f$")
    fig.tight_layout(); fig.savefig(FIG / "2c_ng_fdtd_vs_mpb.png"); plt.close(fig)


def fig_length():
    Ns = np.array([6, 12, 24]); fig, ax = plt.subplots(1, 2, figsize=(8, 3.1))
    for f0, c in ((0.235, BLUE), (0.2564, ORANGE), (0.270, AQUA)):
        i = int(np.argmin(np.abs(f - f0)))
        ph = np.array([0.0] + [np.angle(FINAL[f"S21_N{N}"][i] / FINAL["S21_N6"][i]) for N in Ns[1:]])
        ref = FINAL["fit_slope"][i] * (Ns - 6)
        ph = ph + 2 * np.pi * np.round((ref - ph) / (2 * np.pi))
        ax[0].plot(Ns, ph, "o", color=c, label=f"f = {f[i]:.4f}")
        xs = np.linspace(6, 24, 20); ax[0].plot(xs, FINAL["fit_slope"][i] * (xs - 6), "-", color=c, lw=1)
        ax[1].plot(f0, FINAL["fit_rms"][i], "o", color=c)
    ax[0].set_xlabel("device length N (periods)"); ax[0].set_ylabel(r"arg $S_{21}(N)/S_{21}(6)$ (rad), branch-fixed"); ax[0].legend(frameon=False, fontsize=7)
    ax[1].semilogy(f, FINAL["fit_rms"], color=MUTED, lw=1); shade(ax[1]); ax[1].set_xlabel(r"$f$"); ax[1].set_ylabel("rms residual of linear-in-N fit (rad)")
    fig.tight_layout(); fig.savefig(FIG / "2d_length_scaling.png"); plt.close(fig)


def fig_convergence():
    C = json.loads((RES / "2e_convergence.json").read_text())["variants"]
    names = [k for k in ["base_res20", "base_res30", "base_res40", "dsrc6_res30", "dmon5_res30", "pml_res30", "tol7_res30", "tol8_res30", "tol7_res40", "J8_res30"] if k in C]
    reps = list(next(iter(C.values()))["rep"])
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.4))
    w = 0.2
    for j, r in enumerate(reps):
        ax[0].bar(np.arange(len(names)) + (j - 1.5) * w, [100 * C[n]["rep"][r]["rel_vs_mpb"] for n in names], w, label=r.split(" (")[0], color=[BLUE, ORANGE, AQUA, MAG][j])
    ax[0].axhline(0, color=INK, lw=0.6); ax[0].set_xticks(range(len(names))); ax[0].set_xticklabels(names, rotation=40, ha="right", fontsize=7)
    ax[0].set_ylabel(r"$n_g$ (N = 6, 12 pair) vs MPB, window mean (%)"); ax[0].legend(frameon=False, fontsize=6.5)
    ax[1].bar(range(len(names)), [100 * C[n]["vs_mpb_rms_useful"] for n in names], color=MUTED)
    ax[1].set_xticks(range(len(names))); ax[1].set_xticklabels(names, rotation=40, ha="right", fontsize=7); ax[1].set_ylabel("pointwise rms vs MPB, useful interval (%)")
    fig.tight_layout(); fig.savefig(FIG / "2e_convergence.png"); plt.close(fig)


def fig_sampling():
    S = json.loads((RES / "2f_sampling.json").read_text())["results"]["tol7_res40"]
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    for meth, c in (("central", BLUE), ("poly", ORANGE)):
        ks = [k for k in S if k.endswith(meth)]
        df = [float(k.split("=")[1].split("_")[0]) for k in ks]
        ax.plot(df, [100 * S[k]["rms_vs_mpb"] for k in ks], "o-", color=c, label=f"rms vs MPB, {meth}")
        ax.plot(df, [100 * S[k]["max_vs_mpb"] for k in ks], "s--", color=c, mfc="none", label=f"max vs MPB, {meth}")
    ax.set_xscale("log"); ax.set_xlabel(r"frequency spacing $\Delta f$"); ax.set_ylabel(r"$n_g$ error, useful interval (%)"); ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout(); fig.savefig(FIG / "2f_sampling.png"); plt.close(fig)


if __name__ == "__main__":
    for fn in (fig_mpb, fig_sparams, fig_ng, fig_length, fig_convergence, fig_sampling):
        fn(); print(fn.__name__, "ok")
