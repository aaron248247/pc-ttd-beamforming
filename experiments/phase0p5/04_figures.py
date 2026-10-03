"""Phase 0.5 figures. Reads ONLY saved results (results/phase0p5/*); recomputes nothing from MPB."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _common import *  # noqa
from pcttd import io as pio, bands as B

FIG = ROOT / "figures" / "phase0p5"
FIG.mkdir(parents=True, exist_ok=True)
COL = {"H": "#2a78d6", "R": "#eb6834", "H_r030": "#1baf7a"}     # fixed categorical order
LABEL = {"H": "H: air holes, r=0.40 (rule)", "R": "R: rods, r=0.19 (rule)", "H_r030": "H$_{r030}$: air holes, r=0.30 (suppl.)"}
INK, MUTED = "#0b0b0b", "#898781"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6, "figure.dpi": 130,
                     "axes.axisbelow": True})


def load_design(cfg, f1, s):
    z = np.load(RES / "defect" / cfg / f"f1_{f1:.2f}_s_{s:+.2f}.npz")
    return {k: z[k] for k in z.files}


def crossings(d, nlev=1500):
    gap = d["gap"]
    br = B.class_branches(d["freq"], d["vg"], d["parity"], d["strip"], d["edge"])
    lv = np.linspace(gap[0], gap[1], nlev + 2)[1:-1]
    return lv, B.find_crossings(d["k"], br, lv)


def fig_bulk_gap():
    bulk = pio.load_json(RES / "bulk_gap_scan.json")
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.2), sharey=True)
    for ax, lat, name, title in zip(axs, ("triangular", "square"), ("H", "R"),
                                    ("Triangular air holes, TE (Hz)", "Square rods, TM (Ez)")):
        rows = bulk[lat]["rows"]
        for n, ls in ((1, "-"), (2, "--"), (3, ":")):
            xs = [r["r"] for r in rows if any(g["n"] == n and g["rel"] > 1e-3 for g in r["gaps"])]
            ys = [g["rel"] for r in rows for g in r["gaps"] if g["n"] == n and g["rel"] > 1e-3]
            if xs:
                ax.plot(xs, ys, ls, color=COL[name], lw=1.6, marker="o", ms=3, label=f"gap {n}|{n+1}")
        ax.axvline(bulk[lat]["r_base"], color=MUTED, lw=1)
        ax.annotate(f"r$_{{base}}$={bulk[lat]['r_base']:.3g}\nΔω/ω={bulk[lat]['rel_gap_base']:.3f}",
                    (bulk[lat]["r_base"], 0.05), xytext=(5, 2), textcoords="offset points", color=INK, fontsize=8)
        ax.axhline(0.05, color=MUTED, lw=0.8, ls=":")
        ax.set_title(title, fontsize=9, color=INK)
        ax.set_xlabel("r / a"); ax.legend(frameon=False, fontsize=8, loc="upper left")
    axs[0].set_ylabel("relative gap width Δω/ω$_{mid}$")
    fig.tight_layout(); fig.savefig(FIG / "f1_bulk_gap_vs_radius.png"); plt.close(fig)


def fig_bulk_bands():
    z = np.load(RES / "bulk_bands.npz")
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.2), sharey=True)
    sel = pio.load_json(RES / "selection.json") if (RES / "selection.json").exists() else None
    for ax, n in zip(axs, ("H", "R", "H_r030")):
        p = z[f"{n}_path"]
        for b in range(p.shape[1]):
            ax.plot(p[:, b], color=COL[n], lw=1.3)
        both = np.vstack([p, z[f"{n}_grid"]])
        g = B.find_gaps(both, 3)[0]
        ax.axhspan(g.w_lo, g.w_hi, color=COL[n], alpha=0.15, lw=0)
        ax.set_title(f"{LABEL[n]}\ngap {g.w_lo:.3f}–{g.w_hi:.3f}  (Δω/ω={g.rel:.3f})", fontsize=8)
        ax.set_xlabel("k-path index (Γ→…→Γ)")
    axs[0].set_ylabel("ω a / 2πc")
    fig.tight_layout(); fig.savefig(FIG / "f2_bulk_bands.png"); plt.close(fig)


def fig_plain_w1():
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.6), sharey=False)
    for ax, n in zip(axs, ("H", "R", "H_r030")):
        d = load_design(n, 1.0, 0.0)
        g = d["gap"]
        for b in range(d["freq"].shape[1]):
            for p, ls in ((+1, "-"), (-1, "--")):
                w = np.where(d["parity"][:, b] == p, d["freq"][:, b], np.nan)
                inside = (w > g[0]) & (w < g[1])
                ax.plot(d["k"], np.where(inside, w, np.nan), ls, color=COL[n], lw=1.8)
                ax.plot(d["k"], np.where(~inside, w, np.nan), ls, color=MUTED, lw=0.6, alpha=0.7)
        lv, c = crossings(d)
        sa, se, eo = B.single_mode_masks(c)
        for a, b_ in B.runs(se):
            ax.axhspan(lv[a], lv[b_ - 1], color=COL[n], alpha=0.18, lw=0)
        ax.set_ylim(g[0] - 0.02, g[1] + 0.02)
        ax.axhline(g[0], color=MUTED, lw=0.8, ls=":"); ax.axhline(g[1], color=MUTED, lw=0.8, ls=":")
        ax.set_title(LABEL[n] + " — plain W1", fontsize=8)
        ax.set_xlabel("k a / 2π"); 
    axs[0].set_ylabel("ω a / 2πc")
    axs[0].plot([], [], "-", color=INK, label="even Φ"); axs[0].plot([], [], "--", color=INK, label="odd Φ")
    axs[0].fill_between([], [], color=COL["H"], alpha=0.18, label="single-mode (even) window")
    axs[0].legend(frameon=False, fontsize=7, loc="lower left")
    fig.tight_layout(); fig.savefig(FIG / "f3_plain_w1_dispersion.png"); plt.close(fig)


def fig_heatmaps():
    sw = pio.load_json(RES / "defect_sweep.json")
    keys = [("sm_strict_even_frac", "single-mode (even) fractional bandwidth", 1.0),
            ("ng_max_strict_even", "max n$_g$ in single-mode window", 1.0),
            ("gbp10", "GBP (n$_g$≥10, ±10% flat)", 1.0)]
    fig, axs = plt.subplots(3, 3, figsize=(9.5, 8))
    for j, n in enumerate(("H", "R", "H_r030")):
        ds = sw[n]["designs"]
        for i, (key, title, _) in enumerate(keys):
            Z = np.full((len(S), len(F1)), np.nan)
            for d in ds:
                if "skipped" in d: continue
                Z[S.index(round(d["s"], 2)), F1.index(round(d["f1"], 2))] = d.get(key, 0.0)
            ax = axs[i, j]
            im = ax.imshow(Z, origin="lower", cmap="Blues", aspect="auto")
            ax.grid(False)
            for a in range(len(S)):
                for b in range(len(F1)):
                    t = "×" if np.isnan(Z[a, b]) else f"{Z[a, b]:.3g}"
                    ax.text(b, a, t, ha="center", va="center", fontsize=7, color=INK)
            ax.set_xticks(range(len(F1)), [f"{x:g}" for x in F1]); ax.set_yticks(range(len(S)), [f"{x:+.2f}" for x in S])
            if i == 0: ax.set_title(LABEL[n], fontsize=8)
            if i == 2: ax.set_xlabel("first-row radius scale f1")
            if j == 0: ax.set_ylabel("outward shift s (a)\n" + title, fontsize=8)
    fig.text(0.5, 0.005, "× = unresolved at screening resolution (excluded)", ha="center", fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(FIG / "f4_sweep_heatmaps.png"); plt.close(fig)


def fig_ng():
    sw = pio.load_json(RES / "defect_sweep.json")
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.4))
    for ax, n in zip(axs, ("H", "R", "H_r030")):
        ds = [d for d in sw[n]["designs"] if "skipped" not in d]
        best = max(ds, key=lambda d: d["gbp10"]) if ds else None
        picks = [(1.0, 0.0, "plain W1", "-")]
        if best and (best["f1"], best["s"]) != (1.0, 0.0):
            picks.append((best["f1"], best["s"], f"best GBP: f1={best['f1']:g}, s={best['s']:+.2f}", "--"))
        for f1, s, lab, ls in picks:
            d = load_design(n, f1, s); lv, c = crossings(d)
            sa, se, eo = B.single_mode_masks(c)
            idx = np.nonzero(se[c.level_idx])[0]
            ax.plot(lv[c.level_idx[idx]], 1 / np.abs(c.vg[idx]), ls, color=COL[n], lw=1.6, label=lab)
        ax.set_yscale("log"); ax.set_xlabel("ω a / 2πc"); ax.set_title(LABEL[n], fontsize=8)
        ax.legend(frameon=False, fontsize=7)
    axs[0].set_ylabel("group index n$_g$ (single-mode even window only)")
    fig.tight_layout(); fig.savefig(FIG / "f5_group_index.png"); plt.close(fig)


def fig_criteria():
    sel = pio.load_json(RES / "selection.json")["metrics"]
    keys = [("rel_gap", "C1 gap"), ("c2_max_sm_frac", "C2 SM bandwidth"), ("c3_median_strip", "C3 confinement"),
            ("c4_max_gbp10", "C4 GBP(10)"), ("c5_pixels_per_period", "C5 px/period (lower=better)"),
            ("c6_fraction", "C6 FDTD-suitable fraction")]
    fig, axs = plt.subplots(1, 6, figsize=(11, 2.8))
    for ax, (k, t) in zip(axs, keys):
        vals = [sel[n][k] or 0 for n in ("H", "R", "H_r030")]
        ax.bar(range(3), vals, color=[COL[n] for n in ("H", "R", "H_r030")], width=0.6)
        ax.set_xticks(range(3), ["H", "R", "H$_{r030}$"]); ax.set_title(t, fontsize=8)
        for i, v in enumerate(vals): ax.text(i, v, f"{v:.3g}", ha="center", va="bottom", fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "f6_criteria_summary.png"); plt.close(fig)


if __name__ == "__main__":
    fig_bulk_gap(); fig_bulk_bands(); fig_plain_w1(); fig_heatmaps(); fig_ng(); fig_criteria()
    print("figures written to", FIG)
