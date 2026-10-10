"""Phase 7B: publication figures from existing result files only (no FDTD). Output: figures/final/fig01..fig09 (.png and .pdf)."""
from _common import RES, FIG, jload, nload, json, np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch
from pcttd import geometry as geo

# Okabe-Ito colour-blind-safe palette, fixed assignment per entity across all figures
C = dict(fixed="#D55E00", ideal="#000000", ideal_sub="#7F7F7F", raw="#E69F00", trim="#0072B2", p4="#009E73", D0="#0072B2", D1="#CC79A7", D2="#E69F00", mpb="#7F7F7F")
CH = ["#0072B2", "#009E73", "#E69F00", "#CC79A7"]
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8, "lines.linewidth": 1.4,
                     "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans"})
LO, HI = 0.231, 0.268
F0 = 0.2495


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png"); fig.savefig(FIG / f"{name}.pdf"); plt.close(fig)


def tag(ax, s):
    ax.text(-0.2, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def band(ax):
    ax.axvspan(LO, HI, color="#0072B2", alpha=0.07, lw=0)


p4 = nload("phase4/phase4_handoff.npz"); p5 = nload("phase5/phase5_handoff.npz"); p6 = nload("phase6/phase6_handoff.npz")
p6j = jload("phase6/phase6_handoff.json"); p6r = jload("phase6/6b_robustness.json"); FR = jload("final/final_results.json")
f6 = p6["f"]; th6 = p6["theta_pat"]

# ---------------- Fig 1: concept -----------------------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 6.2)); gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.25])
ax = fig.add_subplot(gs[0]); ax.set_aspect("equal"); ax.grid(False); ax.spines[:].set_visible(False); ax.set_xticks([]); ax.set_yticks([])
H, sc = geo.triangular_w1(6, 0.3)
Np = 11
ax.add_patch(Rectangle((0, -H / 2), Np, H, fc="#dbe9f6", ec="none"))
for n in range(-1, Np + 1):
    for s in sc:
        x = s.x + n
        if -0.5 <= x <= Np + 0.5:
            ax.add_patch(Circle((x, s.y), s.radius, fc="white", ec="#444", lw=0.5))
ax.set_xlim(-0.8, Np + 0.8); ax.set_ylim(-H / 2 - 0.2, H / 2 + 0.2)
ax.annotate("", xy=(Np + 0.6, -0.45), xytext=(-0.7, -0.45), arrowprops=dict(arrowstyle="-|>", color="#D55E00", lw=1.5))
ax.text(Np / 2, 0.0, "guided W1 mode, n$_g$ ≈ 3.9", ha="center", va="center", fontsize=8, color="#D55E00", bbox=dict(fc="#dbe9f6", ec="none", pad=1.5))
ax.set_title("(a) W1 channel: triangular air holes, ε = 12, r/a = 0.30, TE.\n2-D FDTD domain (Phases 2–4), N = 8–11 periods; drawn for N = 11", fontsize=8.5)
ax = fig.add_subplot(gs[1]); ax.set_xlim(0, 16); ax.set_ylim(0, 6.6); ax.axis("off"); ax.grid(False)
ax.set_title("(b) Beamforming model (Phase 6): array factor from the simulated transfer functions $H_m(f)$ (analytic, not FDTD)", fontsize=8.5)
for q in range(4):
    x0 = 0.4 + q * 3.9
    ax.add_patch(Rectangle((x0, 4.1), 3.4, 1.0, fc="#dbe9f6", ec=CH[q], lw=1.4)); ax.text(x0 + 1.7, 4.6, f"PhC channel, N = {8 + q}", ha="center", va="center", fontsize=8)
    ax.annotate("", xy=(x0 + 1.7, 4.1), xytext=(x0 + 1.7, 3.15), arrowprops=dict(arrowstyle="-", color="#666"))
    ax.text(x0 + 1.7, 3.35, "trim φ$_q$", ha="center", fontsize=7, color="#555", bbox=dict(fc="white", ec="none", pad=1))
    for p in range(4):
        xe = x0 + 0.35 + p * 0.9
        ax.plot([xe, xe], [2.0, 2.6], color="k", lw=1); ax.plot([xe - 0.2, xe, xe + 0.2], [2.6, 3.0, 2.6], color="k", lw=1)
        ax.plot([x0 + 1.7, xe], [3.15, 2.6], color="#aaa", lw=0.6)
    ax.text(x0 + 1.7, 1.45, "4 elements, λ$_0$/2 pitch\nfixed phase steering", ha="center", fontsize=7)
    if q < 3: ax.annotate("", xy=(x0 + 3.9 + 1.7, 5.65), xytext=(x0 + 1.7, 5.65), arrowprops=dict(arrowstyle="<->", color="#888"))
ax.text(8, 6.05, "subarray pitch 2λ$_0$; PhC one-period step = 0.981 cycle; required step = 1.000 cycle (θ$_0$ = 30°)", ha="center", fontsize=8)
ax.text(8, 0.55, "Not modeled: loss, coupling, bends, disorder, thermal effects, 3-D", ha="center", fontsize=8, style="italic", color="#444")
save(fig, "fig01_concept")

# ---------------- Fig 2: verification ----------------------------------------------------------------------------------------------------
j1a = jload("phase1/1a_vacuum.json")["resolution"]["rows"]; j1c = jload("phase1/1c_bragg.json")["resolution"]
z2 = nload("phase2/2c_tol7_res40.npz"); z2b = nload("phase2/2c_tol7_res30.npz"); uj = jload("phase2/2c_tol7_res40.json")["useful_interval"]
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.9))
ax = axs[0]; r = np.array([x["value"] for x in j1a]); e = np.array([x["mean_err"] for x in j1a])
ax.loglog(r, e, "o-", color=C["trim"]); ax.loglog(r, e[0] * (r[0] / r) ** 2, "--", color="#999", label="slope −2")
ax.set_xlabel("resolution (px/a)"); ax.set_ylabel("vacuum delay error (a/c)"); ax.legend(); ax.set_title("1-D vacuum delay", fontsize=8.5); tag(ax, "a")
ax = axs[1]; r = np.array([x["value"] for x in j1c]); e = np.array([x["median_rel_complex_t"] for x in j1c])
ax.loglog(r, e, "o-", color=C["trim"]); ax.loglog(r, e[0] * (r[0] / r) ** 2, "--", color="#999", label="slope −2")
ax.set_xlabel("resolution (px/a)"); ax.set_ylabel("median |t$_{FDTD}$ − t$_{TMM}$|/|t|"); ax.legend(); ax.set_title("Bragg stack vs TMM", fontsize=8.5); tag(ax, "b")
ax = axs[2]; fz = z2["f"]; sel = (fz >= uj[0]) & (fz <= uj[1])
ax.plot(fz[sel], z2["ng_mpb"][sel], color=C["mpb"], lw=2.2, label="MPB")
ax.plot(fz[sel], z2["ng_fit"][sel], ".", ms=3, color=C["trim"], label="FDTD (res 40)")
ax.axvspan(LO, HI, color="#0072B2", alpha=0.07, lw=0); ax.set_xlabel("frequency f (c/a)"); ax.set_ylabel("group index n$_g$"); ax.legend(); ax.set_title("W1: FDTD vs MPB", fontsize=8.5); tag(ax, "c")
ax.set_ylim(3.7, 4.4)
save(fig, "fig02_verification")

# ---------------- Fig 3: delay-bandwidth trade ---------------------------------------------------------------------------------------------
mp = nload("phase3/3d_mpb_analysis.npz")
fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9))
ax = axs[0]
for d in ("D0", "D1", "D2"):
    ax.plot(mp["f"], mp[f"{d}_ng"], color=C[d], label=d)
    z3 = nload(f"phase3/3e_{d}_res30.npz"); m = (z3["f"] >= LO) & (z3["f"] <= HI)
    ax.plot(z3["f"][m][::3], z3["tau_fit"][m][::3], ".", ms=3, color=C[d])
ax.axvspan(LO, HI, color="#0072B2", alpha=0.07, lw=0); ax.set_xlim(0.22, 0.275); ax.set_ylim(3.5, 8); ax.set_xlabel("frequency f (c/a)"); ax.set_ylabel("n$_g$ (lines: MPB; dots: FDTD)"); ax.legend(); tag(ax, "a")
ax = axs[1]
bw = [100 * FR[f"B_{d}_res30_frac_bw"]["value"] for d in ("D0", "D1", "D2")]; ng = [FR[f"B_{d}_res30_median_tau"]["value"] for d in ("D0", "D1", "D2")]
for d, x, y in zip(("D0", "D1", "D2"), ng, bw):
    ax.plot(x, y, "o", ms=8, color=C[d]); ax.annotate(d, (x, y), textcoords="offset points", xytext=(6, 4), fontsize=8)
ax.set_xlabel("median delay per period (a/c)"); ax.set_ylabel("±10 % flat-delay bandwidth (%)"); ax.set_xlim(3.8, 4.6); ax.set_ylim(0, 17)
ax.set_title("D1, D2: +7 %, +14 % delay,\n−43 %, −62 % bandwidth", fontsize=8); tag(ax, "b")
ax = axs[2]; g = [abs(FR[f"B_{d}_res30_gdd_band_avg_linfit"]["value"]) for d in ("D0", "D1", "D2")]
ax.bar(["D0", "D1", "D2"], g, color=[C[d] for d in ("D0", "D1", "D2")], width=0.55); ax.set_ylabel("|GDD| per period (a/c)$^2$"); ax.set_title("D0 selected (plain W1)", fontsize=8); tag(ax, "c")
save(fig, "fig03_tradeoff")

# ---------------- Fig 4: four-channel family ----------------------------------------------------------------------------------------------
f4 = p4["f"]; m4 = p4["mask_common"].astype(bool); tgt = float(p4["dtau_target"])
fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.2))
ax = axs[0, 0]
for q in range(4): ax.plot(f4[m4], 20 * np.log10(np.abs(p4["H_raw"][q, m4])), color=CH[q], label=f"N = {8 + q}")
ax.set_ylim(-0.1, 0.1); ax.set_xlabel("f (c/a)"); ax.set_ylabel("|H$_m$| (dB)"); ax.legend(ncol=2); ax.set_title(f"amplitude spread ≤ {FR['D_amp_imbalance_max_db']['value']:.3f} dB", fontsize=8.5); tag(ax, "a")
ax = axs[0, 1]
ph = p4["phase_abs_unwrapped"]; dd = []
for q in range(4):
    c = np.polyfit(f4[m4] - F0, ph[q, m4], 3); dc = np.polyder(c); dd.append(np.polyval(dc, f4[m4] - F0) / (2 * np.pi))
dd = np.array(dd)
for q in range(3): ax.plot(f4[m4], dd[q + 1] - dd[q], color=CH[q + 1], label=f"N={9 + q} − N={8 + q}")
ax.axhline(tgt, color="k", ls="--", lw=0.8, label="target 3.9323"); ax.set_ylim(3.7, 4.2); ax.set_xlabel("f (c/a)"); ax.set_ylabel("delay increment (a/c)"); ax.legend(ncol=2); ax.set_title("increments (cubic fit of phase)", fontsize=8.5); tag(ax, "b")
ax = axs[1, 0]; i0 = int(np.argmin(np.abs(f4 - F0)))
def phase_err(H):
    rel = np.unwrap(np.angle(H / H[0]), axis=1); e = rel - 2 * np.pi * f4[None, :] * (np.arange(4)[:, None] * tgt)
    return e - 2 * np.pi * np.round(e[:, [i0]] / (2 * np.pi))
er = phase_err(p4["H_raw"]); et = phase_err(p4["H_trim"])
for q in range(1, 4): ax.plot(f4[m4], er[q, m4], color=CH[q])
ax.set_xlabel("f (c/a)"); ax.set_ylabel("phase error vs ideal TTD (rad)"); ax.set_title(f"raw: constant intercept up to {FR['E_raw_intercept_error_max_abs_rad']['value']:.2f} rad", fontsize=8.5); tag(ax, "c")
ax = axs[1, 1]
for q in range(1, 4): ax.plot(f4[m4], et[q, m4], color=CH[q], label=f"m = {q}")
ax.set_xlabel("f (c/a)"); ax.set_ylabel("trimmed phase error (rad)"); ax.legend(); ax.set_title(f"trimmed: ≤ {FR['E_trimmed_phase_residual_max_rad']['value']:.3f} rad", fontsize=8.5); tag(ax, "d")
for a in axs.ravel(): band(a)
save(fig, "fig04_channels")

# ---------------- Fig 5: Phase 5 limitation ------------------------------------------------------------------------------------------------
th5 = p5["theta_grid"]; pc = list(p5["pattern_cases"]); pf = p5["pattern_freqs"]
fig, axs = plt.subplots(1, 3, figsize=(7.8, 2.9))
ax = axs[0]
for k, lab in (("ideal", "ideal TTD"), ("trim", "trimmed PhC")):
    P = p5["patterns_res30"][pc.index(k), 1]; ax.plot(th5, 10 * np.log10(np.maximum(P, 1e-6)), color=C[k if k != "ideal" else "ideal"], lw=1.2 if k == "trim" else 0.9, ls="-" if k == "trim" else "--", label=lab)
ax.set_ylim(-30, 2); ax.set_xlim(-90, 90); ax.set_xlabel("θ (deg)"); ax.set_ylabel("power (dB re peak)"); ax.set_title("4 elements, d = 1.96 λ$_0$, f = f$_0$", fontsize=8); ax.legend(loc="lower left"); tag(ax, "a")
ax = axs[1]; f5 = p5["f"]
for k, lab in (("fixed", "fixed phase"), ("ideal", "ideal TTD"), ("raw", "raw PhC"), ("trim", "trimmed PhC")):
    ax.plot(f5, p5[f"res30_{k}_theta_peak"], color=C[k], label=lab, ls="--" if k == "ideal" else "-")
ax.set_xlabel("f (c/a)"); ax.set_ylabel("design-lobe angle (deg)"); ax.legend(); ax.set_title("squint: fixed ±2.7°", fontsize=8); tag(ax, "b")
ax = axs[2]; x = np.arange(4); labs = ["fixed", "ideal", "raw", "trim"]
ax.plot(x, [p5[f"res30_{k}_strongest_db"].max() for k in labs], "o", ms=8, color="k")
ax.set_xticks(x); ax.set_xticklabels(["fixed", "ideal", "raw", "trimmed"], fontsize=7); ax.set_xlim(-0.5, 3.5); ax.set_ylim(-14, 2); ax.axhline(0, color="k", lw=0.6)
ax.set_ylabel("strongest other lobe (dB re peak)"); ax.set_title("all cases ≈ 0 dB: grating lobes", fontsize=8); tag(ax, "c")
save(fig, "fig05_phase5_limit")

# ---------------- Fig 6: Phase 6 architecture --------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(7.8, 2.9))
ax = axs[0]; Ms = [1, 2, 4, 8]; x = np.arange(4)
req = [FR[f"H_msub{M}_required_cycles"]["value"] for M in Ms]; phc = FR["H_phc_delay_quantum_cycles"]["value"]
ax.bar(x - 0.19, req, 0.36, color="#7F7F7F", label="required step"); ax.bar(x + 0.19, [phc] * 4, 0.36, color=C["trim"], label="PhC, one period")
ax.set_xticks(x); ax.set_xticklabels([f"Msub={M}" for M in Ms], fontsize=7); ax.set_ylabel("delay step (carrier cycles at f$_0$)"); ax.legend(loc="upper left"); ax.set_title("Msub = 4: −1.89 % mismatch", fontsize=8); tag(ax, "a")
ax = axs[1]; qs = p6r["quantization_scan"]
ax.semilogy([q["msub"] for q in qs], [q["max_err"] for q in qs], "o-", color=C["trim"], ms=4)
ax.axhline(0.367, color="k", ls="--", lw=0.8); ax.text(12, 0.4, "tracking tolerance 0.37°", ha="right", fontsize=7)
ax.set_xlabel("elements per PhC channel, Msub"); ax.set_ylabel("trimmed max pointing error (deg)"); ax.set_title("one-period steps only", fontsize=8); tag(ax, "b")
ax = axs[2]; nl = [FR[f"K_msub{M}_trim_strongest_other_lobe_db"]["value"] for M in Ms]
ax.bar(x, nl, color=C["trim"], width=0.55); ax.set_xticks(x); ax.set_xticklabels([f"Msub={M}" for M in Ms], fontsize=7); ax.set_ylabel("strongest other lobe (dB)"); ax.set_ylim(-14, 0)
ax.set_title("Phase 5 element-level: 0 dB", fontsize=8); tag(ax, "c")
save(fig, "fig06_architecture")

# ---------------- Fig 7: Phase 6 beam patterns -------------------------------------------------------------------------------------------
pcs = list(p6["pattern_cases"]); pf6 = p6["pattern_freqs"]
fig, axs = plt.subplots(1, 3, figsize=(7.8, 2.8), sharey=True)
for i, ax in enumerate(axs):
    for k, lab, ls, lw in (("fixed", "fixed phase", "-", 1.0), ("ideal_sub", "ideal subarray TTD", "--", 1.2), ("trim", "trimmed PhC", "-", 1.3)):
        P = p6["patterns_res30_M4_" + k][i]; ax.plot(th6, 10 * np.log10(np.maximum(P, 1e-6)), color=C[k if k != "ideal_sub" else "ideal_sub"], ls=ls, lw=lw, label=lab)
    ax.axvline(30, color="k", lw=0.5, alpha=0.5); ax.set_ylim(-35, 2); ax.set_xlim(-90, 90); ax.set_xlabel("θ (deg)"); ax.set_title(f"f = {pf6[i]:.4g}", fontsize=8.5); tag(ax, "abc"[i])
axs[0].set_ylabel("power (dB re ideal peak)"); axs[1].legend(loc="lower left", fontsize=7)
save(fig, "fig07_patterns")

# ---------------- Fig 8: pointing vs frequency -------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw=dict(width_ratios=[1.4, 1]))
ax = axs[0]
for k, lab, ls in (("fixed", "fixed phase", "-"), ("ideal_sub", "ideal subarray TTD", "--"), ("trim", "trimmed PhC", "-"), ("trim_p4", "Phase 4 trim unchanged", ":"), ("raw", "raw PhC", "-")):
    c = {"trim_p4": C["p4"]}.get(k, C.get(k))
    ax.plot(f6, p6[f"res30_M4_{k}_theta_peak"] - 30, color=c, ls=ls, label=lab, lw=1.2)
ax.set_ylim(-7.5, 3.5); ax.set_xlabel("f (c/a)"); ax.set_ylabel("pointing error (deg)"); ax.legend(loc="lower right"); band(ax); tag(ax, "a"); ax.set_title("Msub = 4, 16 elements", fontsize=8.5)
ax = axs[1]
for k, lab, ls in (("ideal_sub", "ideal subarray TTD", "--"), ("trim", "trimmed PhC", "-")):
    ax.plot(f6, p6[f"res30_M4_{k}_theta_peak"] - 30, color=C[k], ls=ls, label=lab)
ax.axhline(0.0, color="k", lw=0.5); ax.axhspan(-0.367, 0.367, color="#009E73", alpha=0.1, lw=0); ax.text(0.2315, 0.31, "tolerance ±0.37°", fontsize=7, color="#006d4f")
ax.set_xlabel("f (c/a)"); ax.set_ylabel("pointing error (deg)"); ax.set_ylim(-0.5, 0.5); ax.set_title("zoom", fontsize=8.5); ax.legend(loc="lower left"); tag(ax, "b")
save(fig, "fig08_pointing")

# ---------------- Fig 9: robustness -------------------------------------------------------------------------------------------------------
mc = p6r["monte_carlo"]
fig, axs = plt.subplots(1, 4, figsize=(8.4, 2.8), sharey=True)
spec = [("phase", "eps_deg", "phase error ± (deg)"), ("delay_uncal", "delay_frac", "delay error ± (frac of q τ)"), ("delay_cal", "delay_frac", "delay error ± (calibrated)"), ("amp", "amp_db", "amplitude error ± (dB)")]
for ax, (k, sk, xl) in zip(axs, spec):
    xs = np.array([e["spec"][sk] for e in mc[k]]); p95 = np.array([e["max_err"]["p95"] for e in mc[k]]); md = np.array([e["max_err"]["median"] for e in mc[k]])
    o = np.argsort(xs); ax.plot(xs[o], p95[o], "o-", ms=3, color=C["trim"], label="P95"); ax.plot(xs[o], md[o], "s--", ms=3, color="#7F7F7F", label="median")
    ax.axhline(0.5, color="k", ls=":", lw=0.8); ax.set_xlabel(xl, fontsize=8)
    if k.startswith("delay"): ax.set_xlim(0, xs.max())
axs[0].set_ylabel("max pointing error (deg)"); axs[0].legend(loc="upper left"); axs[3].text(0.97, 0.52, "0.5° limit", transform=axs[3].get_yaxis_transform(), ha="right", fontsize=7)
for i, ax in enumerate(axs): tag(ax, "abcd"[i])
axs[0].set_ylim(0, 1.8)
save(fig, "fig09_robustness")
print("figures:", sorted(p.name for p in FIG.glob("*.png")))
