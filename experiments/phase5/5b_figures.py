"""Phase 5B: figures from results/phase5/5a_array.{npz,json}.  Everything is inside mask_common (0.231-0.268)."""
from _common import RES, FIG, np, THETA0, CASES, LABEL, COL, LS
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pcttd import beam as B

z = np.load(RES / "5a_array.npz"); j = json.load(open(RES / "5a_array.json"))
f = z["f"]; th = z["theta_grid"]; f0 = float(z["f0"]); d = float(z["d"]); pat = z["patterns_res30"]; pf = z["pat_freqs"]
g = lambda k, q, res=30: z[f"res{res}_{k}_{q}"]
S = j["setup"]; tol = S["tol_pt_deg"]
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130, "axes.spines.top": False, "axes.spines.right": False})
NOTE = "common band 0.231-0.268 (domain-limited)"


def deco(ax, ylab, title):
    ax.axvline(f0, color="k", lw=0.7, alpha=0.5); ax.set_xlabel("frequency f (c/a)"); ax.set_ylabel(ylab); ax.set_title(title, fontsize=9)


def save(fig, name):
    fig.tight_layout(); fig.savefig(FIG / name); plt.close(fig)


def lines(ax, q, cases=CASES, fn=lambda x: x):
    for k in cases:
        ax.plot(f, fn(g(k, q)), color=COL[k], ls=LS[k], label=LABEL[k], lw=1.6 if k != "ideal" else 1.2)


fig, ax = plt.subplots(figsize=(6.4, 4)); lines(ax, "theta_peak"); deco(ax, "design-lobe angle (deg)", f"Design-lobe pointing angle, {NOTE}")
ax.axhline(THETA0, color="gray", lw=0.6, ls=":"); ax.legend(fontsize=8); save(fig, "5a_theta_peak.png")

fig, ax = plt.subplots(figsize=(6.4, 4)); lines(ax, "theta_peak", fn=lambda x: x - THETA0); deco(ax, "pointing error (deg)", f"Pointing error theta_peak - 30 deg, {NOTE}")
ax.legend(fontsize=8); save(fig, "5b_pointing_error.png")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.8))
mx = [j["30"]["max_abs_err"][k] for k in CASES]; sp = [j["30"]["R_squint_spread"][k] * j["30"]["max_abs_err"]["fixed"] for k in CASES]
x = np.arange(4); a1.bar(x - 0.2, mx, 0.4, color="#4c78a8", label="max |theta - 30 deg| (absolute pointing)"); a1.bar(x + 0.2, sp, 0.4, color="#f58518", label="max |theta - theta(f0)| (squint spread)")
a1.set_xticks(x); a1.set_xticklabels(["fixed\nphase", "ideal\nTTD", "raw\nPhC", "trimmed\nPhC"]); a1.set_yscale("symlog", linthresh=0.01); a1.set_ylim(0, 30); a1.set_ylabel("deg (symlog)"); a1.legend(fontsize=7); a1.set_title("Absolute pointing vs squint spread", fontsize=9)
r1 = [j["30"]["R_squint"][k] for k in CASES]; r2 = [j["30"]["R_squint_spread"][k] for k in CASES]
a2.bar(x - 0.2, r1, 0.4, color="#4c78a8", label="R_squint (as defined)"); a2.bar(x + 0.2, r2, 0.4, color="#f58518", label="offset-removed spread ratio"); a2.axhline(0.2, color="gray", ls="--", lw=0.8, label="0.20 stretch target (not a criterion)")
a2.set_xticks(x); a2.set_xticklabels(["fixed\nphase", "ideal\nTTD", "raw\nPhC", "trimmed\nPhC"]); a2.set_yscale("symlog", linthresh=0.01); a2.set_ylim(0, 10); a2.set_ylabel("ratio to fixed phase (symlog)"); a2.legend(fontsize=7); a2.set_title("Squint-reduction ratios", fontsize=9)
save(fig, "5c_squint_comparison.png")

for key, name, zoom in (("full", "5d_patterns_full.png", False), ("zoom", "5e_design_lobe_zoom.png", True)):
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    for i, ax in enumerate(axs):
        for c, k in enumerate(CASES):
            ax.plot(th, 10 * np.log10(np.maximum(pat[c, i], 1e-6)), color=COL[k], ls=LS[k], label=LABEL[k], lw=1.3)
        ax.set_title(f"f = {pf[i]:.4f}" + (" (f0)" if i == 1 else " (band edge)"), fontsize=9); ax.set_xlabel("angle theta (deg)"); ax.set_ylim(-40, 2)
        if zoom:
            ax.set_xlim(10, 50)
        else:
            ax.set_xlim(-90, 90); ax.axvspan(S["window_theta_deg"][0], S["window_theta_deg"][1], color="#2ca02c", alpha=0.06)
    axs[0].set_ylabel("power |AF|^2/16 (dB)"); axs[0].legend(fontsize=7, loc="lower left")
    fig.suptitle("Full visible-region patterns (green: design-lobe tracking window)" if not zoom else "Design-lobe zoom", fontsize=9); save(fig, name)

fig, ax = plt.subplots(figsize=(6.4, 4)); lines(ax, "gain"); deco(ax, "peak power |AF|^2/16", f"Design-lobe peak power (1 = ideal equal-amplitude)"); ax.legend(fontsize=8); ax.set_ylim(0.99, 1.002); save(fig, "5f_peak_gain.png")
fig, ax = plt.subplots(figsize=(6.4, 4)); lines(ax, "hpbw"); deco(ax, "local HPBW of design lobe (deg)", f"Half-power beamwidth (-3 dB, local), {NOTE}"); ax.legend(fontsize=8); save(fig, "5g_hpbw.png")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.8))
lines(a1, "strongest_db"); deco(a1, "strongest other lobe rel. to design peak (dB)", "Strongest sidelobe/grating lobe"); a1.set_ylim(-3, 1); a1.legend(fontsize=7)
for k in CASES:
    a2.plot(f, g(k, "n_3db"), color=COL[k], ls=LS[k], drawstyle="steps-mid", label=LABEL[k])
deco(a2, "interior lobes within 3 dB of design peak", "Principal (grating) lobe count"); a2.set_ylim(0, 6)
save(fig, "5h_strongest_lobe.png")

fig, axs = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)
for ax, key, ttl in zip(axs, ("raw", "trim"), ("raw PhC", "phase-intercept-trimmed PhC")):
    e = g(key, "eps")
    for m in range(1, 4):
        ax.plot(f, e[m], label=f"channel {m} (N = {8 + m})")
    deco(ax, "phase error vs ideal TTD (rad)", f"{ttl}: arg(H_m/H_0) - 2 pi f m dtau"); ax.legend(fontsize=7)
save(fig, "5i_relative_phase_error.png")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.8))
for k in ("raw", "trim"):
    a1.plot(f, g(k, "theta_peak") - g(k + "_unit", "theta_peak"), color=COL[k], label=LABEL[k])
    a2.plot(f, 10 * np.log10(g(k, "gain") / g(k + "_unit", "gain")), color=COL[k], label=LABEL[k])
deco(a1, "theta(actual |H|) - theta(|H| = 1) (deg)", "Pointing: amplitude effect"); deco(a2, "peak gain ratio (dB)", "Peak gain: actual vs unit amplitudes"); a1.legend(fontsize=7)
a1.axhline(tol / 3, color="gray", ls=":", lw=0.8); save(fig, "5j_amplitude_sensitivity.png")

fig, ax = plt.subplots(figsize=(6.6, 4))
for q in range(-3, 4):
    s = np.sin(np.deg2rad(g("ideal", "theta_peak"))) + q / (f * d); ok = np.abs(s) <= 1
    ax.plot(f[ok], np.rad2deg(np.arcsin(s[ok])), color="k" if q == 0 else "#7f7f7f", lw=2 if q == 0 else 1, label="design lobe (q = 0)" if q == 0 else ("grating orders q != 0" if q == -3 else None))
sq = np.sin(np.deg2rad(g("fixed", "theta_peak"))) 
for q in (-1, 1):
    s = sq + q / (f * d); ok = np.abs(s) <= 1; ax.plot(f[ok], np.rad2deg(np.arcsin(s[ok])), color=COL["fixed"], ls=":", lw=1, label="fixed-phase orders" if q == -1 else None)
deco(ax, "lobe angle (deg)", "Visible principal lobes (analytic), d = 1.96 lambda0"); ax.set_ylim(-90, 90); ax.legend(fontsize=7); save(fig, "5k_grating_lobes.png")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.4, 3.9))
for k, c in (("trim", COL["trim"]),):
    e = g(k, "eps"); a1.plot(f, np.abs(e).max(0), color=c, label="max_m |eps_m(f)| (Phase 4 trimmed phase error)")
a1.axhline(0.026, color="gray", ls=":", lw=0.8, label="0.026 rad"); deco(a1, "rad", "Phase 4 residual phase error (trimmed)"); a1.legend(fontsize=7)
a2.plot(f, g("trim", "theta_peak") - g("ideal", "theta_peak"), color=COL["trim"], lw=1.8, label="numerical theta_trim - theta_ideal")
a2.plot(f, g("trim", "pred"), color="k", ls="--", lw=1, label="linear phase-slope prediction")
a2.axhspan(-tol, tol, color="#2ca02c", alpha=0.1, label=f"tolerance +-{tol:.2f} deg (0.05 HPBW)")
deco(a2, "deg", "Phase 5 pointing error from the same channels"); a2.legend(fontsize=7); save(fig, "5l_phase_to_pointing.png")
print(sorted(p.name for p in FIG.iterdir()))
