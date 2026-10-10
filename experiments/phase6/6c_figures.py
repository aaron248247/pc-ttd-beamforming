"""Phase 6C: figures from results/phase6/6a_subarrays, 6b_robustness (and read-only Phase 5 handoff)."""
from _common import RES, P5, FIG, np, THETA0, MSUBS, CASES, LABEL, SHORT, COL, LS, MCOL
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

z = np.load(RES / "6a_subarrays.npz"); j = json.load(open(RES / "6a_subarrays.json")); zb = np.load(RES / "6b_robustness.npz"); jb = json.load(open(RES / "6b_robustness.json"))
p5 = np.load(P5 / "phase5_handoff.npz"); p5j = json.load(open(P5 / "phase5_handoff.json"))
f = z["f"]; f0 = float(z["f0"]); I = z["i_lpf"]; th_pat = z["theta_pat"]; tau_phc = float(z["tau_phc"]); lam0 = float(z["lam0"]); d_e = float(z["d_e"])
best = jb["best"]; BM = best["msub"]; BC = best["case"]
g = lambda res, M, k, q: z[f"res{res}_M{M}_{k}_{q}"]
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130, "axes.spines.top": False, "axes.spines.right": False})
NOTE = "common band 0.231-0.268 (domain-limited)"


def deco(ax, ylab, title, x=True):
    ax.axvline(f0, color="k", lw=0.7, alpha=0.5); ax.set_xlabel("frequency f (c/a)"); ax.set_ylabel(ylab); ax.set_title(title, fontsize=9)


def save(fig, name):
    fig.tight_layout(); fig.savefig(FIG / name); plt.close(fig)


# 1 required delay vs subarray size ------------------------------------------------------------------------------------------------------
Mc = np.linspace(0.5, 12, 200)
fig, ax = plt.subplots(figsize=(6.4, 4))
ax.plot(Mc, f0 * Mc * lam0 / 4, color="k", label="required step  f0 Msub lambda0/4 (cycles)")
ax.axhline(f0 * tau_phc, color="#1f77b4", label=f"PhC one-period step  f0 x 3.9323 = {f0 * tau_phc:.4f} cycle")
for M in MSUBS:
    ax.plot(M, M / 4, "o", color=MCOL[M]); ax.annotate(f"Msub={M}\n{M / 4:.2f} cycle", (M, M / 4), textcoords="offset points", xytext=(4, 6), fontsize=7)
ax.plot(4, f0 * tau_phc, "s", color="#1f77b4"); ax.set_xlabel("elements per PhC channel, Msub"); ax.set_ylabel("delay step between adjacent channels (carrier cycles at f0)")
ax.set_title("Required inter-subarray delay vs subarray size (theta0 = 30 deg, lambda0/2 elements)", fontsize=9); ax.legend(fontsize=8); save(fig, "6_01_required_delay.png")

# 2 delay quantum vs required --------------------------------------------------------------------------------------------------------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
x = np.arange(len(MSUBS)); req = [j["arch"][str(M)]["req_cycles_f0"] for M in MSUBS]
a1.bar(x - 0.2, req, 0.4, color="#7f7f7f", label="required"); a1.bar(x + 0.2, [f0 * tau_phc] * 4, 0.4, color="#1f77b4", label="PhC (N -> N+1 period)")
for i, M in enumerate(MSUBS):
    a1.text(i, max(req[i], f0 * tau_phc) + 0.05, f"rho={j['arch'][str(M)]['rho']:.3f}", ha="center", fontsize=8)
a1.set_xticks(x); a1.set_xticklabels([f"Msub={M}" for M in MSUBS]); a1.set_ylabel("cycles at f0"); a1.legend(fontsize=8); a1.set_title("PhC delay quantum vs required step", fontsize=9)
for k, mk in ((1, "o"), (2, "s"), (3, "^")):
    xs = np.arange(1, 13); a2.plot(xs, [j["quantization"][f"M{m}_k{k}"]["rho"] for m in xs], mk + "-", ms=4, label=f"k = {k} period(s) per step")
a2.axhspan(0.95, 1.05, color="#2ca02c", alpha=0.15, label="|rho - 1| <= 5 %"); a2.axhline(1, color="k", lw=0.6); a2.set_yscale("log"); a2.set_xlabel("Msub"); a2.set_ylabel("rho = k tau_PhC / tau_req"); a2.legend(fontsize=7)
a2.set_title("Integer-period quantization of the delay ratio", fontsize=9); save(fig, "6_02_delay_quantum.png")

# 3 theta_peak & 4 pointing error --------------------------------------------------------------------------------------------------------
for name, fn, yl, ttl in (("6_03_theta_peak.png", lambda x: x, "design-lobe angle (deg)", "Design-lobe angle"), ("6_04_pointing_error.png", lambda x: x - THETA0, "pointing error (deg)", "Pointing error")):
    fig, axs = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for ax, M in zip(axs.ravel(), MSUBS):
        for k in ("fixed", "ideal_sub", "trim", "trim_p4", "raw"):
            y = fn(g(30, M, k, "theta_peak"))
            if k == "raw" and M == 1: y = np.where(np.abs(y - THETA0) > 15, np.nan, y)
            ax.plot(f, y, color=COL[k], ls=LS[k], label=SHORT[k], lw=1.5)
        ax.set_title(f"Msub = {M} ({4 * M} elements)  rho = {j['arch'][str(M)]['rho']:.3f}", fontsize=9); ax.axvline(f0, color="k", lw=0.6, alpha=0.5); ax.set_ylabel(yl)
        if name.startswith("6_04"): ax.set_ylim(-8, 8)
    axs[0, 0].legend(fontsize=7); [a.set_xlabel("f (c/a)") for a in axs[1]]; fig.suptitle(f"{ttl}, {NOTE}  (raw clipped where off-scale)", fontsize=9); save(fig, name)

# 5 full patterns best -------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
for ax, i, lab in zip(axs, range(3), ("lower edge 0.231", "f0 = 0.2495", "upper edge 0.268")):
    for k in ("fixed", "ideal_elem", BC):
        p = z[f"pat_M{BM}_{k}"][i].astype(float)
        ax.plot(th_pat, 10 * np.log10(np.maximum(p, 1e-6)), color=COL[k], ls=LS[k], label=SHORT[k], lw=1.3)
    for a_ in j["subarray_level_replica_angles_f0"][str(BM)]: pass
    ax.set_ylim(-45, 2); ax.set_xlim(-90, 90); ax.set_xlabel("angle (deg)"); ax.set_title(f"Msub={BM}, {lab}", fontsize=9)
axs[0].set_ylabel("normalized power (dB re N_el^2)"); axs[0].legend(fontsize=7); save(fig, "6_05_full_patterns_best.png")

# 6 Phase 5 vs best ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9))
a = axs[0]
a.plot(f, p5["res30_trim_pointing_error"], color="#9467bd", label="Phase 5: element-level trimmed PhC (d = 1.96 lambda0)"); a.plot(f, p5["res30_fixed_pointing_error"], color="#9467bd", ls=":", label="Phase 5 fixed phase")
a.plot(f, g(30, BM, BC, "theta_peak") - THETA0, color=COL["trim"], label=f"Phase 6: Msub={BM} {SHORT[BC]} (lambda0/2 elements)"); a.plot(f, g(30, BM, "fixed", "theta_peak") - THETA0, color=COL["fixed"], ls=":", label=f"Phase 6 fixed phase (Msub={BM} array)")
deco(a, "pointing error (deg)", "Pointing error", False); a.legend(fontsize=6.5)
pc5 = list(p5["pattern_cases"]); k5 = pc5.index("trim"); th5 = p5["theta_grid"]
for ax, i, lab in ((axs[1], 1, "f0"), (axs[2], 2, "upper edge 0.268")):
    ax.plot(th5, 10 * np.log10(np.maximum(p5["patterns_res30"][k5][i], 1e-6)), color="#9467bd", label="Phase 5 trimmed PhC (4 el., 1.96 lambda0)", lw=1.2)
    ax.plot(th_pat, 10 * np.log10(np.maximum(z[f"pat_M{BM}_{BC}"][i].astype(float), 1e-6)), color=COL["trim"], label=f"Phase 6 Msub={BM} {SHORT[BC]} (16 el.)", lw=1.2)
    ax.set_ylim(-40, 2); ax.set_xlim(-90, 90); ax.set_xlabel("angle (deg)"); ax.set_ylabel("normalized power (dB)"); ax.set_title(f"Full pattern at {lab}", fontsize=9); ax.legend(fontsize=6.5)
save(fig, "6_06_phase5_vs_phase6.png")

# 7 strongest lobe, 8 gain, 9 HPBW -------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(10, 3.9))
for M in MSUBS:
    axs[0].plot(f, g(30, M, "trim", "strongest_db"), color=MCOL[M], label=f"Msub={M}"); axs[1].plot(f, g(30, M, "trim", "grating_db"), color=MCOL[M], label=f"Msub={M}")
axs[0].plot(f, p5["res30_trim_strongest_db"], color="#9467bd", ls="--", label="Phase 5 (d = 1.96 lambda0)")
axs[0].axhline(-10, color="gray", ls=":", lw=0.8); axs[1].axhline(-10, color="gray", ls=":", lw=0.8)
deco(axs[0], "strongest other lobe (dB re design peak)", "Strongest non-design lobe (trimmed PhC)"); deco(axs[1], "strongest subarray replica (dB re design peak)", "Strongest grating (replica) lobe, trimmed PhC")
axs[0].legend(fontsize=7); axs[1].legend(fontsize=7); save(fig, "6_07_strongest_lobe.png")
fig, ax = plt.subplots(figsize=(6.4, 4))
for M in MSUBS: ax.plot(f, 10 * np.log10(g(30, M, "trim", "gain")), color=MCOL[M], label=f"Msub={M} trimmed"); ax.plot(f, 10 * np.log10(g(30, M, "raw", "gain")), color=MCOL[M], ls=":", lw=1)
ax.plot(f, 10 * np.log10(p5["res30_trim_gain"]), color="#9467bd", ls="--", label="Phase 5 trimmed"); deco(ax, "peak power (dB re N_el^2)", "Peak power of the design lobe (dotted: raw)"); ax.legend(fontsize=7); save(fig, "6_08_peak_gain.png")
fig, ax = plt.subplots(figsize=(6.4, 4))
for M in MSUBS: ax.plot(f, g(30, M, "trim", "hpbw"), color=MCOL[M], label=f"Msub={M} trimmed"); ax.plot(f, g(30, M, "ideal_elem", "hpbw"), color=MCOL[M], ls=":", lw=1)
ax.plot(f, p5["res30_trim_hpbw"], color="#9467bd", ls="--", label="Phase 5 trimmed"); deco(ax, "HPBW (deg)", "Half-power beamwidth (dotted: ideal element-level TTD)"); ax.set_yscale("log"); ax.legend(fontsize=7); save(fig, "6_09_hpbw.png")

# 10 raw vs trimmed ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(10, 3.9))
for k in ("fixed", "ideal_sub", "raw", "trim", "trim_p4"):
    axs[0].plot(f, g(30, 4, k, "theta_peak") - THETA0, color=COL[k], ls=LS[k], label=LABEL[k] if k in ("trim", "trim_p4") else SHORT[k])
    axs[1].plot(f, g(30, 4, k, "theta_peak") - g(30, 4, k, "theta_peak")[int(I[1])], color=COL[k], ls=LS[k])
deco(axs[0], "pointing error (deg)", "Msub = 4: absolute pointing error"); deco(axs[1], "theta - theta(f0) (deg)", "Msub = 4: squint spread (offset removed)"); axs[0].legend(fontsize=6.5); save(fig, "6_10_raw_vs_trim.png")

# 11 delay mismatch -> pointing ----------------------------------------------------------------------------------------------------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.9))
a1.plot(z["rho_grid"], z["mm4_max_err_retrim"], color=COL["trim"], label="max |error| over band, retrimmed at f0"); a1.plot(z["rho_grid"], np.abs(z["mm4_err_f0_noretrim"]), color="#ff7f0e", label="|error| at f0, NOT re-trimmed (intercept off)")
a1.axvline(tau_phc / (4 * lam0 / 4), color="k", ls=":", label=f"PhC rho = {tau_phc / lam0:.4f}"); a1.set_yscale("symlog", linthresh=0.05); a1.set_xlabel("rho = actual / required delay step (Msub = 4, ideal-delay channels)"); a1.set_ylabel("deg")
a1.set_title("Delay-mismatch -> pointing error", fontsize=9); a1.legend(fontsize=7)
q = jb["quantization_scan"]; xs = [r["msub"] for r in q]
a2.plot(xs, [r["max_err"] for r in q], "o-", color=COL["trim"], label="max |error| (retrimmed PhC, k = 1)"); a2.plot(xs, [r["R_squint"] for r in q], "s--", color="#d62728", label="R_squint"); a2.set_yscale("log"); a2.set_xlabel("Msub (one-period channels)"); a2.legend(fontsize=7)
a2.set_title("Consequence of integer-period delay: all Msub, k = 1", fontsize=9); save(fig, "6_11_delay_mismatch.png")

# 12-14 sensitivity ----------------------------------------------------------------------------------------------------------------------
def sens(kinds, det_key, xl, name, ttl, xfac=1.0):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.9))
    for kind, col, lab in kinds:
        L = [m["level"] for m in jb["monte_carlo"][kind]]
        for st, ls in (("median", "-"), ("p95", "--"), ("worst", ":")):
            axs[0].plot(L, [m["max_err"][st] for m in jb["monte_carlo"][kind]], color=col, ls=ls, label=f"{lab} {st}" if True else None)
        axs[1].plot(L, [m["strongest_db"]["p95"] for m in jb["monte_carlo"][kind]], color=col, label=lab + " P95");
    if det_key:
        d = jb["deterministic"][det_key]; axs[0].plot([r["level"] for r in d], [r["shift_max"] for r in d], "k.-", lw=0.8, label="deterministic ramp (numeric shift)")
        axs[0].plot([r["level"] for r in d], [r["pred_max"] for r in d], "k--", lw=0.8, alpha=0.6, label="first-order analytic")
    for lim in (0.5, 1.0): axs[0].axhline(lim, color="gray", ls="-.", lw=0.6)
    axs[0].set_xscale("log"); axs[0].set_yscale("log"); axs[0].set_xlabel(xl); axs[0].set_ylabel("max |pointing error| over band (deg)"); axs[0].set_title(ttl, fontsize=9); axs[0].legend(fontsize=6)
    axs[1].set_xscale("log"); axs[1].set_xlabel(xl); axs[1].set_ylabel("strongest other lobe, P95 (dB)"); axs[1].set_title("Strongest non-design lobe", fontsize=9); axs[1].legend(fontsize=6)
    save(fig, name)

sens([("phase", "#1f77b4", "phase")], "phase", "constant phase-error half-width (deg)", "6_12_phase_sensitivity.png", f"Phase-trim error (Msub={BM}, {BC}); uniform +-a, 1000 samples")
sens([("delay_uncal", "#d62728", "delay uncalibrated"), ("delay_cal", "#2ca02c", "delay calibrated (slope only)"), ("delay_abs_uncal", "#ff7f0e", "delay, abs. device-length scale, uncal."), ("delay_abs_cal", "#8c564b", "delay, abs. device-length scale, calibrated")], "delay_uncal", "relative delay-error half-width (%)", "6_13_delay_sensitivity.png", "Group-delay error; uniform +-e of channel delay, 1000 samples")
fig, axs = plt.subplots(1, 2, figsize=(10, 3.9))
L = [m["level"] for m in jb["monte_carlo"]["amp"]]
for st, ls in (("median", "-"), ("p95", "--"), ("worst", ":")):
    axs[0].plot(L, [m["max_err"][st] for m in jb["monte_carlo"]["amp"]], color="#9467bd", ls=ls, label=st); axs[1].plot(L, [m["eff_loss_db"][st] for m in jb["monte_carlo"]["amp"]], color="#9467bd", ls=ls, label=st)
axs[0].axhline(jb["nominal_max_err_deg"], color="gray", lw=0.6); axs[0].set_xscale("log"); axs[0].set_xlabel("amplitude error half-width (dB)"); axs[0].set_ylabel("max |pointing error| (deg)"); axs[0].legend(fontsize=7); axs[0].set_title("Amplitude mismatch: pointing", fontsize=9)
axs[1].set_xscale("log"); axs[1].set_xlabel("amplitude error half-width (dB)"); axs[1].set_ylabel("worst-frequency efficiency loss at fixed input power (dB)"); axs[1].set_title("Amplitude mismatch: peak-power efficiency loss", fontsize=9); save(fig, "6_14_amplitude_sensitivity.png")

# 15 Monte Carlo percentiles (combined) --------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(2, 3, figsize=(12, 6.5)); cm = jb["monte_carlo"]["combined"]; xl = [m["level"] for m in cm]
for ax, key, yl in zip(axs.ravel(), ("max_err", "rms_err", "gain_loss_db", "strongest_db", "r_squint"), ("max |error| (deg)", "RMS error (deg)", "worst-f peak-power loss (dB; negative = gain)", "strongest other lobe (dB)", "R_squint")):
    for st, mk in (("median", "o-"), ("p90", "s-"), ("p95", "^-"), ("worst", "x:")): ax.plot(xl, [m[key][st] for m in cm], mk, label=st, ms=4)
    ax.set_title(yl, fontsize=9); ax.set_xlabel("combined level (L1 1deg/0.25%/0.1dB ... L4 10deg/2%/1dB)", fontsize=7)
axs[0, 0].legend(fontsize=7); axs.ravel()[-1].axis("off"); save(fig, "6_15_monte_carlo_percentiles.png")

# 16 tolerance budget --------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 5, figsize=(17, 3.6)); B_ = jb["budget"]
for ax, kind, ttl, un in ((axs[0], "phase", "phase-trim error", "deg"), (axs[1], "delay_uncal", "delay, uncal. (of q tau)", "%"), (axs[2], "delay_cal", "delay, cal. (of q tau)", "%"), (axs[3], "delay_abs_cal", "delay, cal. (of abs. delay)", "%"), (axs[4], "amp", "amplitude error", "dB")):
    vals = []; labs = []
    for lim in ("0.5", "1.0"):
        for st in ("median", "p95"):
            v = B_[kind][st][lim]; vals.append(np.nan if v is None else (np.inf if v == float("inf") else v)); labs.append(f"{st}\n<={lim} deg")
    top = np.nanmax([v for v in vals if np.isfinite(v)] + [1]) * 1.2
    ax.bar(range(4), [top if np.isinf(v) else v for v in vals], color=["#9ecae1", "#3182bd", "#9ecae1", "#3182bd"]); ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=7); ax.set_ylabel(un); ax.set_title(f"{ttl}\nlargest half-width", fontsize=7.5)
    for i, v in enumerate(vals): ax.text(i, (top if np.isinf(v) else v), ">swept max" if np.isinf(v) else ("n/a" if np.isnan(v) else f"{v:.2g}"), ha="center", va="bottom", fontsize=7)
save(fig, "6_16_tolerance_budget.png")

# 17 res30 vs res40 ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 4, figsize=(15, 3.7))
for res, ls in ((30, "-"), (40, "--")):
    axs[0].plot(f, g(res, BM, BC, "theta_peak") - THETA0, ls=ls, label=f"res {res}"); axs[2].plot(f, 10 * np.log10(g(res, BM, BC, "gain")), ls=ls, label=f"res {res}"); axs[3].plot(f, g(res, BM, BC, "strongest_db"), ls=ls, label=f"res {res}")
axs[1].plot(f, g(40, BM, BC, "theta_peak") - g(30, BM, BC, "theta_peak"), color="k"); axs[1].axhline(j["30"][str(BM)]["tol_pt"] / 3, color="gray", ls=":"); axs[1].axhline(-j["30"][str(BM)]["tol_pt"] / 3, color="gray", ls=":")
for a_, (yl, t) in zip(axs, (("deg", "pointing error"), ("deg", "theta(res40) - theta(res30)\n(dotted: +-tol_pt/3)"), ("dB", "peak power"), ("dB", "strongest other lobe"))): deco(a_, yl, t)
axs[0].legend(fontsize=7); save(fig, "6_17_res30_vs_res40.png")
print(sorted(p.name for p in FIG.glob("6_*.png")))
