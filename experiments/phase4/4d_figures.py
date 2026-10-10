"""Phase 4D: figures from results/phase4/4b_res30.{npz,json} (+ res 40 overlay where marked)."""
from _common import RES, ROOT, np, F0, DTAU_TARGET, VALID, LENGTHS
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pcttd import channels as C

FIG = ROOT / "figures" / "phase4"; FIG.mkdir(parents=True, exist_ok=True)
z = dict(np.load(RES / "4b_res30.npz")); j = json.load(open(RES / "4b_res30.json"))
try:
    z40 = dict(np.load(RES / "4b_res40.npz"))
except FileNotFoundError:
    z40 = None
f = z["f"]
band = j["bands"]["primary_poly_tol10"]
cb = (band["f_low"], band["f_high"])
COL = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
LS = ["-", "--", "-.", ":"]
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130, "axes.spines.top": False, "axes.spines.right": False})
lab = lambda m: f"ch {m} (N = {LENGTHS[m]})"


def deco(ax, band_shade=True, xlim=(0.22, 0.275)):
    ax.axvline(F0, color="k", lw=0.8, alpha=0.6)
    for v in VALID:
        ax.axvline(v, color="gray", lw=0.8, ls=":")
    if band_shade:
        ax.axvspan(*cb, color="#2ca02c", alpha=0.07)
    ax.set_xlim(*xlim); ax.set_xlabel("frequency f (c/a)")


def save(fig, name):
    fig.tight_layout(); fig.savefig(FIG / name); plt.close(fig)


# 1-3 T, R, |H|
for key, name, yl, ttl in (("T", "4a_transmission.png", "T_m = |S21|^2", "Channel transmission (power), raw FDTD res 30"),
                           ("R", "4b_reflection.png", "R_m = |S11|^2", "Channel reflection (power), raw FDTD res 30")):
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for m in range(4):
        ax.plot(f, z[key][m], color=COL[m], ls=LS[m], label=lab(m))
    deco(ax); ax.set_ylabel(yl); ax.set_title(ttl)
    if key == "R":
        ax.set_yscale("log")
    ax.legend(ncol=2, fontsize=8); save(fig, name)
fig, ax = plt.subplots(figsize=(6.2, 3.6))
for m in range(4):
    ax.plot(f, np.abs(z["H"][m]), color=COL[m], ls=LS[m], label=lab(m))
deco(ax); ax.set_ylabel("|H_m| (amplitude)"); ax.set_title("Complex transfer-function amplitude |H_m| = |S21|"); ax.set_ylim(0.97, 1.03); ax.legend(ncol=2, fontsize=8)
save(fig, "4c_amplitude.png")

# 4 relative phase (cycles) minus the ideal 2 pi f m dtau at... show phase itself in cycles
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
for m in range(4):
    axs[0].plot(f, z["phi_rel"][m] / (2 * np.pi), color=COL[m], ls=LS[m], label=lab(m) if m else "ch 0 (reference)")
deco(axs[0]); axs[0].set_ylabel("unwrapped arg[H_m / H_0] / 2 pi (cycles)"); axs[0].set_title("Relative phase"); axs[0].legend(fontsize=8)
for m in range(1, 4):
    axs[1].plot(f, z["phi_rel"][m] - 2 * np.pi * f * m * DTAU_TARGET - 2 * np.pi * np.round((z["phi_rel"][m][np.argmin(abs(f - F0))] - 2 * np.pi * F0 * m * DTAU_TARGET) / (2 * np.pi)),
                color=COL[m], ls=LS[m], label=lab(m))
deco(axs[1]); axs[1].set_ylabel("arg[H_m/H_0] - 2 pi f m dtau_target (rad, branch at f0)"); axs[1].set_title("Relative phase minus ideal TTD (= raw phase error)"); axs[1].legend(fontsize=8)
save(fig, "4d_relative_phase.png")

# 5 measured delay by channel
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
for m in range(4):
    axs[0].plot(f, z["tau_abs"][m], color=COL[m], ls=LS[m], lw=0.8, label=lab(m))
deco(axs[0]); axs[0].set_ylabel("absolute delay (a/c), central difference"); axs[0].set_title("Total channel delay (N_m + 6 periods)"); axs[0].legend(fontsize=8)
for m in range(1, 4):
    axs[1].plot(f, z["Dtau_poly"][m], color=COL[m], ls=LS[m], label=lab(m) + " minus ch 0 (local quadratic)")
    axs[1].axhline(m * DTAU_TARGET, color=COL[m], lw=0.6, alpha=0.6)
deco(axs[1]); axs[1].set_ylabel("relative delay Dtau_m (a/c)"); axs[1].set_title("Delay relative to ch 0 (thin lines: m dtau_target)"); axs[1].legend(fontsize=7)
save(fig, "4e_delay_by_channel.png")

# 6 adjacent incremental delay
fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
for m in range(3):
    ax = axs[m]
    ax.plot(f, z["inc"][m], color="0.7", lw=0.7, label="central difference (pre-declared)")
    ax.plot(f, z["inc_poly"][m], color=COL[m + 1], lw=1.5, label="local quadratic (primary)")
    ax.axhline(DTAU_TARGET, color="k", lw=0.8); ax.axhspan(0.9 * DTAU_TARGET, 1.1 * DTAU_TARGET, color="k", alpha=0.06)
    deco(ax); ax.set_title(f"increment {m} -> {m + 1}"); ax.set_ylim(3.0, 5.0)
axs[0].set_ylabel("adjacent delay step (a/c)"); axs[0].legend(fontsize=7)
save(fig, "4f_incremental_delay.png")

# 7 incremental delay error
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
for m in range(3):
    axs[0].plot(f, (z["inc_poly"][m] / DTAU_TARGET - 1) * 100, color=COL[m + 1], ls=LS[m + 1], label=f"{m} -> {m + 1}")
axs[0].axhline(10, color="k", lw=0.6); axs[0].axhline(-10, color="k", lw=0.6)
deco(axs[0]); axs[0].set_ylabel("increment error vs target (%) [local quadratic]"); axs[0].set_title("Adjacent-increment error"); axs[0].legend(fontsize=8)
for m in range(1, 4):
    axs[1].plot(f, z["eps_tau_rel"][m] * 100 if False else (z["Dtau_poly"][m] / (m * DTAU_TARGET) - 1) * 100, color=COL[m], ls=LS[m], label=lab(m))
deco(axs[1]); axs[1].set_ylabel("(Dtau_m - m dtau_target) / (m dtau_target) (%)"); axs[1].set_title("Delay error relative to ch 0"); axs[1].legend(fontsize=8)
save(fig, "4g_incremental_delay_error.png")

# 8-9 phase errors
for key, name, ttl in (("eps_phi_raw", "4h_phase_error_raw.png", "Raw PhC phase error vs ideal TTD"), ("eps_phi_trim", "4i_phase_error_trimmed.png", "Phase-intercept-trimmed PhC phase error vs ideal TTD")):
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for m in range(1, 4):
        ax.plot(f, z[key][m], color=COL[m], ls=LS[m], label=lab(m))
    deco(ax); ax.set_ylabel("phase error (rad)"); ax.set_title(ttl)
    if key == "eps_phi_trim":
        ax.axhspan(-np.pi / 8, np.pi / 8, color="k", alpha=0.05); ax.text(0.2205, np.pi / 8 * 1.05, "pi/8 reporting threshold", fontsize=7)
        ax.set_ylim(-0.45, 0.45)
        if z40 is not None:
            for m in range(1, 4):
                ax.plot(f, z40[key][m], color=COL[m], lw=0.6, alpha=0.5)
    ax.legend(fontsize=8); save(fig, name)

# 10 amplitude imbalance
fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.plot(f, z["imbalance_db"], color=COL[0], label="res 30")
if z40 is not None:
    ax.plot(f, z40["imbalance_db"], color=COL[1], lw=0.8, label="res 40")
ax.axhline(1.0, color="k", lw=0.8, ls="--", label="1 dB criterion")
deco(ax); ax.set_ylabel("20 log10(max_m|H_m| / min_m|H_m|) (dB, amplitude)"); ax.set_title("Channel amplitude imbalance"); ax.set_yscale("log"); ax.set_ylim(1e-3, 3); ax.legend(fontsize=8)
save(fig, "4j_amplitude_imbalance.png")

# 11 common usable band
fig, ax = plt.subplots(figsize=(7, 3.6))
names = [("domain", "validated domain [0.231, 0.268]"), ("T", "T_m >= 0.5"), ("unitarity", "|defect| <= 0.02"), ("residual", "projection residual <= 0.10"),
         ("increment", "increments within +-10% (local quadratic)"), ("amplitude", "imbalance <= 1 dB")]
m_all, _ = (None, None)
import importlib.util, sys
sys.argv = ["x"]
spec = importlib.util.spec_from_file_location("b4", str(ROOT / "experiments" / "phase4" / "4b_analyze.py")); b4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(b4)
mk, allm = b4.masks(z, f, inc_key="inc_poly", inc_tol=0.10)
mk_c, allm_c = b4.masks(z, f, inc_key="inc", inc_tol=0.10)
rows = [(lab_, mk[k]) for k, lab_ in names] + [("increments +-10% (central, pre-declared)", mk_c["increment"]), ("COMMON band (primary)", allm), ("COMMON band (central difference)", allm_c)]
for i, (lab_, mask) in enumerate(rows[::-1]):
    ax.fill_between(f, i - 0.35, i + 0.35, where=mask, step="mid", color="#2ca02c" if "COMMON" in lab_ else "#1f77b4", alpha=0.8)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=7)
deco(ax, band_shade=False); ax.set_title("Per-criterion pass masks and common usable band")
save(fig, "4k_common_band.png")

# 12 cascade
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
for k, (key, col) in enumerate(zip(j["cascade"], COL)):
    c = j["cascade"][key]["target_N"]; ci = LENGTHS.index(c)
    S = z[f"casc_S21_{key}"]; Hd = z["H"][ci]
    axs[0].plot(f, np.abs(S) / np.abs(Hd) - 1, color=col, label=f"{key.replace('_to_', ' -> ')}")
    axs[1].plot(f, np.angle(S / Hd), color=col, label=f"{key.replace('_to_', ' -> ')}")
for ax in axs:
    deco(ax); ax.legend(fontsize=8)
axs[0].set_ylabel("|S21 cascade| / |S21 direct| - 1"); axs[1].set_ylabel("arg(S21 cascade / S21 direct) (rad)")
fig.suptitle("Cascade (T-matrix) prediction vs directly simulated device", fontsize=10)
save(fig, "4l_cascade_validation.png")
print("figures written", sorted(p.name for p in FIG.glob("*.png")))
