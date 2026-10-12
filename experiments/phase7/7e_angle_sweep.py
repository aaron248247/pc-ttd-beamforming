"""Analytic additions A: steering-angle sweep (Msub = 4) with the EXISTING Phase 4 transfer functions and the Phase 6 hierarchical model, plus the analytic
(Msub, k) match-angle table.  No FDTD.  Writes results/final/angle_sweep.{json,npz} and figures/final/fig10_angle_sweep.*, fig11_match_map.*.

Cases per steering angle theta0 (fixed Phase 4 delay slope, only the constant intercept trim is reset):
  fixed      constant-phase shifters set at f0 (baseline)
  ideal_sub  ideal TTD between subarrays, fixed phase inside (the delay-quantum-free reference)
  trim       Phase 4 raw H_m, one constant trim per channel chosen at f0 for the REQUIRED step Msub d_e sin(theta0)  (Phase 6 'trim')
  ideal_elem ideal element-level TTD (HPBW reference for the Phase 6 5 %-HPBW tolerance)
"""
from _common import RES, FIG, nload, np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pcttd import subarray as S
from pcttd.io import save_json

MSUB, NQ = 4, 4
z = nload("phase4/phase4_handoff.npz")
fa = z["f"]; mk = z["mask_common"]; f = fa[mk]
f0 = float(z["f0"]); tau_phc = float(z["dtau_target"]); lam0 = 1 / f0; d_e = lam0 / 2
i0 = int(np.argmin(abs(f - f0))); fmax = f.max(); I = [0, i0, len(f) - 1]
assert abs(f[i0] - f0) < 1e-9
Hraw = z["H_raw"][:, mk]
q_phc = f0 * tau_phc                                  # PhC delay per period in carrier cycles at f0

def q_req(msub, th):                                  # required cycles per adjacent subarray step at f0: f0 msub d_e sin(th) = msub sin(th) / 2
    return f0 * S.tau_required(msub, d_e, th)

def match_sin(msub, k):                               # sin(theta_match) = 2 k q_phc / msub
    return 2 * k * q_phc / msub

th_match = float(np.rad2deg(np.arcsin(match_sin(MSUB, 1))))
grid = np.unique(np.round(np.r_[np.arange(5.0, 70.5, 1.0), 30.0, th_match], 6))
eng = S.Hier(f, NQ, MSUB, d_e)                        # engine is independent of theta0
win = S.window_halfwidth(fmax, MSUB, d_e)

keys = ["fixed", "ideal_sub", "trim", "ideal_elem"]
A = {k: [] for k in ["theta0", "req_cycles", "phc_cycles", "rho", "mismatch_pct", "max_err_trim", "rms_err_trim", "max_err_fixed", "rms_err_fixed", "max_err_ideal_sub",
                     "R_squint", "R_squint_spread", "spread_trim", "vs_ideal_sub_max", "gain_min", "gain_min_db", "strongest_db_max", "grating_db_max", "hpbw_f0_trim",
                     "hpbw_ideal_f0", "tol_5pct_hpbw", "pred_dev_max", "pred_err_max", "num_err_max_signed", "pred_valid", "window_use_frac", "n3db_max",
                     "err_trim_low", "err_trim_f0", "err_trim_high", "strongest_db_fixed", "strongest_db_ideal_sub", "bias_f0_trim"]}
TH = {}
for th0 in grid:
    s0 = np.sin(np.deg2rad(th0)); tr = S.tau_required(MSUB, d_e, th0)
    rho, mm, _ = S.mismatch(tau_phc, tr, f0)
    Gtrim, _ = S.retrim(Hraw, f, tr, f0)
    G = dict(fixed=S.fixed_channels(f, NQ, tr, f0), ideal_sub=S.ideal_channels(f, NQ, tr), ideal_elem=S.ideal_channels(f, NQ, tr), trim=Gtrim)
    Uf = S.intra_fixed(MSUB, d_e, th0, f0, len(f))
    U = dict(fixed=Uf, ideal_sub=Uf, trim=Uf, ideal_elem=S.intra_ttd(MSUB, d_e, th0, f))
    R = {k: eng.analyze(G[k], U[k], th0, win) for k in keys}
    err = {k: R[k]["theta_peak"] - th0 for k in keys}
    fixed_max = np.abs(err["fixed"]).max()
    # analytic subarray-level squint law (as used in Phase 6 6a_subarrays.py): sin(theta_pk) - sin(theta0) = (f0/f - 1) sin(theta0) (1 - rho)
    arg = s0 + (f0 / f - 1) * s0 * (1 - rho)
    valid = bool(np.all(np.abs(arg) < 1))
    pred = np.rad2deg(np.arcsin(np.clip(arg, -1, 1)))
    row = dict(theta0=th0, req_cycles=f0 * tr, phc_cycles=q_phc, rho=rho, mismatch_pct=100 * mm,
               max_err_trim=np.abs(err["trim"]).max(), rms_err_trim=np.sqrt(np.mean(err["trim"] ** 2)),
               max_err_fixed=fixed_max, rms_err_fixed=np.sqrt(np.mean(err["fixed"] ** 2)), max_err_ideal_sub=np.abs(err["ideal_sub"]).max(),
               R_squint=np.abs(err["trim"]).max() / fixed_max,
               R_squint_spread=np.abs(R["trim"]["theta_peak"] - R["trim"]["theta_peak"][i0]).max() / fixed_max,
               spread_trim=np.ptp(R["trim"]["theta_peak"]),
               vs_ideal_sub_max=np.abs(R["trim"]["theta_peak"] - R["ideal_sub"]["theta_peak"]).max(),
               gain_min=R["trim"]["gain"].min(), gain_min_db=10 * np.log10(R["trim"]["gain"].min()),
               strongest_db_max=R["trim"]["strongest_db"].max(), grating_db_max=R["trim"]["grating_db"].max(),
               hpbw_f0_trim=R["trim"]["hpbw"][i0], hpbw_ideal_f0=R["ideal_elem"]["hpbw"][i0], tol_5pct_hpbw=0.05 * R["ideal_elem"]["hpbw"][i0],
               pred_dev_max=np.abs(pred - R["trim"]["theta_peak"]).max(), pred_err_max=np.abs(pred - th0).max(),
               num_err_max_signed=err["trim"][np.argmax(np.abs(err["trim"]))], pred_valid=valid,
               window_use_frac=np.abs(np.sin(np.deg2rad(R["trim"]["theta_peak"])) - s0).max() / win, n3db_max=R["trim"]["n_3db"].max(),
               err_trim_low=err["trim"][I[0]], err_trim_f0=err["trim"][I[1]], err_trim_high=err["trim"][I[2]],
               strongest_db_fixed=R["fixed"]["strongest_db"].max(), strongest_db_ideal_sub=R["ideal_sub"]["strongest_db"].max(), bias_f0_trim=err["trim"][i0])
    for k, v in row.items():
        A[k].append(v)
    TH[float(th0)] = dict(trim=R["trim"]["theta_peak"], fixed=R["fixed"]["theta_peak"], ideal_sub=R["ideal_sub"]["theta_peak"], pred=pred)
A = {k: np.array(v) for k, v in A.items()}
th = A["theta0"]

# ------------------------------------------------------------------ windows (post hoc interpretive thresholds, NOT gates)
def intervals(mask):
    out = []; i = 0
    while i < len(mask):
        if mask[i]:
            j = i
            while j + 1 < len(mask) and mask[j + 1]:
                j += 1
            out.append([float(th[i]), float(th[j])]); i = j + 1
        else:
            i += 1
    return out
crit = {"A_Rsquint_le_0.20": A["R_squint"] <= 0.20,
        "B_vs_ideal_sub_le_5pct_HPBW": A["vs_ideal_sub_max"] <= A["tol_5pct_hpbw"],
        "C_max_err_le_0.5deg": A["max_err_trim"] <= 0.5,
        "D_Rsquint_lt_1_beats_fixed": A["R_squint"] < 1.0}
windows = {k: dict(intervals=intervals(v), n_pass=int(v.sum()), n_total=len(v)) for k, v in crit.items()}
both = crit["A_Rsquint_le_0.20"] & crit["B_vs_ideal_sub_le_5pct_HPBW"]
windows["A_and_B"] = dict(intervals=intervals(both), n_pass=int(both.sum()), n_total=len(both))
allc = both & crit["C_max_err_le_0.5deg"]
windows["A_and_B_and_C"] = dict(intervals=intervals(allc), n_pass=int(allc.sum()), n_total=len(allc))

# minimum-mismatch angle on the numerical grid (consistency with analytic exact match)
i_min_mm = int(np.argmin(np.abs(A["mismatch_pct"]))); i_min_err = int(np.argmin(A["max_err_trim"])); i_min_vs = int(np.argmin(A["vs_ideal_sub_max"]))

# ------------------------------------------------------------------ match table
MS = [1, 2, 4, 8, 16]
table = []
for m in MS:
    ks = []; k = 1
    while match_sin(m, k) <= 1 + 1e-12:
        sn = match_sin(m, k)
        ks.append(dict(k=k, sin_theta=sn, theta_deg=float(np.rad2deg(np.arcsin(min(sn, 1.0)))), subarray_pitch_lambda0=m / 2,
                       q_req_at_match=m * sn / 2, delay_cycles_per_step=k * q_phc, n_subarray_replicas_visible_f0=int(2 * np.floor((1 + abs(sn)) * m / 2 + 0) if False else 0)))
        k += 1
    # first k that would not fit: required sin > 1
    table.append(dict(msub=m, k_max=len(ks), matches=ks, sin_k1=match_sin(m, 1)))
# (re-computed replica count below to keep the table entries simple)
for t in table:
    for e in t["matches"]:
        sp = e["sin_theta"]; pitch = t["msub"] / 2                       # in lambda0 at f0
        qq = np.arange(-40, 41); qq = qq[qq != 0]
        e["n_replicas_visible_f0"] = int(np.sum(np.abs(sp + qq / pitch) <= 1 + 1e-12))
        e.pop("n_subarray_replicas_visible_f0", None)

# ------------------------------------------------------------------ summary
i30 = int(np.argmin(abs(th - 30.0)))
summ = dict(theta_match_exact_deg=th_match, sin_theta_match=match_sin(MSUB, 1), q_phc=q_phc, tau_phc=tau_phc, f0=f0,
            grid_deg=[float(th[0]), float(th[-1])], n_grid=len(th),
            theta_min_abs_mismatch_grid=float(th[i_min_mm]), theta_min_max_err_grid=float(th[i_min_err]), theta_min_vs_ideal_sub_grid=float(th[i_min_vs]),
            max_err_min=float(A["max_err_trim"][i_min_err]), at_30=dict(max_err_trim=float(A["max_err_trim"][i30]), R_squint=float(A["R_squint"][i30]), rho=float(A["rho"][i30])),
            max_err_trim_range=[float(A["max_err_trim"].min()), float(A["max_err_trim"].max())], rms_err_trim_range=[float(A["rms_err_trim"].min()), float(A["rms_err_trim"].max())],
            analytic_law=dict(max_dev_vs_numeric_deg=float(np.nanmax(A["pred_dev_max"])), theta_at_max_dev=float(th[np.argmax(A["pred_dev_max"])]),
                              median_dev_deg=float(np.median(A["pred_dev_max"])), dev_le_0p05_deg_interval=intervals(A["pred_dev_max"] <= 0.05),
                              dev_le_0p05_over_err_ratio_note="pred_dev_max vs max_err_trim in arrays"),
            windows=windows, thresholds=dict(A="R_squint <= 0.20", B="max|theta_PhC - theta_ideal_subarray| <= 0.05 HPBW_ideal_elem(f0)", C="max |theta_PhC - theta0| <= 0.5 deg",
                                             status="post hoc interpretive thresholds, not gates"),
            bandwidth=dict(common_band=[float(f[0]), float(fmax)], n_freq=len(f)))
js = dict(summary=summ, match_table=table,
          sweep={k: A[k] for k in A},
          definitions=dict(trim="constant phase per channel at f0 for the REQUIRED step q*tau_req(theta0) (pcttd.subarray.retrim) applied to Phase 4 H_raw; delay slope untouched",
                           ideal_sub="ideal TTD between subarrays, fixed-phase intra-subarray weights (Phase 6 ideal_sub)",
                           analytic_law="sin(theta_pk) - sin(theta0) = (f0/f - 1) sin(theta0) (1 - rho), rho = q_phc/q_req (same as Phase 6 6a_subarrays.py)",
                           match="Msub sin(theta)/2 = k f0 tau_PhC; k = number of PhC periods per adjacent subarray delay increment (k > 1 not simulated)",
                           power="|AF|^2/N_el^2; window = 1/(2 fmax Msub d_e); 40001-pt sin grid; band 0.231-0.268"))
save_json(RES / "angle_sweep.json", js)
np.savez(RES / "angle_sweep.npz", f=f, f0=f0, tau_phc=tau_phc, q_phc=q_phc, **{f"sw_{k}": v for k, v in A.items()},
         **{f"theta_peak_{c}_{int(round(t * 1000)):05d}": TH[t][c] for t in TH for c in ("trim", "fixed", "ideal_sub", "pred")})

# ------------------------------------------------------------------ figures
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
fig, ax = plt.subplots(2, 2, figsize=(10, 7.2))
a = ax[0, 0]; a.plot(th, A["req_cycles"], "k-", label=r"required $M_{sub}\sin\theta_0/2$"); a.axhline(q_phc, color="C0", label=f"PhC 1 period = {q_phc:.4f} cycle"); a.axvline(th_match, color="C3", ls=":", label=f"exact match {th_match:.2f}°")
a.set_xlabel(r"steering angle $\theta_0$ (deg)"); a.set_ylabel("delay step at $f_0$ (carrier cycles)"); a.set_title("(a) required vs available delay ($M_{sub}=4$)"); a.legend(fontsize=7, loc="lower right", frameon=False)
a = ax[0, 1]; a.plot(th, A["mismatch_pct"], "C0-", label=r"mismatch $(\rho-1)$ (%)"); a.axhline(0, color="k", lw=0.6); a.axvline(th_match, color="C3", ls=":"); a.axvline(30, color="gray", ls="--", lw=0.8)
a.set_ylim(-60, 60); a.set_xlabel(r"$\theta_0$ (deg)"); a.set_ylabel("delay mismatch (%)"); a.set_title(r"(b) mismatch, $\rho=q_{PhC}/q_{req}$"); a2 = a.twinx(); a2.plot(th, A["rho"], "C1--", label=r"$\rho$"); a2.set_ylim(0, 6); a2.set_ylabel(r"$\rho$", color="C1"); a2.grid(False)
a = ax[1, 0]; a.semilogy(th, A["max_err_fixed"], "C3-", label="fixed phase"); a.semilogy(th, A["max_err_trim"], "C0-", label="re-trimmed PhC, max |err| vs $\\theta_0$"); a.semilogy(th, A["vs_ideal_sub_max"], "C2--", label="re-trimmed PhC vs ideal subarray TTD")
a.semilogy(th, A["max_err_ideal_sub"], "k-.", lw=0.9, label="ideal subarray TTD vs $\\theta_0$"); a.semilogy(th, A["tol_5pct_hpbw"], "C4:", label="5 % HPBW"); a.axhline(0.5, color="gray", ls=":", lw=0.8); a.text(5.5, 0.55, "0.5°", fontsize=7)
a.set_xlabel(r"$\theta_0$ (deg)"); a.set_ylabel("max pointing error over 0.231–0.268 (deg)"); a.set_title("(c) pointing error"); a.legend(fontsize=6.6, loc="upper right", frameon=False)
a = ax[1, 1]; a.semilogy(th, A["R_squint"], "C0-", label=r"$R_{squint}$ (re-trimmed PhC)"); a.semilogy(th, A["R_squint_spread"], "C0--", lw=0.9, label="offset-removed spread"); a.axhline(0.2, color="gray", ls=":"); a.axhline(1.0, color="C3", ls=":", lw=0.8)
a.set_xlabel(r"$\theta_0$ (deg)"); a.set_ylabel("$R_{squint}$ (rel. fixed phase)"); a.set_title(r"(d) squint ratio; dotted: 0.2 and 1"); a.legend(fontsize=7, loc="upper right", frameon=False)
for x in ax.flat:
    x.set_xlim(5, 70)
fig.tight_layout(); fig.savefig(FIG / "fig10_angle_sweep.png", dpi=200); fig.savefig(FIG / "fig10_angle_sweep.pdf"); plt.close(fig)

fig, a = plt.subplots(figsize=(6.4, 3.6))
for j, t in enumerate(table):
    for e in t["matches"]:
        a.plot(e["theta_deg"], j, "o", color=f"C{e['k'] % 10}", ms=6)
        if e["theta_deg"] < 85:
            a.annotate(f"k={e['k']}", (e["theta_deg"], j), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=6.5)
a.axvspan(0, 70, color="0.93", zorder=0); a.axvline(30, color="gray", ls="--", lw=0.8); a.set_yticks(range(len(table))); a.set_yticklabels([f"{t['msub']}" for t in table])
a.set_xlabel(r"exact-match steering angle $\theta_{match}$ (deg)"); a.set_ylabel(r"$M_{sub}$ (elements per subarray)"); a.set_xlim(0, 90); a.set_ylim(-0.5, len(table) - 0.3)
a.set_title(r"$\sin\theta_{match}=2k f_0\tau_{PhC}/M_{sub}$ (shaded: swept range)", fontsize=9)
fig.tight_layout(); fig.savefig(FIG / "fig11_match_map.png", dpi=200); fig.savefig(FIG / "fig11_match_map.pdf"); plt.close(fig)

# ------------------------------------------------------------------ console summary
print("exact match theta", th_match, "sin", match_sin(MSUB, 1), "q_phc", q_phc)
print("grid min |mismatch|", th[i_min_mm], " min max_err", th[i_min_err], A["max_err_trim"][i_min_err], " min vs ideal sub", th[i_min_vs])
print(json.dumps({k: v for k, v in summ.items() if k in ("windows", "max_err_trim_range", "rms_err_trim_range", "analytic_law", "at_30")}, indent=1))
for t in table:
    print(t["msub"], [(e["k"], round(e["theta_deg"], 2), e["n_replicas_visible_f0"]) for e in t["matches"]])
print(" th  req   rho  mm%   maxE   rms   Rsq  vsIdeal  tol  gmin_dB  strong  grat  hp   predDev")
for i in range(0, len(th)):
    if i % 3 == 0 or abs(th[i] - th_match) < 1e-3 or th[i] == 30:
        print(f"{th[i]:6.2f} {A['req_cycles'][i]:.3f} {A['rho'][i]:5.2f} {A['mismatch_pct'][i]:7.1f} {A['max_err_trim'][i]:7.3f} {A['rms_err_trim'][i]:6.3f} {A['R_squint'][i]:6.3f} {A['vs_ideal_sub_max'][i]:7.3f} {A['tol_5pct_hpbw'][i]:5.2f} {A['gain_min_db'][i]:6.2f} {A['strongest_db_max'][i]:6.1f} {A['grating_db_max'][i]:6.1f} {A['hpbw_f0_trim'][i]:5.1f} {A['pred_dev_max'][i]:6.3f}")
