"""Phase 6B: robustness of the best Phase 6 architecture (selected in 6a by the Revision 7 rule): deterministic perturbation sweeps with first-order
analytic predictions, Monte Carlo (seed 20261010, 1000 samples per tolerance point), tolerance budget, period-quantization scan.
Writes results/phase6/6b_robustness.{npz,json}."""
from _common import P4, RES, np, THETA0, NQ
import json, multiprocessing as mp
from pcttd import array as A, beam as B, subarray as S, robustness as R
from pcttd.io import save_json

SEED = 20261010; NS = 1000; NG = 8001
j6a = json.load(open(RES / "6a_subarrays.json")); best = j6a["best"]; M = best["msub"]; case = best["case"]
assert case in ("trim", "trim_p4", "raw")
z = np.load(P4 / "phase4_handoff.npz", allow_pickle=True)
mk = z["mask_common"]; f = z["f"][mk]; f0 = float(z["f0"]); tau_phc = float(z["dtau_target"]); d_e = 1 / (2 * f0); Nch = z["N"].astype(float)
i0 = int(np.argmin(abs(f - f0))); s0 = np.sin(np.deg2rad(THETA0))
tr = S.tau_required(M, d_e, THETA0); win = S.window_halfwidth(f.max(), M, d_e)
Gtrim, _ = S.retrim(z["H_raw"][:, mk], f, tr, f0)
G_nom = {"trim": Gtrim, "trim_p4": z["H_trim"][:, mk], "raw": z["H_raw"][:, mk]}[case]
U = S.intra_fixed(M, d_e, THETA0, f0, len(f))
fixed_max = j6a["30"][str(M)]["max_abs_err"]["fixed"]

eng_fine = S.Hier(f, NQ, M, d_e, n=40001); E_fine = eng_fine.element_factor(U)
eng = S.Hier(f, NQ, M, d_e, n=NG); E = eng.element_factor(U)
p_fine = eng_fine.power(G_nom, E_fine); th_fine, g_fine, _ = eng_fine.peak(p_fine, THETA0, win)
p_nom = eng.power(G_nom, E); th_nom, g_nom, idx_nom = eng.peak(p_nom, THETA0, win)
js = dict(best=best, seed=SEED, n_samples=NS, grid_points_mc=NG,
          engine_validation=dict(max_dtheta_deg=float(np.abs(th_nom - th_fine).max()), max_dgain=float(np.abs(g_nom - g_fine).max()),
                                 nominal_max_err_deg=float(np.abs(th_fine - THETA0).max())))
assert js["engine_validation"]["max_dtheta_deg"] <= 0.01 and js["engine_validation"]["max_dgain"] <= 0.01
nominal_err = float(np.abs(th_nom - THETA0).max())
q = np.arange(NQ)

# ------------------------------------------------------------------ deterministic sweeps + analytic link (nominal-referenced shifts)
def shift(Gp):
    p = eng.power(Gp, E); th, g, idx = eng.peak(p, THETA0, win)
    return th - th_nom, g, eng.strongest_fast(p, idx)

det = {}
PH = [0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30]; DL = [0.1, 0.25, 0.5, 1, 2, 3, 5, 10]; AM = [0.1, 0.25, 0.5, 1, 2, 3, 5, 6]
for a in PH:                                       # linear phase ramp across the four channels, amplitude +-a deg end to end, LS slope 2a/3 per channel
    eps = np.deg2rad(a) * (2 * q - 3) / 3
    sh, g, st = shift(R.perturb(G_nom, f, eps=eps))
    pred, _ = B.phase_slope_pointing(np.repeat(eps[:, None], len(f), axis=1), f, M * d_e, THETA0)
    det.setdefault("phase", []).append(dict(level=a, shift_max=float(np.abs(sh).max()), shift_f0=float(sh[i0]), pred_max=float(np.abs(pred).max()), pred_f0=float(pred[i0]),
                                           max_pred_minus_num=float(np.abs(pred - sh).max()), gain_min=float(g.min())))
for a in DL:
    for cal, fref in (("uncal", 0.0), ("cal", f0)):
        dtau = a / 100 * q * tau_phc                # every channel late by the same fraction of its relative delay: a pure slope (progressive) error
        sh, g, st = shift(R.perturb(G_nom, f, dtau=dtau, f_ref=fref))
        eps_f = 2 * np.pi * (f[None, :] - fref) * dtau[:, None]
        pred, _ = B.phase_slope_pointing(eps_f, f, M * d_e, THETA0)
        det.setdefault("delay_" + cal, []).append(dict(level=a, shift_max=float(np.abs(sh).max()), shift_f0=float(sh[i0]), shift_low_high=[float(sh[0]), float(sh[-1])],
                                                   pred_max=float(np.abs(pred).max()), pred_low_high=[float(pred[0]), float(pred[-1])],
                                                   max_pred_minus_num=float(np.abs(pred - sh).max()), gain_min=float(g.min())))
for a in AM:                                        # amplitude ramp +-a dB across the channels
    adb = a * (2 * q - 3) / 3
    sh, g, st = shift(R.perturb(G_nom, f, amp_db=adb))
    det.setdefault("amp", []).append(dict(level=a, shift_max=float(np.abs(sh).max()), gain_min=float(g.min()), gain_loss_db=float(10 * np.log10(g / g_nom).min()), strongest_db_max=float(st.max())))
js["deterministic"] = det

# ------------------------------------------------------------------ Monte Carlo tolerance points (child streams of one SeedSequence, fixed order)
tau_rel = q * tau_phc; tau_abs = (Nch + 6) * tau_phc
points = []
for a in PH: points.append(("phase", a, dict(eps_deg=a), 0.0, tau_rel))
for a in DL: points.append(("delay_uncal", a, dict(delay_frac=a / 100), 0.0, tau_rel))
for a in DL: points.append(("delay_cal", a, dict(delay_frac=a / 100), f0, tau_rel))
for a in [0.1, 0.25, 0.5, 1, 2]: points.append(("delay_abs_uncal", a, dict(delay_frac=a / 100), 0.0, tau_abs))
for a in AM: points.append(("amp", a, dict(amp_db=a), 0.0, tau_rel))
COMB = {"L1": (1, 0.25, 0.1), "L2": (2, 0.5, 0.25), "L3": (5, 1.0, 0.5), "L4": (10, 2.0, 1.0)}
for k, (e, dl, am) in COMB.items(): points.append(("combined", k, dict(eps_deg=e, delay_frac=dl / 100, amp_db=am), 0.0, tau_rel))
for a in [0.25, 0.5, 1, 2, 5, 10]: points.append(("delay_abs_cal", a, dict(delay_frac=a / 100), f0, tau_abs))        # appended last: earlier child streams keep their indices
ss = np.random.SeedSequence(SEED).spawn(len(points))


def run(i):
    kind, lev, level, fref, trel = points[i]
    r = R.monte_carlo(eng, E, G_nom, f, THETA0, win, g_nom, fixed_max, ss[i], NS, level, trel, f_ref=fref)
    return i, r


if __name__ == "__main__":
    import time; t0 = time.time()
    with mp.get_context("fork").Pool(min(10, mp.cpu_count())) as pool:
        res = dict(pool.imap_unordered(run, range(len(points))))
    mc = {}
    for i, (kind, lev, level, fref, trel) in enumerate(points):
        r = res[i]
        mc.setdefault(kind, []).append(dict(level=lev, spec={k: float(v) for k, v in level.items()},
                                           **{m: R.percentiles(r[m]) for m in r}))
    js["monte_carlo"] = mc; js["mc_runtime_s"] = time.time() - t0
    js["nominal_max_err_deg"] = nominal_err
    # tolerance budget
    bud = {}
    for kind in ("phase", "delay_uncal", "delay_cal", "delay_abs_uncal", "delay_abs_cal", "amp"):
        L = [m["level"] for m in mc[kind]]
        bud[kind] = {}
        for stat in ("median", "p95"):
            y = [m["max_err"][stat] for m in mc[kind]]
            bud[kind][stat] = {str(lim): R.budget(L, y, lim) for lim in (0.5, 1.0)}
    js["budget"] = bud
    js["combined_L2_p95_max_err"] = [m for m in mc["combined"] if m["level"] == "L2"][0]["max_err"]["p95"]
    js["combined_L2_ok"] = bool(js["combined_L2_p95_max_err"] <= 1.0)
    # Monte Carlo reproducibility spot check (same seed, point 0)
    js["repro_check"] = bool(np.allclose(R.monte_carlo(eng, E, G_nom, f, THETA0, win, g_nom, fixed_max, np.random.SeedSequence(SEED).spawn(len(points))[0], 50, points[0][2], tau_rel=points[0][4])["max_err"],
                                         R.monte_carlo(eng, E, G_nom, f, THETA0, win, g_nom, fixed_max, np.random.SeedSequence(SEED).spawn(len(points))[0], 50, points[0][2], tau_rel=points[0][4])["max_err"]))

    # ------------------------------------------------------------------ period-quantization scan: retrimmed actual PhC channels, k = 1, Msub = 1..12 (analytic rho in 6a)
    quant = []
    for Mq in range(1, 13):
        trq = S.tau_required(Mq, d_e, THETA0); Gq, _ = S.retrim(z["H_raw"][:, mk], f, trq, f0); Uq = S.intra_fixed(Mq, d_e, THETA0, f0, len(f))
        h = S.Hier(f, NQ, Mq, d_e, n=20001); wq = S.window_halfwidth(f.max(), Mq, d_e)
        Gf = S.fixed_channels(f, NQ, trq, f0)
        a_ = h.analyze(Gq, Uq, THETA0, wq); b_ = h.analyze(Gf, Uq, THETA0, wq)
        quant.append(dict(msub=Mq, rho=tau_phc / trq, max_err=float(np.abs(a_["theta_peak"] - THETA0).max()), spread=float(np.ptp(a_["theta_peak"])),
                          R_squint=float(np.abs(a_["theta_peak"] - THETA0).max() / np.abs(b_["theta_peak"] - THETA0).max()),
                          R_spread=float(np.abs(a_["theta_peak"] - a_["theta_peak"][i0]).max() / np.abs(b_["theta_peak"] - THETA0).max()),
                          strongest_db_max=float(a_["strongest_db"].max()), grating_db_max=float(a_["grating_db"].max())))
    js["quantization_scan"] = quant
    # intercept-vs-slope split for Msub=4 (analytic first-order): Phase 4 trim carries f0 (tau_phc - tau_req) per channel step
    dphi = 2 * np.pi * f0 * (tau_phc - tr)
    js["p4_trim_intercept_pred_deg"] = float(np.rad2deg(dphi / (2 * np.pi * f0 * M * d_e) / np.cos(np.deg2rad(THETA0))))
    js["retrim_slope_pred_low_high_deg"] = [float(np.rad2deg(((f0 / ff - 1) * s0 * (1 - tau_phc / tr)) / np.cos(np.deg2rad(THETA0)))) for ff in (f[0], f[-1])]
    js["definitions"] = dict(perturb="G_q -> G_q a_q exp(i eps_q) exp(+i 2 pi (f - f_ref) dtau_q), f_ref = 0 uncalibrated, f0 calibrated",
                             draws="independent uniform per channel (all four for phase and amplitude; delay error e_q * q * tau_PhC so channel 0 has none, 'abs' variant e_q (N_q + 6) tau_PhC)",
                             deterministic="phase/amplitude: linear ramp +-level end to end; delay: every channel scaled by (1 + level) in relative delay",
                             metrics="max/rms |theta_peak - 30 deg| over mask_common; gain loss (positive = loss) = -min_f 10 log10(p/p_nom); eff_loss = same at fixed input power (post hoc); strongest = max over band of the strongest other lobe (dB re peak); R_squint = max_err / max_err(fixed phase)",
                             budget="largest half-width with P95 (or median) of max |pointing error| (incl. nominal) <= limit, linear interpolation; system-level tolerance estimate, not a fabrication specification")
    save_json(RES / "6b_robustness.json", js)
    np.savez(RES / "6b_robustness.npz", f=f, theta_nominal=th_nom, gain_nominal=g_nom)
    print(json.dumps({k: js[k] for k in ("engine_validation", "nominal_max_err_deg", "mc_runtime_s", "budget", "combined_L2_p95_max_err", "combined_L2_ok", "repro_check", "p4_trim_intercept_pred_deg", "retrim_slope_pred_low_high_deg")}, indent=1))
    for kind in det:
        print(kind); [print("  ", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}) for r in det[kind]]
    for kind in mc:
        print(kind); [print("  ", m["level"], "maxerr med/p95/worst", round(m["max_err"]["median"], 3), round(m["max_err"]["p95"], 3), round(m["max_err"]["worst"], 3), "loss p95", round(m["gain_loss_db"]["p95"], 3), "eff", round(m["eff_loss_db"]["p95"], 3), "strong p95", round(m["strongest_db"]["p95"], 2), "R p95", round(m["r_squint"]["p95"], 3)) for m in mc[kind]]
    for r in quant: print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
