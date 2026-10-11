"""Analytic additions B: one-carrier-cycle propagation length for (1) vacuum/air, (2) an IDEAL homogeneous eps = 12 medium (nondispersive analytic
benchmark; NOT a modeled dielectric waveguide), (3) the selected plain W1 PhC.  Pure arithmetic on project values; no simulation.
Writes results/final/conventional_reference.json."""
from _common import RES, nload, np
from pcttd.io import save_json

z = nload("phase4/phase4_handoff.npz")
f0 = float(z["f0"]); tau = float(z["dtau_target"])          # tau in a/c per period; c = a = 1
lam0 = 1.0 / f0                                             # in units of a
n_g_phc = tau                                               # group index = c tau / a (per-period delay / period)

def one_cycle(n_g):
    L_over_lam0 = 1.0 / n_g                                 # L n_g / c = 1/f0  ->  L = lambda0 / n_g
    return dict(n_g=n_g, L_over_lambda0=L_over_lam0, L_over_a=L_over_lam0 * lam0)

ref = {"air": one_cycle(1.0), "uniform_eps12": one_cycle(np.sqrt(12.0)), "W1_PhC": one_cycle(n_g_phc)}
# PhC cross-check by the period route: delay per period = f0 tau cycles -> periods for one cycle = 1/(f0 tau); length = that many a
periods_one_cycle = 1.0 / (f0 * tau)
ref["W1_PhC"]["periods_for_one_cycle"] = periods_one_cycle
ref["W1_PhC"]["L_over_a_period_route"] = periods_one_cycle
ref["W1_PhC"]["delay_per_period_cycles"] = f0 * tau
ref["W1_PhC"]["a_over_lambda0"] = f0
for k, r in ref.items():
    r["relative_length_vs_W1"] = r["L_over_lambda0"] / ref["W1_PhC"]["L_over_lambda0"]
    r["length_reduction_of_W1_vs_this_pct"] = 100 * (1 - ref["W1_PhC"]["L_over_lambda0"] / r["L_over_lambda0"])
ref["air"].update(intrinsic_dispersion="none (assumed)", needs_constant_phase_trim=False, delay_selectable_by_length="continuous", integer_period_quantum=False)
ref["uniform_eps12"].update(intrinsic_dispersion="none (assumed: ideal nondispersive homogeneous medium, n_g = n_phase = sqrt(12))", needs_constant_phase_trim=False,
                            delay_selectable_by_length="continuous", integer_period_quantum=False)
ref["W1_PhC"].update(intrinsic_dispersion="n_g flat to +-10 % over 0.231-0.268 c/a (14.8 % fractional bandwidth); nonzero GVD", needs_constant_phase_trim=True,
                     delay_selectable_by_length="only in whole lattice periods in the simulated channels (fractional-period lengths not simulated)", integer_period_quantum=True)
out = dict(f0=f0, lambda0_over_a=lam0, tau_phc_a_over_c=tau, references=ref,
           caveats=["uniform eps = 12 is an IDEAL HOMOGENEOUS-MEDIUM propagation reference, not a modeled dielectric waveguide",
                    "a real dielectric waveguide's effective and group index are set by its geometry; n_g != sqrt(eps) in general",
                    "length only: lateral width, cladding, couplers, loss and bends are not compared; no FDTD of any conventional line was run",
                    "the comparison limits the footprint framing only; it does not affect the beamforming result"],
           consistency=dict(L_phc_two_routes_abs_diff=abs(ref["W1_PhC"]["L_over_a"] - periods_one_cycle),
                            uniform_eps12_L_over_lambda0_times_sqrt12=ref["uniform_eps12"]["L_over_lambda0"] * np.sqrt(12.0)))
save_json(RES / "conventional_reference.json", out)
for k, r in ref.items():
    print(f"{k:14s} n_g {r['n_g']:.4f}  L/lam0 {r['L_over_lambda0']:.4f}  L/a {r['L_over_a']:.4f}  rel to W1 {r['relative_length_vs_W1']:.3f}  W1 shorter by {r['length_reduction_of_W1_vs_this_pct']:.1f} %")
print("periods for one cycle", periods_one_cycle, "abs diff routes", out["consistency"]["L_phc_two_routes_abs_diff"])
