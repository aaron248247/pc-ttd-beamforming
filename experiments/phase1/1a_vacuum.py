"""Phase 1A: vacuum propagation. Expected: H = E(z+L)/E(z) = exp(+i k L), group delay L/c (+ known Yee dispersion)."""
from dataclasses import replace
import numpy as np
from _common import *  # noqa
from pcttd import fdtd1d as F1, extract as X, io as pio

BAND = (0.16, 0.44)          # analysis band (edges trimmed from the simulated 0.15-0.45 to avoid one-sided derivative ends)
BASE = F1.Cfg1D()            # res 40, dpml 2, L = 10 (monitors 0 and 10), source at -5, decay 1e-10, decimation 1


def analyze(r, L, n=1.0):
    H = r["E"][1] / r["E"][0]
    _, tau_y = X.yee_1d_delay_phase(r["f"], L, n, r["dx"], r["dt"])
    m = F1.delay_metrics(r["f"], H, L * n, tau_y, BAND)
    amp = np.abs(r["E"][0])
    m.update(dt=r["dt"], t_end=r["t_end"], hit_max_time=r["hit_max_time"],
             coverage=float(amp[X.band_mask(r["f"], *BAND)].min() / amp.max()),
             abs_H_dev=float(np.max(np.abs(np.abs(H[X.band_mask(r["f"], *BAND)]) - 1.0))),
             freq_axis_max_dev=float(np.max(np.abs(r["f"] - np.array(BASE.freqs)))) if len(r["f"]) == len(BASE.freqs) else None)
    return m


def study(name, param, values, make_cfg, L_of=lambda v: 10.0):
    cfgs = [make_cfg(v) for v in values]
    rs = pool_map(F1.run_1d, cfgs)
    rows = []
    for v, c, r in zip(values, cfgs, rs):
        row = dict(param=param, value=v)
        row.update(analyze(r, L_of(v)))
        rows.append(row)
    return rows


if __name__ == "__main__":
    out = dict(config=cfg_dict(BASE), band=BAND, expected="tau = L/c (+ Yee numerical dispersion, second order)")
    base = pool_map(F1.run_1d, [BASE])[0]
    out["base"] = analyze(base, 10.0)
    np.savez_compressed(RES / "1a_base_curves.npz", f=base["f"], E=base["E"], dt=base["dt"], dx=base["dx"], L=10.0)

    # sign convention, tested numerically: the correct sign gives a POSITIVE delay equal to L/c
    H = base["E"][1] / base["E"][0]
    t_pos, _ = X.group_delay_phase(base["f"], H, sign=+1)
    t_neg, _ = X.group_delay_phase(base["f"], H, sign=-1)
    mk = X.band_mask(base["f"], *BAND)
    out["sign_test"] = dict(mean_tau_sign_plus=float(np.mean(t_pos[mk])), mean_tau_sign_minus=float(np.mean(t_neg[mk])),
                            expected_L=10.0, established_sign=+1 if abs(np.mean(t_pos[mk]) - 10.0) < abs(np.mean(t_neg[mk]) - 10.0) else -1,
                            meep_sign_constant=X.MEEP_SIGN)

    # 1. resolution
    res_rows = study("res", "res", [10, 20, 40, 80, 160], lambda v: replace(BASE, res=v))
    h = [1.0 / r["value"] for r in res_rows]
    out["resolution"] = dict(rows=res_rows,
                             order_mean_err=X.observed_order(h, [abs(r["mean_err"]) for r in res_rows])[0],
                             order_max_err=X.observed_order(h, [r["max_abs_err"] for r in res_rows])[0])
    # 2. PML thickness (interior domain fixed)
    out["pml"] = study("pml", "dpml", [0.25, 0.5, 1.0, 2.0, 4.0], lambda v: replace(BASE, dpml=v))
    # 3. decay criterion
    out["decay"] = study("decay", "decay_tol", [1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12], lambda v: replace(BASE, decay_tol=v))
    # 4. source spectral coverage and DFT decimation
    out["source_width"] = study("fw", "src_fwidth", [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6], lambda v: replace(BASE, src_fwidth=v))
    out["decimation"] = study("dec", "decimation(0=auto)", [0, 1, 4], lambda v: replace(BASE, decimation=v))
    # 5. monitor / reference-plane placement (L fixed at 10)
    out["placement_source_distance"] = study("ds", "source->monitor1 distance", [0.25, 0.5, 1.0, 2.0, 4.0, 8.0],
                                             lambda v: replace(BASE, z_src=-v, z_lo=min(-8.0, -v - 2.0)))
    out["placement_pml_distance"] = study("dr", "monitor2->PML-inner distance", [0.25, 0.5, 1.0, 2.0, 4.0],
                                          lambda v: replace(BASE, z_hi=10.0 + v))
    out["placement_shift"] = study("sh", "both monitors shifted", [-2.0, 0.0, 2.0, 3.025],
                                   lambda v: replace(BASE, monitors=(0.0 + v, 10.0 + v), z_lo=-8.0 + min(0, v), z_hi=18.0 + max(0, v)))
    # 6. delay vs length
    Ls = [2.5, 5.0, 10.0, 20.0, 40.0]
    len_rows = study("L", "L", Ls, lambda v: replace(BASE, monitors=(0.0, v), z_hi=v + 8.0), L_of=lambda v: v)
    slope, icpt, resid = X.fit_delay_vs_length(Ls, [r["mean_tau"] for r in len_rows])
    out["length"] = dict(rows=len_rows, slope=slope, intercept=icpt, max_residual=resid, expected_slope=1.0)
    out["provenance"] = pio.provenance()
    pio.save_json(RES / "1a_vacuum.json", out, with_provenance=False)
    print("1A done; base mean_err", out["base"]["mean_err"], "res order", out["resolution"]["order_mean_err"],
          "length slope", slope)
