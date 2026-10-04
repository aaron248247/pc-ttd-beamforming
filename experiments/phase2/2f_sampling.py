"""Phase 2F: frequency-sampling and derivative-estimator sensitivity of n_g, from the converged 3-length phase-slope (tol 1e-7 res 40, and res 30).
The fine grid (df = 5e-4) is subsampled to 1e-3, 2e-3, 4e-3; n_g is differentiated with the central difference and with the local quadratic (+-3 samples; +-1 for the
coarsest grids).  Reported: rms/max over the useful interval of n_g relative to MPB, and the change relative to the df = 5e-4 central estimate
(pre-declared criterion: change <= 1% from 5e-4 to 2e-3).  usage: python 2f_sampling.py -> results/phase2/2f_sampling.json
"""
from _common import RES, np
import json
from pcttd import w1analysis as WA, guided as GD
from pcttd.io import save_json

ints = json.loads((RES / "2a_mpb_plain.json").read_text())["useful_interval"]
out = {}
for tag, res in (("tol7", 40), ("tol7", 30)):
    z = np.load(RES / f"2c_{tag}_res{res}.npz")
    f, sl, ng_mpb = z["f"], z["fit_slope"], z["ng_mpb"]
    ref = WA.ng_from_theta(f, sl, "central")
    res_ = {}
    for step in (1, 2, 4, 8):
        fs, (ss, nm) = WA.subsample(f, [sl, ng_mpb], step)
        m = (fs >= ints[0]) & (fs <= ints[1])
        for method, half in (("central", 0), ("poly", 3 if step <= 2 else 1)):
            n = WA.ng_from_theta(fs, ss, method, half=half)
            ok = m & np.isfinite(n)                       # the local quadratic is undefined at the end samples of a coarse grid
            rel_mpb = n[ok] / nm[ok] - 1
            chg = n[ok] / np.interp(fs[ok], f, ref) - 1
            res_[f"df={fs[1] - fs[0]:.4g}_{method}"] = dict(rms_vs_mpb=float(np.sqrt(np.mean(rel_mpb ** 2))), max_vs_mpb=float(np.abs(rel_mpb).max()),
                                                           median_vs_mpb=float(np.median(rel_mpb)),
                                                           rms_change_vs_ref=float(np.sqrt(np.mean(chg ** 2))), max_change_vs_ref=float(np.abs(chg).max()))
    out[f"{tag}_res{res}"] = res_
    print(tag, res)
    for k, v in res_.items():
        print(f"  {k:22s} rms vs MPB {v['rms_vs_mpb']:.4f} max {v['max_vs_mpb']:.4f} | change vs ref rms {v['rms_change_vs_ref']:.4f} max {v['max_change_vs_ref']:.4f}")
save_json(RES / "2f_sampling.json", dict(useful_interval=ints, results=out))
