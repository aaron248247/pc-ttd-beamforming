"""Numerical convergence of the Phase 0.5 metrics for chosen designs.

usage: 05_convergence.py CONFIG [CONFIG ...]     (CONFIG in H, R, H_r030)
Cases per design: base (as swept) | res48 tol1e-8 | supercell M+1 | tol1e-9.  Same analysis code as the sweep.
"""
import multiprocessing as mp_
import sys
import numpy as np
from _common import *  # noqa
from pcttd import io as pio

KEYS = ["sm_strict_even_frac", "ng_max_strict_even", "gbp10", "gap_rel", "gap_w_lo", "gap_w_hi", "c6_pass"]


def case(args):
    name, lattice, r, n_gap, f1, s, tag, res, Mval, tol = args
    from pcttd import mpb_bands as MB, bands as B
    gaps, _ = MB.bulk_gaps(lattice, r, EPS, res, n_max=n_gap)
    g = [g for g in gaps if g.n == n_gap][0]
    d = MB.defect_bands(lattice, Mval, r, f1, s, EPS, res, 41, n_gap, 8, tol)
    out = MB.analyze_design(d, g)
    out.update(tag=tag, res=res, M=Mval, tol=tol, f1=f1, s=s, wall_s=d["wall_s"])
    np.savez_compressed(RES / "defect" / name / f"conv_{tag}_f1_{f1:.2f}_s_{s:+.2f}.npz",
                        **{k: v for k, v in d.items() if isinstance(v, np.ndarray)}, gap=np.array([g.w_lo, g.w_hi]))
    return dict(name=name, **out)


if __name__ == "__main__":
    sw = pio.load_json(RES / "defect_sweep.json")
    tasks = []
    for name in sys.argv[1:]:
        cfg = sw[name]
        ds = [d for d in cfg["designs"] if "skipped" not in d]
        designs = {(1.0, 0.0)}
        if ds:
            b = max(ds, key=lambda d: d["gbp10"])
            designs.add((b["f1"], b["s"]))
        for f1, s in sorted(designs):
            base = (name, cfg["lattice"], cfg["r_base"], cfg["n_gap"], f1, s)
            tasks += [base + ("res48_tol1e-8", 48, cfg["M"], 1e-8),
                      base + ("M+1", SWEEP_RES, cfg["M"] + 1, SWEEP_TOL),
                      base + ("tol1e-9", SWEEP_RES, cfg["M"], 1e-9)]
    with mp_.get_context("spawn").Pool(min(10, len(tasks))) as pool:
        res = pool.map(case, tasks, chunksize=1)
    pio.save_json(RES / ("convergence_" + "_".join(sys.argv[1:]) + ".json"),
                  dict(cases=res, base_from="defect_sweep.json (res 24, tol 1e-7, M as declared)"))
    print("done")
