"""Step 2: 25-design defect sweep per candidate (first-neighbour-row radius scale f1 x outward shift s)."""
import multiprocessing as mp_
import time
import numpy as np
from _common import *  # noqa
from pcttd import io as pio, geometry as G, bands as B


def gap_task(args):
    name, lattice, r, n_gap = args
    from pcttd import mpb_bands as MB
    gaps, _ = MB.bulk_gaps(lattice, r, EPS, SWEEP_RES, n_max=max(n_gap, 1))
    g = [g for g in gaps if g.n == n_gap]
    if not g:
        raise RuntimeError(f"gap n={n_gap} vanished at res {SWEEP_RES}")
    return name, dict(n=g[0].n, w_lo=g[0].w_lo, w_hi=g[0].w_hi)


def design_task(args):
    name, lattice, r, n_gap, gap, f1, s, Mval = args
    from pcttd import mpb_bands as MB
    H, sc = (G.triangular_w1 if lattice == "triangular" else G.square_w1)(Mval, r, f1, s)
    mg = G.min_gap(sc, H)
    rec = dict(config=name, lattice=lattice, f1=f1, s=s, min_edge_gap=mg, unresolved=bool(mg < 3.0 / SWEEP_RES))
    if rec["unresolved"]:
        rec["skipped"] = "unresolved at screening resolution (min wall < 3 px) or overlapping"
        return rec
    t0 = time.time()
    d = MB.defect_bands(lattice, Mval, r, f1, s, EPS, SWEEP_RES, 41, n_gap, 8, SWEEP_TOL)
    g = B.Gap(gap["n"], gap["w_lo"], gap["w_hi"])
    rec.update(MB.analyze_design(d, g))
    rec.update(wall_s=d["wall_s"], H=d["H"], r_min=d["r_min"], eps_mirror_asym=d["eps_mirror_asym"])
    out = RES / "defect" / name
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / f"f1_{f1:.2f}_s_{s:+.2f}.npz", **{k: v for k, v in d.items()
                        if isinstance(v, np.ndarray)}, gap=np.array([gap["w_lo"], gap["w_hi"]]))
    return rec


if __name__ == "__main__":
    import sys
    supp = "supp" in sys.argv[1:]
    bulk = pio.load_json(RES / "bulk_gap_scan.json")
    if not supp:
        configs = {  # name -> (lattice, baseline radius, gap index, supercell M)
            "H": ("triangular", bulk["triangular"]["r_base"], bulk["triangular"]["n_gap"], M["triangular"]),
            "R": ("square", bulk["square"]["r_base"], bulk["square"]["n_gap"], M["square"]),
            "H_r030": ("triangular", 0.30, bulk["triangular"]["n_gap"], M["triangular"]),     # supplementary
        }
        outfile = "defect_sweep.json"
    else:   # supplementary supercell-width sensitivity (G3 edge-energy convergence)
        configs = {"R_M6": ("square", bulk["square"]["r_base"], bulk["square"]["n_gap"], 6)}
        outfile = "defect_sweep_supp.json"
    with mp_.get_context("spawn").Pool(10) as pool:
        gaps = dict(pool.map(gap_task, [(n, c[0], c[1], c[2]) for n, c in configs.items()]))
        tasks = [(n, c[0], c[1], c[2], gaps[n], f1, s, c[3]) for n, c in configs.items() for f1 in F1 for s in S]
        recs = pool.map(design_task, tasks, chunksize=1)
    out = {}
    for n, c in configs.items():
        out[n] = dict(lattice=c[0], r_base=c[1], n_gap=c[2], gap_at_sweep_res=gaps[n], M=c[3],
                      res=SWEEP_RES, tolerance=SWEEP_TOL, eps=EPS, decision_bearing=(n in ("H", "R")),
                      designs=[r for r in recs if r["config"] == n])
    pio.save_json(RES / outfile, out)
    print("done", {n: len(out[n]["designs"]) for n in out})
