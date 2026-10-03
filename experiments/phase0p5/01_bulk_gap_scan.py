"""Step 1 of the pre-declared rule: radius scan -> candidate gap and baseline radius, for both lattices."""
import multiprocessing as mp_
import numpy as np
from _common import *  # noqa
from pcttd import io as pio


def row(args):
    lattice, r = args
    from pcttd import mpb_bands as MB
    gaps, _ = MB.bulk_gaps(lattice, r, EPS, BULK_RES, n_max=3)
    return dict(lattice=lattice, r=r, gaps=[dict(n=g.n, w_lo=g.w_lo, w_hi=g.w_hi, rel=g.rel) for g in gaps])


def scan(pool, lattice, radii):
    return pool.map(row, [(lattice, float(r)) for r in radii])


def best_by_n(rows):
    best = {}
    for rw in rows:
        for g in rw["gaps"]:
            if g["n"] not in best or g["rel"] > best[g["n"]]["rel"]:
                best[g["n"]] = dict(rel=g["rel"], r=rw["r"])
    return best


if __name__ == "__main__":
    ranges = {"triangular": np.arange(0.20, 0.4001, 0.02), "square": np.arange(0.10, 0.3501, 0.02)}
    out = {}
    with mp_.get_context("spawn").Pool(10) as pool:
        for lat, radii in ranges.items():
            rows = scan(pool, lat, np.round(radii, 4))
            best = best_by_n(rows)
            n_gap = next((n for n in sorted(best) if best[n]["rel"] >= 0.05), None)
            refined = []
            if n_gap is not None:
                r0 = best[n_gap]["r"]
                lo_, hi_ = radii.min() - 1e-9, radii.max() + 1e-9      # stay inside the declared scan range
                fine = [r for r in np.round(np.arange(r0 - 0.02, r0 + 0.0201, 0.005), 4) if lo_ <= r <= hi_]
                done = {rw["r"] for rw in rows}
                refined = scan(pool, lat, [r for r in fine if float(r) not in done])
            allrows = sorted(rows + refined, key=lambda d: d["r"])
            best = best_by_n(allrows)
            if n_gap is not None:
                n_gap = next((n for n in sorted(best) if best[n]["rel"] >= 0.05), None)
            out[lat] = dict(rows=allrows, best_by_n=best, n_gap=n_gap,
                            r_base=best[n_gap]["r"] if n_gap else None,
                            rel_gap_base=best[n_gap]["rel"] if n_gap else None,
                            params=dict(eps=EPS, res=BULK_RES))
            print(lat, "n_gap", n_gap, "r_base", out[lat]["r_base"], "rel", out[lat]["rel_gap_base"], flush=True)
    pio.save_json(RES / "bulk_gap_scan.json", out)
