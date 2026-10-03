"""Phase 1D (MPB side): directional bands of the triangular air-hole crystal along x at k_y = 0 (rectangular two-row cell).

A y-uniform Hz plane wave at normal incidence can only excite y-even (Hz) modes at k_y = 0, so the finite-crystal stop band
is compared with the frequency intervals that contain NO even-Phi band for any k_x in [0, 0.5] (all-parity gap shown as well).
"""
import numpy as np
from _common import *  # noqa
from pcttd import bands as B, geometry as G, io as pio, mpb_bands as MB

R_BASE = 0.30
EPS = 12.0


def directional(res, nk=81, extra=9, tol=1e-9, r=R_BASE):
    H = G.SQRT3
    sc = [G.Scatterer(0.0, 0.0, r), G.Scatterer(0.5, H / 2, r)]
    d = MB.defect_bands("triangular", 2, r, 1.0, 0.0, EPS, res, nk, 1, extra, tol, _override=(H, sc))
    levels = np.linspace(0.05, 0.45, 8001)
    br = B.class_branches(d["freq"], d["vg"], d["parity"], d["strip"], d["edge"])
    ok = B.branches_sufficient(br, levels[-1])
    c = B.find_crossings(d["k"], br, levels)
    n_even = c.counts(c.parity == +1)
    n_odd = c.counts(c.parity == -1)
    def zero_intervals(cnt):
        out = []
        for a, b in B.runs(cnt == 0):
            out.append((float(levels[a]), float(levels[b - 1]), float(levels[b - 1] - levels[a])))
        return out
    return d, dict(res=res, nk=nk, tol=tol, r=r, sufficient=bool(ok), top_band_min=float(d["freq"][:, -1].min()),
                   even_gap=zero_intervals(n_even), any_gap=zero_intervals(n_even + n_odd),
                   odd_gap=zero_intervals(n_odd), wall_s=d["wall_s"])


if __name__ == "__main__":
    out = dict(description="k_y = 0 directional gaps (rectangular cell 1 x sqrt(3), 2 holes), r/a = 0.30, eps = 12, TE (Hz)")
    rows = []
    for res in (32, 64, 128):
        d, info = directional(res)
        rows.append(info)
        if res == 64:
            np.savez_compressed(RES / "1d_mpb_directional_bands.npz", k=d["k"], freq=d["freq"], parity=d["parity"], vg=d["vg"])
        print(res, "sufficient", info["sufficient"], "even-gap:", [(round(a, 5), round(b, 5)) for a, b, _ in info["even_gap"]],
              "any-gap:", [(round(a, 5), round(b, 5)) for a, b, _ in info["any_gap"]], "wall", round(info["wall_s"], 1))
    out["rows"] = rows
    # bulk omnidirectional gap at the same radius for reference (full-BZ sample), from the Phase 0.5 machinery
    gaps, _ = MB.bulk_gaps("triangular", R_BASE, EPS, 64, n_max=2)
    out["bulk_omnidirectional_gap"] = [dict(n=g.n, w_lo=g.w_lo, w_hi=g.w_hi, rel=g.rel) for g in gaps]
    out["provenance"] = pio.provenance()
    pio.save_json(RES / "1d_mpb_directional.json", out, with_provenance=False)
