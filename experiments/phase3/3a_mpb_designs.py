"""Phase 3A: MPB dispersion (Hz-even guided band 8) of the three selected H_r030 designs and five more sweep points (screening),
plus the conventional square-lattice R air-channel reference. Saves results/phase3/3a_mpb_designs.{npz,json}.

All geometries are points of the existing Phase 0.5 H_r030 sweep (f1 = radius factor, s = outward shift of the first rows); nothing is invented.
"""
from _common import RES, np
import os
os.environ["OMP_NUM_THREADS"] = "2"
from pcttd import w1modes as WM, mpb_bands as MB
from pcttd.io import save_json

M = 6
DESIGNS = {"D0": (1.0, 0.0), "D1": (1.2, -0.05), "D2": (1.2, -0.10)}
SCREEN = [(1.1, -0.1), (1.1, -0.05), (0.9, -0.1), (1.0, -0.1), (1.0, -0.05), (1.1, 0.0)]
BAND = 8


def branch(d, band=BAND):
    keep = d["vg"][:, band - 1] < -2e-3
    f = d["freq"][keep, band - 1]; ng = 1.0 / np.abs(d["vg"][keep, band - 1])
    o = np.argsort(f)
    return f[o], ng[o]


out, arrays = {}, {}
for name, (f1, s) in DESIGNS.items():
    d = WM.even_bands(M=M, f1=f1, s=s, res=32, nk=101, nb=9)
    f, ng = branch(d)
    arrays[name + "_k"] = d["k"]; arrays[name + "_freq"] = d["freq"]; arrays[name + "_vg"] = d["vg"]
    out[name] = dict(f1=f1, s=s, res=32, band=BAND, band_edge=float(d["freq"][-1, BAND - 1]), band_top=float(d["freq"][0, BAND - 1]))
for f1, s in SCREEN:
    d = WM.even_bands(M=M, f1=f1, s=s, res=24, nk=61, nb=9)
    out[f"screen_{f1:.2f}_{s:+.2f}"] = dict(f1=f1, s=s, res=24, band=BAND, band_edge=float(d["freq"][-1, BAND - 1]), band_top=float(d["freq"][0, BAND - 1]),
                                           ng_table={f"{x:.3f}": float(np.interp(x, *branch(d))) for x in np.arange(0.23, 0.281, 0.01)
                                                     if branch(d)[0].min() < x < branch(d)[0].max()})

# conventional reference: square-lattice TM air-channel waveguide (Phase 0.5 candidate R: r/a = 0.19 rods, eps = 12, one rod row removed)
dr = MB.defect_bands("square", M=6, r=0.19, res=24, nk=61, gap_n=1, extra=8)
k = dr["k"]
best = None
for b in range(dr["freq"].shape[1]):
    if np.all(dr["parity"][:, b] == +1):
        fr = dr["freq"][:, b]; vg = dr["vg"][:, b]
        inside = (fr > 0.30) & (fr < 0.42) & (vg > 2e-3)
        if inside.sum() >= 5 and (best is None or inside.sum() > best[1]):
            best = (b, int(inside.sum()))
b = best[0]
fr = dr["freq"][:, b]; vg = dr["vg"][:, b]
sel = (vg > 2e-3) & (fr > 0.30) & (fr < 0.42)
arrays["R_k"] = k[sel]; arrays["R_freq"] = fr[sel]; arrays["R_vg"] = vg[sel]
arrays["R_strip"] = dr["strip"][sel, b]; arrays["R_edge"] = dr["edge"][sel, b]
out["R"] = dict(lattice="square", pol="TM", r=0.19, M=6, res=24, band_index_0based=int(b), fr_range=[float(fr[sel].min()), float(fr[sel].max())],
                edge_energy_max=float(dr["edge"][sel, b].max()))
np.savez(RES / "3a_mpb_designs.npz", **arrays)
save_json(RES / "3a_mpb_designs.json", out)
print({k_: (v["band_edge"] if "band_edge" in v else v) for k_, v in out.items()})
