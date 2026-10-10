"""Phase 3C: MPB resolution sensitivity of D1 and D2 (res 32 in 3a vs res 48 here), same supercell M = 6, nk = 101.
Saves results/phase3/3c_mpb_res48.{npz,json}; compared against 3a in 3d_mpb_analyze.py."""
from _common import RES, np
import os
os.environ["OMP_NUM_THREADS"] = "2"
from pcttd import w1modes as WM
from pcttd.io import save_json

DESIGNS = {"D1": (1.2, -0.05), "D2": (1.2, -0.10)}
arrays, out = {}, {}
for name, (f1, s) in DESIGNS.items():
    d = WM.even_bands(M=6, f1=f1, s=s, res=48, nk=101, nb=9)
    arrays[name + "_k"] = d["k"]; arrays[name + "_freq"] = d["freq"]; arrays[name + "_vg"] = d["vg"]
    out[name] = dict(f1=f1, s=s, res=48, nk=101, band_edge=float(d["freq"][-1, 7]), band_top=float(d["freq"][0, 7]))
np.savez(RES / "3c_mpb_res48.npz", **arrays)
save_json(RES / "3c_mpb_res48.json", out)
print(out)
