"""Phase 2E: convergence summary. Compares each variant's n_g (cascade pair given by --pair, default 6 12) with the base res-30 run and with MPB,
at representative frequencies (window-averaged over +-0.0025 to suppress derivative noise) and as rms/max over the useful interval.
usage: python 2e_convergence.py   -> results/phase2/2e_convergence.json
"""
from _common import RES, np
import json
from pcttd.io import save_json

REP = {"ordinary 0.2564 (n_g minimum)": 0.2564, "moderate slow 0.235": 0.235, "difficult 0.2245 (n_g rising to band edge)": 0.2245, "upper 0.270": 0.270}
VARIANTS = [("base", 20), ("base", 30), ("base", 40), ("dsrc6", 30), ("dmon5", 30), ("pml", 30), ("tol7", 30), ("tol8", 30), ("tol7", 40)]
PAIR = "ng_6_12"
ints = json.loads((RES / "2a_mpb_plain.json").read_text())["useful_interval"]


def win(f, y, f0, w=0.0025):
    m = np.abs(f - f0) <= w
    return float(np.mean(y[m]))


base = np.load(RES / "2c_base_res30.npz")
out = {}
for tag, res in VARIANTS:
    p = RES / f"2c_{tag}_res{res}.npz"
    if not p.exists():
        continue
    z = np.load(p); f = z["f"]
    m_int = (f >= ints[0]) & (f <= ints[1])
    d = dict(rep={}, vs_base_rms_useful=float(np.sqrt(np.mean((z[PAIR][m_int] / base[PAIR][m_int] - 1) ** 2))),
             vs_mpb_rms_useful=float(np.sqrt(np.mean((z[PAIR][m_int] / z["ng_mpb"][m_int] - 1) ** 2))),
             vs_mpb_median_useful=float(np.median(z[PAIR][m_int] / z["ng_mpb"][m_int] - 1)))
    for name, f0 in REP.items():
        d["rep"][name] = dict(ng=win(f, z[PAIR], f0), ng_mpb=win(f, z["ng_mpb"], f0), rel_vs_mpb=win(f, z[PAIR], f0) / win(f, z["ng_mpb"], f0) - 1)
    out[f"{tag}_res{res}"] = d
    print(f"{tag:6s} res{res}: rms vs MPB {d['vs_mpb_rms_useful']:.4f}  median {d['vs_mpb_median_useful']:+.4f}  rms vs base30 {d['vs_base_rms_useful']:.4f}  ",
          {k.split()[0] + k.split()[1][:6]: round(v["rel_vs_mpb"], 4) for k, v in d["rep"].items()})
save_json(RES / "2e_convergence.json", dict(pair=PAIR, window=0.0025, useful_interval=ints, variants=out))
