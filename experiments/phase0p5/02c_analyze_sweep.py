"""Re-run the analysis (bands.py / mpb_bands.analyze_design) on the saved raw band data and rewrite the sweep JSON.

Raw MPB output (results/phase0p5/defect/<config>/*.npz) is the source of truth; metrics in the sweep JSON are derived.
usage: 02c_analyze_sweep.py [defect_sweep.json | defect_sweep_supp.json]
"""
import sys
import numpy as np
from _common import *  # noqa
from pcttd import io as pio, bands as B
from pcttd import mpb_bands as MB

if __name__ == "__main__":
    fn = sys.argv[1] if len(sys.argv) > 1 else "defect_sweep.json"
    sw = pio.load_json(RES / fn)
    n_bad = 0
    for name, cfg in sw.items():
        if name == "provenance":
            continue
        gap = cfg["gap_at_sweep_res"]
        g = B.Gap(gap["n"], gap["w_lo"], gap["w_hi"])
        for rec in cfg["designs"]:
            if "skipped" in rec:
                continue
            z = np.load(RES / "defect" / name / f"f1_{rec['f1']:.2f}_s_{rec['s']:+.2f}.npz")
            d = dict(k=z["k"], freq=z["freq"], vg=z["vg"], parity=z["parity"], strip=z["strip"], edge=z["edge"])
            rec.update(MB.analyze_design(d, g))
            n_bad += (not rec["bands_sufficient"])
        cfg["n_bands_insufficient"] = sum(1 for r in cfg["designs"] if "skipped" not in r and not r["bands_sufficient"])
        print(name, "designs", len(cfg["designs"]), "bands_insufficient:", cfg["n_bands_insufficient"])
    pio.save_json(RES / fn, {k: v for k, v in sw.items() if k != "provenance"})
