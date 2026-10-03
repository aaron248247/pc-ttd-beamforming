"""Bulk band structures along the high-symmetry path at each baseline radius (for figures)."""
import numpy as np
from _common import *  # noqa
from pcttd import io as pio
from pcttd import mpb_bands as MB

if __name__ == "__main__":
    bulk = pio.load_json(RES / "bulk_gap_scan.json")
    cfg = {"H": ("triangular", bulk["triangular"]["r_base"]), "R": ("square", bulk["square"]["r_base"]),
           "H_r030": ("triangular", 0.30)}
    arrs = {}
    for name, (lat, r) in cfg.items():
        b = MB.bulk_bands(lat, r, EPS, BULK_RES, nbands=6, path_per_segment=16, grid_n=14)
        arrs[f"{name}_path"] = b["path"]
        arrs[f"{name}_grid"] = b["grid"]
        arrs[f"{name}_r"] = np.array(r)
    np.savez_compressed(RES / "bulk_bands.npz", **arrs)
    print("saved")
