"""Phase 1D (FDTD side): finite triangular air-hole crystal, normal incidence, vs the MPB directional band structure."""
from dataclasses import replace
import numpy as np
from _common import *  # noqa
from pcttd import fdtd2d as F2, io as pio

BASE = F2.Cfg2D(N=8, ny=40)
SEED = 0.245                       # a frequency inside the predicted stop band
GAP_BAND = (0.22, 0.27)            # deep-gap window for attenuation statistics


def crossing(f, y, level, seed_idx, direction, log=False):
    """First crossing of y = level going outward from the seed index (direction -1 left, +1 right), linear interpolation."""
    i = seed_idx
    while 0 < i < len(f) - 1 and y[i] < level:
        i += direction
    j = i - direction
    ya, yb = (np.log10(y[j]), np.log10(y[i])) if log else (y[j], y[i])
    lv = np.log10(level) if log else level
    if yb == ya:
        return float(f[i])
    return float(f[j] + (lv - ya) * (f[i] - f[j]) / (yb - ya))


def metrics(r):
    f, T, R = r["f"], r["T"], r["R"]
    s = int(np.argmin(np.abs(f - SEED)))
    m_gap = (f >= GAP_BAND[0]) & (f <= GAP_BAND[1])
    return dict(
        edge_lo_3dB=crossing(f, T, 0.5, s, -1), edge_hi_3dB=crossing(f, T, 0.5, s, +1),
        edge_lo_20dB=crossing(f, T, 1e-2, s, -1, True), edge_hi_20dB=crossing(f, T, 1e-2, s, +1, True),
        T_gap_geomean=float(np.exp(np.mean(np.log(T[m_gap])))), T_gap_max=float(T[m_gap].max()),
        max_abs_RplusT_minus_1_gap=float(np.max(np.abs(T + R - 1)[m_gap])),
        max_abs_RplusT_minus_1_all=float(np.max(np.abs(T + R - 1))),
        t_end=r["t_end"], hit_max_time=r["hit_max_time"], dt=r["dt"])


def study(param, values, make):
    cfgs = [make(v) for v in values]
    rs = pool_map(F2.run_2d, cfgs)
    rows = []
    for v, c, r in zip(values, cfgs, rs):
        row = dict(param=param, value=v, N=c.N, ny=c.ny, res=c.res, **metrics(r))
        row["curve"] = dict(f=r["f"], T=r["T"], R=r["R"])
        rows.append(row)
    return rows


if __name__ == "__main__":
    mpb = pio.load_json(RES / "1d_mpb_directional.json")
    out = dict(config=cfg_dict(BASE), seed=SEED, gap_band=GAP_BAND,
               mpb_even_gap_res128=[mpb["rows"][-1]["even_gap"][0][0], mpb["rows"][-1]["even_gap"][0][1]],
               mpb_even_gap_res64=[mpb["rows"][1]["even_gap"][0][0], mpb["rows"][1]["even_gap"][0][1]],
               mpb_any_gap_res128=[mpb["rows"][-1]["any_gap"][0][0], mpb["rows"][-1]["any_gap"][0][1]])
    out["base"] = study("base", [0], lambda v: BASE)[0]
    out["length"] = study("N", [2, 3, 4, 6, 8, 12, 16], lambda v: replace(BASE, N=v))
    out["resolution"] = study("ny", [24, 32, 40, 48, 64], lambda v: replace(BASE, ny=v))
    out["pml"] = study("dpml", [1.0, 2.0, 3.0], lambda v: replace(BASE, dpml=v))
    out["decay"] = study("decay_tol", [1e-6, 1e-8, 1e-10], lambda v: replace(BASE, decay_tol=v))
    out["source_width"] = study("src_fwidth", [0.2, 0.3, 0.4, 0.6], lambda v: replace(BASE, src_fwidth=v))
    out["placement_monitor"] = study("mon_dist", [0.5, 1.0, 2.0], lambda v: replace(BASE, mon_dist=v))
    out["placement_source"] = study("src_dist", [1.0, 2.0, 3.0, 5.0], lambda v: replace(BASE, src_dist=v))
    out["provenance"] = pio.provenance()
    pio.save_json(RES / "1d_pc2d.json", out, with_provenance=False)
    b = out["base"]
    print("1D done: base edges 3dB %.4f %.4f  20dB %.4f %.4f  T_gap_geomean %.2e" % (b["edge_lo_3dB"], b["edge_hi_3dB"], b["edge_lo_20dB"], b["edge_hi_20dB"], b["T_gap_geomean"]))
