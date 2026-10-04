"""Phase 2A: MPB dispersion of the plain r/a = 0.30 W1 (the lead and the Phase 2 device): guided-band identity, v_g, n_g, confinement,
resolution check, and the pre-declared useful interval.  Saves results/phase2/2a_mpb_plain.{npz,json}."""
from _common import RES, np
import os
os.environ["OMP_NUM_THREADS"] = "4"
from pcttd import w1modes as WM, mpb_bands as MB
from pcttd.io import save_json

BAND = 8                     # 1-based index of the guided band in the Hz-even class (checked below against Phase 0.5 parity machinery)
M = 6

hi = WM.even_bands(M=M, res=32, nk=101, nb=9)
lo = WM.even_bands(M=M, res=24, nk=51, nb=9)
vhi = WM.even_bands(M=M, res=48, nk=41, nb=9)
np.savez(RES / "2a_mpb_plain.npz", **{k: v for k, v in hi.items() if isinstance(v, np.ndarray)})

keep = hi["vg"][:, BAND - 1] < -2e-3
f = hi["freq"][keep, BAND - 1]; ng = 1.0 / np.abs(hi["vg"][keep, BAND - 1]); k = hi["k"][keep]
i0 = int(np.argmin(ng))
# pre-declared useful interval: contiguous interval around the n_g minimum with n_g <= 1.10 n_g,min
ok = ng <= 1.10 * ng[i0]
a = i0
while a > 0 and ok[a - 1]: a -= 1
b = i0
while b < len(ng) - 1 and ok[b + 1]: b += 1
interval = (float(f[b]), float(f[a]))           # f decreases with index along kappa

def ng_at(d, ff):
    kp = d["vg"][:, BAND - 1] < -2e-3
    return float(np.interp(ff, d["freq"][kp, BAND - 1][::-1], (1 / np.abs(d["vg"][kp, BAND - 1]))[::-1]))

fgrid = np.array([0.23, 0.24, 0.25, 0.26, 0.27])
conv_table = {str(d["res"]): [ng_at(d, x) for x in fgrid] for d in (lo, hi, vhi)}

# parity and confinement from the Phase 0.5 machinery (independent code path: parity by field overlap, energy fractions)
d05 = MB.defect_bands("triangular", M=M, r=0.30, res=24, nk=41, gap_n=1, extra=10)
kk = d05["k"]
sel = (d05["parity"] == +1)
conf = []
for ff in fgrid:
    # the Hz-even band whose frequency passes through ff: pick the band index (0-based over all bands) and k by interpolation
    best = None
    for bnd in range(d05["freq"].shape[1]):
        if sel[:, bnd].all():
            fr = d05["freq"][:, bnd]
            if fr.min() < ff < fr.max():
                j = int(np.argmin(np.abs(fr - ff)))
                best = (bnd, j, float(d05["strip"][j, bnd]), float(d05["edge"][j, bnd]), int(d05["parity"][j, bnd]))
    conf.append(dict(f=float(ff), match=best))
summary = dict(
    band_index_even_class=BAND, parity="Hz even under y -> -y (MPB run_te_yodd == Phase 0.5 parity +1)", M=M, r=0.30, eps=12.0,
    mpb_res=32, nk=101,
    ng_min=float(ng[i0]), f_at_ng_min=float(f[i0]), useful_interval_rule="contiguous n_g <= 1.10 n_g,min around the minimum",
    useful_interval=interval, useful_fractional_bandwidth=float((interval[1] - interval[0]) / (0.5 * (interval[0] + interval[1]))),
    ng_range_in_interval=[float(ng[a:b + 1].min()), float(ng[a:b + 1].max())],
    guided_band_frequency_range=[float(hi["freq"][:, BAND - 1].min()), float(hi["freq"][:, BAND - 1].max())],
    omnidirectional_gap_phase0p5_res24=[0.20721447488539896, 0.27455168679705266],
    ng_resolution_check=dict(freqs=fgrid, table=conv_table),
    confinement_phase0p5_machinery_res24=conf,
    edge_energy_max_all_bands=float(d05["edge"].max()), n_scatterer_rows=d05["n_scatterers"],
)
save_json(RES / "2a_mpb_plain.json", summary)
print({k: summary[k] for k in ("ng_min", "f_at_ng_min", "useful_interval", "useful_fractional_bandwidth", "ng_range_in_interval")})
print(conv_table)
for c in conf: print(c)
