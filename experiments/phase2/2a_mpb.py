"""Phase 2A: MPB dispersion of the plain r/a = 0.30 W1 (the lead and the Phase 2 device).

Produces (results/phase2/2a_mpb_plain.{npz,json}):
  * guided-band identity: the Hz-even class from ``run_te_yodd`` (w1modes) is compared band by band, at every k, with the Phase 0.5
    machinery (mpb_bands.defect_bands: parity by field overlap, same supercell) -> max frequency difference, parity and strip-confinement of band 8;
  * v_g, n_g from Hellmann-Feynman velocities, and the pre-declared useful interval (n_g <= 1.10 n_g,min around the minimum);
  * derivative-grid convergence: n_g from a cubic-spline derivative of f(kappa) on nk = 26, 51, 101, 201 vs the Hellmann-Feynman n_g (exact derivative),
    reported over the useful interval and over the whole retained branch (where the derivative is ill conditioned near the band ends);
  * spatial-resolution check of n_g (MPB res 24, 32, 48).
"""
from _common import RES, np
import os
os.environ["OMP_NUM_THREADS"] = "4"
from scipy.interpolate import CubicSpline
from pcttd import w1modes as WM, mpb_bands as MB
from pcttd.io import save_json

BAND = 8                     # 1-based index of the guided band in the Hz-even class
M = 6
VG_FLOOR = 2e-3

fine = WM.even_bands(M=M, res=32, nk=201, nb=9)
hi = WM.even_bands(M=M, res=32, nk=101, nb=9)
lo = WM.even_bands(M=M, res=24, nk=51, nb=9)
vhi = WM.even_bands(M=M, res=48, nk=41, nb=9)
np.savez(RES / "2a_mpb_plain.npz", **{k: v for k, v in hi.items() if isinstance(v, np.ndarray)})
np.savez(RES / "2a_mpb_plain_fine.npz", **{k: v for k, v in fine.items() if isinstance(v, np.ndarray)})


def branch(d, band=BAND):
    keep = d["vg"][:, band - 1] < -VG_FLOOR
    return d["k"][keep], d["freq"][keep, band - 1], 1.0 / np.abs(d["vg"][keep, band - 1])


# ---- band identity vs Phase 0.5 machinery (same supercell, same k grid, same resolution)
ref = MB.defect_bands("triangular", M=M, r=0.30, res=32, nk=101, gap_n=1, extra=10)
nb = 9
dfreq = np.zeros((101, nb)); n_wrong_parity = 0
for ik in range(101):
    even = np.sort(ref["freq"][ik, ref["parity"][ik] == +1])[:nb]
    dfreq[ik] = np.abs(even - hi["freq"][ik, :nb])
ident = dict(max_freq_difference_all_9_even_bands=float(dfreq.max()), max_freq_difference_band8=float(dfreq[:, BAND - 1].max()))
# strip confinement of the guided band: the band of the Phase 0.5 solver that carries parity +1 and the band-8 frequency
conf_k, conf_strip, conf_edge, conf_f = [], [], [], []
for ik in range(101):
    idx = np.nonzero(ref["parity"][ik] == +1)[0]
    idx = idx[np.argsort(ref["freq"][ik, idx])][BAND - 1]
    conf_k.append(ref["k"][ik]); conf_f.append(ref["freq"][ik, idx]); conf_strip.append(ref["strip"][ik, idx]); conf_edge.append(ref["edge"][ik, idx])
conf_k, conf_f, conf_strip, conf_edge = map(np.array, (conf_k, conf_f, conf_strip, conf_edge))
dec = hi["vg"][:, BAND - 1] < -VG_FLOOR

# ---- useful interval (rule fixed in docs/phase2_5_criteria.md)
k, f, ng = branch(hi)
i0 = int(np.argmin(ng))
ok = ng <= 1.10 * ng[i0]
a = b = i0
while a > 0 and ok[a - 1]: a -= 1
while b < len(ng) - 1 and ok[b + 1]: b += 1
interval = (float(f[b]), float(f[a]))

# ---- derivative-grid convergence
kf, ff, ngf = branch(fine)                                   # reference: Hellmann-Feynman on the finest grid
f_lo, f_hi = max(ff.min(), 0.2145), min(ff.max(), 0.2805)
test_f = np.linspace(f_lo + 0.0005, f_hi - 0.0005, 400)
ref_ng = CubicSpline(ff[::-1], ngf[::-1])(test_f)


def spline_ng(nk):
    d = WM.even_bands(M=M, res=32, nk=nk, nb=9) if nk not in (101, 201) else (hi if nk == 101 else fine)
    kk = d["k"]; fr = d["freq"][:, BAND - 1]
    sel = d["vg"][:, BAND - 1] < -VG_FLOOR
    kk, fr = kk[sel], fr[sel]
    cs = CubicSpline(kk, fr)
    kap = CubicSpline(fr[::-1], kk[::-1])(test_f)
    return 1.0 / np.abs(cs(kap, 1))


in_int = (test_f >= interval[0]) & (test_f <= interval[1])
conv = {}
for nk in (26, 51, 101, 201):
    n = spline_ng(nk)
    rel = np.abs(n / ref_ng - 1)
    conv[str(nk)] = dict(max_rel_err_in_useful_interval=float(rel[in_int].max()), max_rel_err_whole_branch=float(rel.max()),
                         frac_of_branch_within_1pct=float(np.mean(rel < 0.01)))
# the Hellmann-Feynman n_g itself, nk = 101 vs 201, on the same frequencies
ng101 = CubicSpline(f[::-1], ng[::-1])(test_f)
hf_nk = dict(max_rel_err_in_useful_interval=float(np.abs(ng101 / ref_ng - 1)[in_int].max()),
             max_rel_err_whole_branch=float(np.abs(ng101 / ref_ng - 1).max()))

fgrid = np.array([0.23, 0.24, 0.25, 0.26, 0.27])


def ng_at(d, x):
    kk, fr, n = branch(d)
    return float(np.interp(x, fr[::-1], n[::-1]))


conv_res = {str(d["res"]): [ng_at(d, x) for x in fgrid] for d in (lo, hi, vhi)}

summary = dict(
    band_index_even_class=BAND, parity="Hz even under y -> -y (MPB run_te_yodd == Phase 0.5 parity +1)", M=M, r=0.30, eps=12.0,
    mpb_res=32, nk=101, vg_floor=VG_FLOOR,
    band_identity=ident,
    guided_band_confinement=dict(strip_energy_fraction_min=float(conf_strip[dec].min()), strip_energy_fraction_median=float(np.median(conf_strip[dec])),
                                 edge_energy_fraction_max=float(conf_edge[dec].max()), n_scatterer_rows=ref["n_scatterers"]),
    ng_min=float(ng[i0]), f_at_ng_min=float(f[i0]), useful_interval_rule="contiguous n_g <= 1.10 n_g,min around the minimum",
    useful_interval=interval, useful_fractional_bandwidth=float((interval[1] - interval[0]) / (0.5 * (interval[0] + interval[1]))),
    ng_range_in_interval=[float(ng[a:b + 1].min()), float(ng[a:b + 1].max())],
    guided_band_frequency_range=[float(hi["freq"][:, BAND - 1].min()), float(hi["freq"][:, BAND - 1].max())],
    omnidirectional_gap_phase0p5_res24=[0.20721447488539896, 0.27455168679705266],
    derivative_grid_convergence=dict(reference="Hellmann-Feynman n_g at nk = 201, res 32", test_freq_range=[float(test_f[0]), float(test_f[-1])],
                                     spline_derivative_of_band_vs_nk=conv, hellmann_feynman_nk101_vs_nk201=hf_nk),
    ng_resolution_check=dict(freqs=fgrid, table=conv_res),
)
save_json(RES / "2a_mpb_plain.json", summary)
print({k_: summary[k_] for k_ in ("band_identity", "guided_band_confinement", "ng_min", "useful_interval", "useful_fractional_bandwidth")})
print(conv); print(hf_nk); print(conv_res)
