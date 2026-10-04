"""Phase 2D: time-domain pulse check of the guided delay. Plain W1, a narrowband Gaussian pulse (source shaped as the MPB forward mode at fcen),
Hz(t) recorded at two on-axis points one integer number of periods apart (so the Bloch envelope is identical at both).

Observed:   dt_obs = centroid(|Hz(x2,t)|^2) - centroid(|Hz(x1,t)|^2)       (first moment of the energy trace)
Predicted:  dt_pred = L * sum_f |Hz_hat(x1,f)|^2 n_g(f) / sum_f |Hz_hat(x1,f)|^2    (spectral-power-weighted group delay, n_g from MPB)
(for a pulse the centroid delay equals the power-weighted group delay to first order; this removes the pulse-bandwidth bias from the comparison)
usage: python 2d_pulse.py <fcen> <fwidth> <res>   -> results/phase2/2d_pulse_f<fcen>_res<res>.json
"""
from _common import RES, np
import json, sys, time
import meep as mp
from scipy.interpolate import CubicSpline
from pcttd import w1fdtd as W, w1modes as WM, w1analysis as WA
from pcttd.io import save_json


def run(fcen, fwidth, res, L=20, x1=2.0, N=30):
    cfg = W.W1Cfg(N=N, res=res, d_src=10, fcen=fcen, fwidth=fwidth, d_mon=3)
    b = dict(np.load(RES / "2a_mpb_plain_fine.npz"))
    kap = WM.kappa_of_f(b, 8, [fcen]); P = WM.forward_profiles(kap, 8, M=6, res=32)
    hz = P["Hz"][0]; hz = np.real(hz * np.exp(-1j * np.angle(hz[np.argmax(np.abs(hz))])))
    sim = W.build(cfg, (P["y"], hz))
    # source here is a pulse that must stay narrow: build() uses GaussianSource(fcen, fwidth) from cfg
    pts = [mp.Vector3(x1, 0.0), mp.Vector3(x1 + L, 0.0)]
    tr = [[], []]
    def rec(s):
        for i, p in enumerate(pts):
            tr[i].append(s.get_field_point(mp.Hz, p).real)
    with W.quiet():
        sim.run(mp.at_every(0.5, rec), until_after_sources=(N + 10) * 8.0)
    return np.arange(len(tr[0])) * 0.5, np.array(tr[0]), np.array(tr[1]), cfg


if __name__ == "__main__":
    fcen, fwidth, res = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    t0 = time.time()
    t, a, bb, cfg = run(fcen, fwidth, res)
    dt_obs = WA.pulse_centroid_delay(t, a, bb)
    # spectral weights from the first trace (dt = 0.5; periodogram, Hann-free: the trace is a finite transient that has decayed)
    A = np.abs(np.fft.rfft(a)) ** 2; fr = np.fft.rfftfreq(len(a), 0.5)
    m = (fr > fcen - 3 * fwidth / 2) & (fr < fcen + 3 * fwidth / 2) & (fr > 0.2145) & (fr < 0.2805)
    b = dict(np.load(RES / "2a_mpb_plain_fine.npz")); keep = b["vg"][:, 7] < -2e-3
    ngf = CubicSpline(b["freq"][keep, 7][::-1], (1 / np.abs(b["vg"][keep, 7]))[::-1])
    ng_w = float(np.sum(A[m] * ngf(fr[m])) / np.sum(A[m]))
    out = dict(fcen=fcen, fwidth=fwidth, res=res, L=20, dt_obs=dt_obs, ng_obs=dt_obs / 20, ng_pred_mpb_weighted=ng_w, ng_mpb_at_fcen=float(ngf(fcen)),
               rel_diff=dt_obs / 20 / ng_w - 1, wall_s=time.time() - t0, t_end=float(t[-1]))
    save_json(RES / f"2d_pulse_f{fcen:.4f}_res{res}.json", out)
    print(out)
