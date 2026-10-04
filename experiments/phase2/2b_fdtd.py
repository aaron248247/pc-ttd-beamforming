"""Phase 2B (and reused for Phase 3/4): run finite W1 devices and save the raw (Ey, Hz) plane fields at the two reference planes.

usage: python 2b_fdtd.py <tag> <res> <f1> <s> <N> [<N> ...]    -> results/phase2/raw/<tag>_res<res>_N<N>.npz (raw, git-ignored)
Each run is a separate single-thread process (parallel across runs; no run depends on another, so concurrency cannot change a result).
"""
from _common import RES, RAW, np, ROOT
import json, sys, time
from pcttd import w1fdtd as W, w1modes as WM
from pcttd.io import provenance

FREQS = tuple(np.round(np.arange(0.220, 0.275 + 1e-9, 0.0005), 6))


def source_profile(f_src=0.2475):
    b = dict(np.load(RES / "2a_mpb_plain.npz"))
    kap = WM.kappa_of_f(b, 8, [f_src])
    P = WM.forward_profiles(kap, 8)
    hz = P["Hz"][0]
    hz = hz * np.exp(-1j * np.angle(hz[np.argmax(np.abs(hz))]))
    return P["y"], np.real(hz)


def one(job):
    tag, res, f1, s, N, outdir = job
    cfg = W.W1Cfg(N=N, res=res, f1=f1, s=s, d_src=10, freqs=FREQS, fcen=0.2475, fwidth=0.25, decay_tol=1e-6)
    t = time.time()
    r = W.run_w1(cfg, source_profile())
    r["wall_s"] = time.time() - t
    cfgd = {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.__dict__.items() if k != "freqs"}
    cfgd["freqs"] = [float(FREQS[0]), float(FREQS[-1]), len(FREQS)]
    r["cfg_json"] = json.dumps(cfgd)
    r["prov_json"] = json.dumps(provenance())
    np.savez_compressed(RAW / f"{tag}_res{res}_N{N}.npz", **r)
    return tag, res, N, r["wall_s"], r["t_end"], r["hit_max_time"]


if __name__ == "__main__":
    from _common import pool_map
    tag, res, f1, s = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    Ns = [int(x) for x in sys.argv[5:]]
    print(pool_map(one, [(tag, res, f1, s, N, None) for N in Ns], n=len(Ns)))
