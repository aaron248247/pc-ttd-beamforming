"""Phase 3B: run finite engineered-W1 devices (device section between plain-W1 leads) and save the raw (Ey, Hz) plane fields at the two reference planes.
Same method and defaults as Phase 2B (experiments/phase2/2b_fdtd.py); production decay_tol is 1e-7 and looser values are refused unless --exploratory=1.

usage: python 3b_fdtd.py <tag> <res> <f1> <s> <N> [<N> ...] [--key=value ...]    -> results/phase3/raw/<tag>_res<res>_N<N>.npz (raw, git-ignored)
Overrides (--d_src=6 --d_mon=5 --decay_tol=1e-7 --dpml_y=2.5 --J_pml=4 --dpml_x=5 --lead_pad=4 ...) change the W1Cfg defaults below; the
effective configuration is stored in each file (cfg_json).
Each run is a separate single-thread process (parallel across runs; no run depends on another, so concurrency cannot change a result).
"""
from _common import RES, RAW, np, ROOT
import json, sys, time
from pcttd import w1fdtd as W, w1modes as WM
from pcttd.io import provenance

P2 = ROOT / "results" / "phase2"            # plain-W1 (lead) MPB bands computed in Phase 2
FREQS = tuple(np.round(np.arange(0.220, 0.275 + 1e-9, 0.0005), 6))


def bands_for(J):
    """Hz-even bands of the plain W1 for supercell M = J (guided band index J + 2); M = 6 is the 2a result, others are cached on first use."""
    p = P2 / ("2a_mpb_plain.npz" if J == 6 else f"2a_mpb_M{J}.npz")
    if not p.exists():
        d = WM.even_bands(M=J, res=32, nk=101, nb=J + 3)
        np.savez(p, **{k: v for k, v in d.items() if isinstance(v, np.ndarray)})
    return dict(np.load(p))


def source_profile(f_src=0.2475, J=6):
    kap = WM.kappa_of_f(bands_for(J), J + 2, [f_src])
    P = WM.forward_profiles(kap, J + 2, M=J)
    hz = P["Hz"][0]
    hz = hz * np.exp(-1j * np.angle(hz[np.argmax(np.abs(hz))]))
    return P["y"], np.real(hz)


def one(job):
    tag, res, f1, s, N, over = job
    kw = dict(d_src=10, decay_tol=1e-7); kw.update(over)
    expl = kw.pop("exploratory", 0)
    if kw["decay_tol"] > 1e-7 and not expl:
        raise ValueError("production Phase 3 runs require decay_tol <= 1e-7 (Phase 2 lesson); pass --exploratory=1 to override")
    cfg = W.W1Cfg(N=N, res=res, f1=f1, s=s, freqs=FREQS, fcen=0.2475, fwidth=0.25, **kw)
    t = time.time()
    r = W.run_w1(cfg, source_profile(J=cfg.J))
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
    Ns = [int(x) for x in sys.argv[5:] if not x.startswith("--")]
    over = {}
    bands_for(int(dict(a[2:].split('=') for a in sys.argv[5:] if a.startswith('--')).get('J', 6)))   # cache MPB bands before the workers start
    for x in sys.argv[5:]:
        if x.startswith("--"):
            k, v = x[2:].split("=")
            over[k] = int(v) if v.lstrip("-").isdigit() else float(v)
    print(pool_map(one, [(tag, res, f1, s, N, over) for N in Ns], n=len(Ns)))
