"""Phase 4A: run the four directly simulated D0 channels (N = N0 + m dN) and save the raw (Ey, Hz) plane fields at the two reference planes.
Same method and defaults as Phase 3B (experiments/phase3/3b_fdtd.py): d_src = 10, d_mon = 3, decay_tol = 1e-7 (looser values refused unless --exploratory=1).

usage: python 4a_fdtd.py <res> [N ...]   -> results/phase4/raw/D0_res<res>_N<N>.npz   (raw, git-ignored; default N = pre-declared lengths)
"""
from _common import RAW, ROOT, LENGTHS, np
import json, sys, time
from pcttd import w1fdtd as W, w1modes as WM
from pcttd.io import provenance

P2 = ROOT / "results" / "phase2"
FREQS = tuple(np.round(np.arange(0.220, 0.275 + 1e-9, 0.0005), 6))


def source_profile(f_src=0.2475, J=6):
    bands = dict(np.load(P2 / "2a_mpb_plain.npz"))
    kap = WM.kappa_of_f(bands, J + 2, [f_src])
    P = WM.forward_profiles(kap, J + 2, M=J)
    hz = P["Hz"][0]
    hz = hz * np.exp(-1j * np.angle(hz[np.argmax(np.abs(hz))]))
    return P["y"], np.real(hz)


def one(job):
    res, N, over = job
    kw = dict(d_src=10, decay_tol=1e-7); kw.update(over)
    if kw["decay_tol"] > 1e-7 and not kw.pop("exploratory", 0):
        raise ValueError("production runs require decay_tol <= 1e-7")
    kw.pop("exploratory", None)
    cfg = W.W1Cfg(N=N, res=res, f1=1.0, s=0.0, freqs=FREQS, fcen=0.2475, fwidth=0.25, **kw)
    t = time.time()
    r = W.run_w1(cfg, source_profile(J=cfg.J))
    r["wall_s"] = time.time() - t
    cfgd = {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.__dict__.items() if k != "freqs"}
    cfgd["freqs"] = [float(FREQS[0]), float(FREQS[-1]), len(FREQS)]
    r["cfg_json"] = json.dumps(cfgd)
    r["prov_json"] = json.dumps(provenance())
    np.savez_compressed(RAW / f"D0_res{res}_N{N}.npz", **r)
    return res, N, r["wall_s"], r["t_end"], r["hit_max_time"]


if __name__ == "__main__":
    from _common import pool_map
    res = int(sys.argv[1])
    Ns = [int(x) for x in sys.argv[2:] if not x.startswith("--")] or LENGTHS
    over = {k[2:].split("=")[0]: float(k.split("=")[1]) for k in sys.argv[2:] if k.startswith("--")}
    print(pool_map(one, [(res, N, over) for N in Ns], n=len(Ns)))
