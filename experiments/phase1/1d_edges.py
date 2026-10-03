"""Stop-band edges of the finite 2-D crystal from the saved N-section transmissions; compare with MPB. No simulation.

Definition (robust to Fabry-Perot fringes and to the narrow termination resonance inside the gap): with Tmax(f) = max over a set of
section lengths of T_N(f), the stop band is the contiguous interval around mid-gap where Tmax < thr. Using the maximum over several
lengths means a pass-band fringe that dips below thr for one length does not end the interval early. Edges are interpolated in log T.
Values below ~1e-7 are numerical floor (decay tolerance 1e-8) and never enter a decision: thr >= 1e-3 here.
"""
import numpy as np
from _common import *  # noqa
from pcttd import io as pio

SEED = 0.245


def interval(f, Tmax, thr):
    s = int(np.argmin(np.abs(f - SEED)))
    if Tmax[s] >= thr:
        return float("nan"), float("nan")
    lo = s
    while lo > 0 and Tmax[lo] < thr:
        lo -= 1
    hi = s
    while hi < len(f) - 1 and Tmax[hi] < thr:
        hi += 1
    def cross(a, b):                       # crossing of thr between samples a (below) and b (above), in log T
        la, lb, lt = np.log(max(Tmax[a], 1e-300)), np.log(Tmax[b]), np.log(thr)
        return float(f[a] + (lt - la) * (f[b] - f[a]) / (lb - la))
    return cross(lo + 1, lo), cross(hi - 1, hi)


if __name__ == "__main__":
    A = pio.load_json(RES / "1d_attenuation.json")
    M = pio.load_json(RES / "1d_mpb_directional.json")
    mpb128 = M["rows"][-1]["even_gap"][0][:2]
    mpb64 = M["rows"][1]["even_gap"][0][:2]
    out = dict(mpb_even_gap_res128=mpb128, mpb_even_gap_res64=mpb64, seed=SEED, rows=[], curves={})
    for ny in A["nys"]:
        for ns in ([8, 10, 12], [10, 12, 16], [8, 10, 12, 16], [6, 8, 10, 12, 16], [6], [8], [10], [12], [16]):
            runs = {r["N"]: r for r in A["runs"] if r["ny"] == ny and r["N"] in ns}
            f = np.array(runs[ns[0]]["f"])
            Tmax = np.max([np.array(runs[n]["T"]) for n in ns], axis=0)
            for thr in (1e-2, 1e-3):
                lo, hi = interval(f, Tmax, thr)
                out["rows"].append(dict(ny=ny, Ns=ns, thr=thr, edge_lo=lo, edge_hi=hi, d_lo=lo - mpb128[0], d_hi=hi - mpb128[1]))
        out["curves"][str(ny)] = {str(r["N"]): dict(f=r["f"], T=r["T"]) for r in A["runs"] if r["ny"] == ny}
    # variability bands
    fin = [r for r in out["rows"] if r["ny"] == max(A["nys"]) and np.isfinite(r["edge_lo"]) and np.isfinite(r["edge_hi"])]
    allr = [r for r in out["rows"] if np.isfinite(r["edge_lo"]) and np.isfinite(r["edge_hi"])]
    def span(rows, k):
        v = [r[k] for r in rows]
        return dict(min=min(v), max=max(v), mean=float(np.mean(v)))
    out["summary"] = dict(
        finest_ny=max(A["nys"]),
        finest=dict(lo=span(fin, "edge_lo"), hi=span(fin, "edge_hi"), d_lo=span(fin, "d_lo"), d_hi=span(fin, "d_hi"), n=len(fin)),
        all=dict(lo=span(allr, "edge_lo"), hi=span(allr, "edge_hi"), d_lo=span(allr, "d_lo"), d_hi=span(allr, "d_hi"), n=len(allr)))
    # deep-gap decay with N using only reliable values (T > 3e-7 for every length used)
    rel = []
    for ny in A["nys"]:
        runs = {r["N"]: r for r in A["runs"] if r["ny"] == ny}
        Ns_ = [6, 8, 10]
        Tm = np.array([runs[n]["T"] for n in Ns_]); f = np.array(runs[6]["f"])
        w = (f > 0.222) & (f < 0.268) & np.all(Tm > 3e-7, axis=0)
        lnT = np.array([np.mean(np.log(Tm[i][w])) for i in range(3)])
        rel.append(dict(ny=ny, n_freq=int(w.sum()), N=Ns_, mean_lnT=[float(x) for x in lnT],
                        decay_per_column=float(-np.polyfit(Ns_, lnT, 1)[0])))
    out["deep_gap_decay"] = rel
    # (1) discretization contribution: change of each single-length edge between the coarsest and finest ny
    nys_ = sorted(A["nys"])
    res_sp = []
    for N in (8, 10, 12, 16):
        a = [r for r in out["rows"] if r["ny"] == nys_[0] and r["Ns"] == [N] and r["thr"] == 1e-2][0]
        b = [r for r in out["rows"] if r["ny"] == nys_[-1] and r["Ns"] == [N] and r["thr"] == 1e-2][0]
        res_sp.append(dict(N=N, d_lo=b["edge_lo"] - a["edge_lo"], d_hi=b["edge_hi"] - a["edge_hi"]))
    out["resolution_change_coarsest_to_finest"] = dict(ny=[nys_[0], nys_[-1]], rows=res_sp,
                                                    max_abs=max(max(abs(r["d_lo"]), abs(r["d_hi"])) for r in res_sp))
    # (2) finite-length / threshold contribution at the finest resolution (single lengths N >= 8)
    fin_single = [r for r in out["rows"] if r["ny"] == nys_[-1] and len(r["Ns"]) == 1 and r["Ns"][0] >= 8 and r["thr"] == 1e-2]
    out["single_length_offsets_finest"] = dict(
        lower=[min(r["d_lo"] for r in fin_single), max(r["d_lo"] for r in fin_single)],
        upper=[min(r["d_hi"] for r in fin_single), max(r["d_hi"] for r in fin_single)], rows=[dict(N=r["Ns"][0], d_lo=r["d_lo"], d_hi=r["d_hi"]) for r in fin_single])
    # (3) symmetry test: T inside the interval where MPB has an odd-parity band but no even one (0.2744-0.289)
    anyg = M["rows"][-1]["any_gap"][0][1]
    w0, w1 = anyg + 0.001, mpb128[1] - 0.001
    odd = []
    for ny in A["nys"]:
        for r in [x for x in A["runs"] if x["ny"] == ny]:
            f_ = np.array(r["f"]); T_ = np.array(r["T"]); m_ = (f_ >= w0) & (f_ <= w1)
            odd.append(dict(ny=ny, N=r["N"], f_lo=float(f_[m_].min()), f_hi=float(f_[m_].max()), T_max=float(T_[m_].max()), T_geomean=float(np.exp(np.mean(np.log(np.maximum(T_[m_], 1e-9)))))))
    out["odd_band_window"] = dict(f_window=[w0, w1], rows=odd)
    # termination dependence of the edges (crystal faces moved by +-offset)
    trows = []
    for fn, nsel in (("1d_termination.json", (8, 10)), ("1d_termination_N12_16.json", (12, 16))):
        if not (RES / fn).exists():
            continue
        Tm_ = pio.load_json(RES / fn)
        for o in Tm_["offsets"]:
            runs = {r["N"]: r for r in Tm_["runs"] if abs(r["face_offset"] - o) < 1e-9}
            f = np.array(runs[nsel[0]]["f"])
            Tmax = np.max([np.array(runs[n]["T"]) for n in nsel], axis=0)
            g = (f >= .222) & (f <= .268)
            peak = f[g][int(np.argmax(Tmax[g]))]
            for thr in (1e-2, 1e-3):
                lo, hi = interval(f, Tmax, thr)
                trows.append(dict(Ns=list(nsel), face_offset=o, thr=thr, edge_lo=lo, edge_hi=hi, ingap_peak_f=float(peak),
                                  ingap_peak_T=float(Tmax[g].max()), ingap_Tmax_overall=float(Tmax[(f >= 0.215) & (f <= 0.283)].max()),
                                  Tgeo={str(n): float(np.exp(np.mean(np.log(np.maximum(np.array(runs[n]["T"])[g], 1e-9))))) for n in nsel}))
    out["termination"] = trows
    out["provenance"] = pio.provenance()
    pio.save_json(RES / "1d_edges.json", out, with_provenance=False)
    print("MPB even gap res128", mpb128, "res64", mpb64)
    for k in ("finest", "all"):
        s = out["summary"][k]
        print(k, "n=%d" % s["n"], "lo [%.5f, %.5f] hi [%.5f, %.5f]  d_lo [%.1e, %.1e] d_hi [%.1e, %.1e]" % (s["lo"]["min"], s["lo"]["max"], s["hi"]["min"], s["hi"]["max"], s["d_lo"]["min"], s["d_lo"]["max"], s["d_hi"]["min"], s["d_hi"]["max"]))
    for r in out["rows"]:
        if r["thr"] == 1e-2 and r["Ns"] == [8, 10, 12]:
            print("  ny %d thr %.0e: lo %.5f hi %.5f (d %.1e, %.1e)" % (r["ny"], r["thr"], r["edge_lo"], r["edge_hi"], r["d_lo"], r["d_hi"]))
    for r in out.get("termination", []):
        print("  termination Ns %s offset %+.2f thr %.0e: lo %.5f hi %.5f in-gap peak f=%.3f T=%.1e" % (r["Ns"], r["face_offset"], r["thr"], r["edge_lo"], r["edge_hi"], r["ingap_peak_f"], r["ingap_peak_T"]))
    for r in rel: print("  decay/col ny %d: %.3f over %d freqs" % (r["ny"], r["decay_per_column"], r["n_freq"]))
