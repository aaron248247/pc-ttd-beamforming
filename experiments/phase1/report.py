"""Generate docs/phase1_report_tables.md from results/phase1/*.json. No number is typed by hand."""
import json
import numpy as np
from _common import *  # noqa
from pcttd import io as pio

OUT = ROOT / "docs" / "phase1_report_tables.md"


def g(x, n=3):
    if x is None:
        return "—"
    if isinstance(x, float):
        return f"{x:.{n}g}"
    return str(x)


def table(head, rows):
    esc = lambda c: str(c).replace("|", "\\|")
    head = [esc(h) for h in head]
    s = "| " + " | ".join(head) + " |\n|" + "|".join("---" for _ in head) + "|\n"
    return s + "\n".join("| " + " | ".join(map(str, r)) + " |" for r in rows) + "\n"


def sweep_rows_1a(rows, keys):
    return [[g(r["value"])] + [g(r[k]) for k in keys] for r in rows]


if __name__ == "__main__":
    o = []
    A = pio.load_json(RES / "1a_vacuum.json")
    B = pio.load_json(RES / "1b_dielectric.json")
    C = pio.load_json(RES / "1c_bragg.json")
    o.append("## T1A. Vacuum propagation (L = 10, band %s)\n" % (A["band"],))
    o.append("Sign test (mean delay over band): sign +1 -> %s, sign -1 -> %s; expected +L/c = %s; established sign = **%s**.\n"
             % (g(A["sign_test"]["mean_tau_sign_plus"], 6), g(A["sign_test"]["mean_tau_sign_minus"], 6), g(A["sign_test"]["expected_L"]),
                A["sign_test"]["established_sign"]))
    o.append("\nBase (res 40): mean delay error %s, max |error| %s, max |measured - Yee prediction| %s, |H| deviation from 1: %s.\n"
             % (g(A["base"]["mean_err"], 4), g(A["base"]["max_abs_err"], 4), g(A["base"]["max_abs_resid_vs_yee"], 3), g(A["base"]["abs_H_dev"], 3)))
    o.append("\n**Resolution** (observed order: mean-error %s, max-error %s)\n\n" % (g(A["resolution"]["order_mean_err"], 4), g(A["resolution"]["order_max_err"], 4)))
    o.append(table(["res", "mean err", "max |err|", "max |meas - Yee|", "dt"], [[r["value"], g(r["mean_err"]), g(r["max_abs_err"]), g(r["max_abs_resid_vs_yee"]), g(r["dt"])] for r in A["resolution"]["rows"]]))
    o.append("\n**Delay vs length**: fit slope %s (expected 1 + Yee dispersion), intercept %s, max residual %s\n\n" % (g(A["length"]["slope"], 6), g(A["length"]["intercept"], 3), g(A["length"]["max_residual"], 3)))
    o.append(table(["L", "mean delay", "mean err"], [[g(r["value"]), g(r["mean_tau"], 8), g(r["mean_err"])] for r in A["length"]["rows"]]))
    for key, title, pname in (("pml", "PML thickness", "dpml"), ("decay", "Decay tolerance", "decay_tol"), ("source_width", "Source fwidth (spectral sigma = fwidth/2pi)", "fwidth"),
                              ("decimation", "DFT decimation (0 = automatic)", "decim"), ("placement_source_distance", "Source-to-monitor distance", "d"),
                              ("placement_pml_distance", "Monitor-to-PML distance", "d"), ("placement_shift", "Both monitors shifted (L fixed)", "shift")):
        o.append(f"\n**{title}**\n\n")
        o.append(table([pname, "mean err", "max |meas - Yee|", "|H| dev", "coverage", "t_end"],
                       [[g(r["value"]), g(r["mean_err"]), g(r["max_abs_resid_vs_yee"]), g(r["abs_H_dev"]), g(r["coverage"]), g(r["t_end"], 4)] for r in A[key]]))

    o.append("\n## T1B. Homogeneous dielectric and single slab vs TMM\n")
    o.append(table(["n", "slope", "slope / n", "intercept", "max residual"], [[g(d["n"], 5), g(d["slope"], 6), g(d["slope_over_n"], 6), g(d["intercept"]), g(d["max_residual"])] for d in B["delay_vs_length"]]))
    o.append("\nResolution, n = sqrt(12), L = 10: observed order (mean error) %s, (max error) %s\n\n" % (g(B["resolution_n_sqrt12"]["order_mean_err"], 4), g(B["resolution_n_sqrt12"]["order_max_err"], 4)))
    o.append(table(["res", "mean err", "max |err|", "max |meas - Yee|"], [[r["res"], g(r["mean_err"]), g(r["max_abs_err"]), g(r["max_abs_resid_vs_yee"])] for r in B["resolution_n_sqrt12"]["rows"]]))
    o.append("\nSlab (thickness 5, monitors at 0 and 10 in vacuum) vs TMM; delay compared with the SAME estimator applied to TMM data\n\n")
    o.append(table(["n", "res", "max |dt|", "max |dT|", "max |dR|", "max |R+T-1|", "max |excess-delay diff|", "estimator truncation", "max TMM excess delay"],
                   [[g(r["n"], 4), r["res"], g(r["max_abs_dt"]), g(r["max_abs_dT"]), g(r["max_abs_dR"]), g(r["max_abs_R_plus_T_minus_1"]), g(r["max_abs_excess_delay_err"]), g(r["estimator_truncation"]), g(r["max_excess_delay_tmm"])] for r in B["slab_vs_tmm"]]))
    for k in B:
        if k.startswith("slab_order"):
            o.append(f"\n{k}: observed order |dt| = {g(B[k]['order_dt'], 3)}, excess delay = {g(B[k]['order_delay'], 3)}\n")
    o.append("\nPML (n = 2) and decay tolerance (n = 2)\n\n")
    o.append(table(["dpml", "mean err", "max |meas - Yee|"], [[g(r["dpml"]), g(r["mean_err"]), g(r["max_abs_resid_vs_yee"])] for r in B["pml_n2"]]))
    o.append("\n")
    o.append(table(["decay tol", "mean err", "max |meas - Yee|", "|H| dev"], [[g(r["decay_tol"]), g(r["mean_err"]), g(r["max_abs_resid_vs_yee"]), g(r["abs_H_dev"])] for r in B["decay_n2"]]))

    o.append("\n## T1C. 1-D Bragg stack (n = 1 / sqrt(12) quarter-wave, N periods) vs TMM\n")
    o.append("Analytic infinite-lattice gap edges %s; quarter-wave frequency f0 = %s; TMM max delay (N=12) %s; frequency spacing %s.\n\n"
             % (tuple(round(x, 5) for x in C["analytic_gap_edges_infinite"]), g(C["f0"], 5), g(C["tmm_max_delay_N12"]), g(C["unwrap_guard_df"])))
    o.append("**Resolution** (N = 8). Observed orders: " + ", ".join(f"{k}: {g(v, 3)}" for k, v in C["resolution_orders"].items()) + "\n\n")
    o.append(table(["res", "s-1", "rms dT", "rms dT (after s)", "edge err lo", "edge err hi", "median rel err t", "median |delay err|", "max |delay err|", "max |R+T-1|"],
                   [[r["value"], g(r["freq_scale_minus1"]), g(r["rms_dT_unscaled"]), g(r["rms_dT_after_scale"]), g(r["edge_err_lo"]), g(r["edge_err_hi"]),
                     g(r["median_rel_complex_t"]), g(r["median_abs_delay_err"]), g(r["max_abs_delay_err"]), g(r["max_abs_R_plus_T_minus_1"])] for r in C["resolution"]]))
    o.append("\n**Number of periods** (res 80; -3 dB stop-band edges, MEEP vs TMM)\n\n")
    o.append(table(["N", "edge lo TMM", "edge lo MEEP", "edge hi TMM", "edge hi MEEP", "rms dT", "max |R+T-1|"],
                   [[r["N"], g(r["edge_lo_tmm"], 5), g(r["edge_lo_meep"], 5), g(r["edge_hi_tmm"], 5), g(r["edge_hi_meep"], 5), g(r["rms_dT_unscaled"]), g(r["max_abs_R_plus_T_minus_1"])] for r in C["periods"]]))
    for key, title in (("pml", "PML thickness"), ("decay", "Decay tolerance"), ("source_width", "Source fwidth"), ("placement_left", "Left monitor distance"),
                       ("placement_right", "Right monitor distance"), ("placement_source", "Source distance")):
        o.append(f"\n**{title}** (self-comparison with the last, most-converged setting; and vs TMM)\n\n")
        o.append(table(["value", "self max |dT|", "self max |dR|", "self max |delay|", "self max rel t", "vs-TMM rms dT", "max |R+T-1|", "source coverage"],
                       [[g(r["value"]), g(r["self_max_abs_dT"]), g(r["self_max_abs_dR"]), g(r["self_max_abs_delay"]), g(r["self_max_rel_complex_t"]), g(r["rms_dT_unscaled"]),
                         g(r["max_abs_R_plus_T_minus_1"]), g(r["source_coverage_min_over_max"])] for r in C[key]]))
    o.append("\n**Bloch decay in the stop band** (slope of ln T vs N, N = %s)\n\n" % (C["bloch_decay"]["N"],))
    o.append(table(["f", "MEEP decay / period", "analytic Bloch decay / period", "relative diff"],
                   [[g(r["f"]), g(r["meep_decay_per_period"], 5), g(r["analytic_bloch_decay"], 5), g((r["meep_decay_per_period"] - r["analytic_bloch_decay"]) / r["analytic_bloch_decay"], 3)] for r in C["bloch_decay"]["rows"]]))

    p2 = RES / "1d_pc2d.json"
    if p2.exists() and (RES / "1d_edges.json").exists():
        D = pio.load_json(p2)
        M = pio.load_json(RES / "1d_mpb_directional.json")
        Ed = pio.load_json(RES / "1d_edges.json")
        At = pio.load_json(RES / "1d_attenuation.json")
        o.append("\n## T1D. Finite 2-D triangular air-hole crystal (r/a = 0.30, eps = 12, TE) vs MPB\n")
        o.append("MPB directional (k_y = 0) even-Phi gap at res 32/64/128: " + "; ".join(f"[{g(r['even_gap'][0][0], 5)}, {g(r['even_gap'][0][1], 5)}]" for r in M["rows"]) +
                 ". All-parity gap, res 128: [%s, %s]. Bulk omnidirectional gap (n=1, res 64): [%s, %s].\n\n"
                 % (g(M["rows"][-1]["any_gap"][0][0], 5), g(M["rows"][-1]["any_gap"][0][1], 5), g(M["bulk_omnidirectional_gap"][0]["w_lo"], 5), g(M["bulk_omnidirectional_gap"][0]["w_hi"], 5)))
        o.append("**Stop band from Tmax(f) = max over section lengths N of T_N(f), contiguous interval around mid-gap with Tmax < thr** (long runs, minimum run time 20000, decay tolerance 1e-8)\n\n")
        o.append(table(["ny", "N set", "thr", "edge lo", "edge hi", "lo - MPB", "hi - MPB"],
                       [[r["ny"], r["Ns"], g(r["thr"]), g(r["edge_lo"], 5), g(r["edge_hi"], 5), g(r["d_lo"], 3), g(r["d_hi"], 3)] for r in Ed["rows"]]))
        s_ = Ed["summary"]
        o.append("\nSpread over all variants at the finest resolution (ny = %d, %d variants with a defined interval): lower edge [%s, %s], upper edge [%s, %s]; offsets from MPB: lower [%s, %s], upper [%s, %s]. Over all resolutions (%d variants): lower offset [%s, %s], upper offset [%s, %s].\n"
                 % (s_["finest_ny"], s_["finest"]["n"], g(s_["finest"]["lo"]["min"], 5), g(s_["finest"]["lo"]["max"], 5), g(s_["finest"]["hi"]["min"], 5), g(s_["finest"]["hi"]["max"], 5),
                    g(s_["finest"]["d_lo"]["min"]), g(s_["finest"]["d_lo"]["max"]), g(s_["finest"]["d_hi"]["min"]), g(s_["finest"]["d_hi"]["max"]),
                    s_["all"]["n"], g(s_["all"]["d_lo"]["min"]), g(s_["all"]["d_lo"]["max"]), g(s_["all"]["d_hi"]["min"]), g(s_["all"]["d_hi"]["max"])))
        o.append("\n**Deep-gap attenuation** (mean ln T over 0.222-0.268 using only frequencies where T > 3e-7 for N = 6, 8, 10)\n\n")
        o.append(table(["ny", "frequencies used", "mean ln T (N=6, 8, 10)", "decay per column"],
                       [[r["ny"], r["n_freq"], ", ".join(g(x, 5) for x in r["mean_lnT"]), g(r["decay_per_column"], 4)] for r in Ed["deep_gap_decay"]]))
        rc = Ed["resolution_change_coarsest_to_finest"]
        o.append("\n**Contributions to the edge uncertainty (thr = 1e-2)**\n\nDiscretization: change of each single-length edge from ny = %d to %d: %s; largest |change| = %s.\n\n"
                 % (rc["ny"][0], rc["ny"][1], "; ".join("N=%d: lo %s, hi %s" % (x["N"], g(x["d_lo"], 3), g(x["d_hi"], 3)) for x in rc["rows"]), g(rc["max_abs"], 3)))
        so = Ed["single_length_offsets_finest"]
        o.append("Finite length / fringes / threshold (single lengths N = 8..16, ny = %d): offset from MPB lower edge in [%s, %s], upper edge in [%s, %s].\n\n"
                 % (rc["ny"][1], g(so["lower"][0], 3), g(so["lower"][1], 3), g(so["upper"][0], 3), g(so["upper"][1], 3)))
        o.append(table(["N", "lo - MPB", "hi - MPB"], [[x["N"], g(x["d_lo"], 3), g(x["d_hi"], 3)] for x in so["rows"]]))
        ob = Ed["odd_band_window"]
        o.append("\n**Symmetry test.** Window [%s, %s] lies above the all-parity gap top and below the even-parity gap top (an odd-Phi band exists in MPB there but cannot couple to a mirror-symmetric plane wave):\n\n"
                 % (g(ob["f_window"][0], 5), g(ob["f_window"][1], 5)))
        o.append(table(["ny", "N", "max T in window", "geometric-mean T"], [[x["ny"], x["N"], g(x["T_max"]), g(x["T_geomean"])] for x in ob["rows"] if x["ny"] in (32, 48)]))
        o.append("\n**Termination dependence** (crystal faces moved outward by face_offset; ny = 32)\n\n")
        o.append(table(["N set", "face offset", "thr", "edge lo", "edge hi", "in-gap peak f", "in-gap peak T", "max T over 0.215-0.283"],
                       [[r["Ns"], g(r["face_offset"]), g(r["thr"]), g(r["edge_lo"], 5), g(r["edge_hi"], 5), g(r["ingap_peak_f"], 4), g(r["ingap_peak_T"]), g(r["ingap_Tmax_overall"])] for r in Ed.get("termination", [])]))
        o.append("\n**Energy conservation R + T - 1 in the long runs** (worst frequency in the band 0.13-0.34 and in the gap 0.22-0.27)\n\n")
        rows_ = []
        for r in At["runs"]:
            f_ = np.array(r["f"]); s_ = np.abs(np.array(r["T"]) + np.array(r["R"]) - 1)
            mb, mg = (f_ >= 0.13) & (f_ <= 0.34), (f_ >= 0.22) & (f_ <= 0.27)
            rows_.append([r["ny"], r["N"], g(r["t_end"], 5), g(float(s_[mb].max())), g(float(f_[mb][np.argmax(s_[mb])]), 4), g(float(np.median(s_[mb]))), g(float(s_[mg].max()))])
        o.append(table(["ny", "N", "t_end", "max band", "at f", "median band", "max gap"], rows_))
        hdr = ["value", "N", "ny", "edge lo -3dB", "edge hi -3dB", "edge lo -20dB", "edge hi -20dB", "T gap geomean", "T gap max", "max |R+T-1| gap", "max |R+T-1| all", "t_end"]
        def rows2(rr):
            return [[g(r["value"]), r["N"], r["ny"], g(r["edge_lo_3dB"], 5), g(r["edge_hi_3dB"], 5), g(r["edge_lo_20dB"], 5), g(r["edge_hi_20dB"], 5),
                     g(r["T_gap_geomean"]), g(r["T_gap_max"]), g(r["max_abs_RplusT_minus_1_gap"]), g(r["max_abs_RplusT_minus_1_all"]), g(r["t_end"], 4)] for r in rr]
        o.append("\n### Single-length (N = 8) numerical-setting sensitivity (edges = first crossing of T = 0.5 / 1e-2 outward from mid-gap; default minimum run time)\n")
        for key, title in (("resolution", "Resolution (ny pixels across sqrt(3))"), ("pml", "PML thickness"), ("decay", "Decay tolerance"),
                           ("source_width", "Source fwidth"), ("placement_monitor", "Monitor distance from faces"), ("placement_source", "Source distance from face")):
            o.append(f"\n**{title}**\n\n")
            o.append(table(hdr, rows2(D[key])))
    o.append("\n## Test run\n\n```\n" + (RES / "pytest_output.txt").read_text()[-1800:] + "\n```\n" if (RES / "pytest_output.txt").exists() else "")
    OUT.write_text("\n".join(o))
    print("wrote", OUT)
