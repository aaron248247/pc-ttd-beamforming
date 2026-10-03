"""Generate docs/phase0p5_report.md. Every number is read from results/phase0p5/*.json; none is typed by hand."""
import json
from _common import *  # noqa
from pcttd import io as pio

OUT = ROOT / "docs" / "phase0p5_report_tables.md"


def f(x, n=3):
    return "—" if x is None else (f"{x:.{n}f}" if isinstance(x, float) else str(x))


def table(head, rows):
    s = "| " + " | ".join(head) + " |\n|" + "|".join("---" for _ in head) + "|\n"
    return s + "\n".join("| " + " | ".join(map(str, r)) + " |" for r in rows) + "\n"


if __name__ == "__main__":
    bulk = pio.load_json(RES / "bulk_gap_scan.json")
    sel = pio.load_json(RES / "selection.json")
    sw = pio.load_json(RES / "defect_sweep.json")
    sw.update({k: v for k, v in pio.load_json(RES / "defect_sweep_supp.json").items() if k != "provenance"})
    env = pio.load_json(ROOT / "environment" / "env_versions.json")
    o = []
    o.append("## T1. Bulk gap (declared rule, res %d)\n" % BULK_RES)
    rows = []
    for lat, pol in (("triangular", "TE (Hz)"), ("square", "TM (Ez)")):
        b = bulk[lat]
        rw = next(r for r in b["rows"] if abs(r["r"] - b["r_base"]) < 1e-9)
        g = next(g for g in rw["gaps"] if g["n"] == b["n_gap"])
        rows.append([lat, pol, b["n_gap"], f(b["r_base"], 3), f(g["w_lo"]), f(g["w_hi"]), f(g["rel"])])
    o.append(table(["lattice", "pol.", "gap n|n+1", "r_base/a", "w_lo", "w_hi", "rel. width"], rows))
    o.append("\nRelative gap vs radius, triangular (n=1): " + ", ".join(
        f"r={r['r']:.3f}:{next((g['rel'] for g in r['gaps'] if g['n']==1), 0):.3f}" for r in bulk['triangular']['rows'][::2]) + "\n")

    o.append("\n## T2. Gates and sweep coverage\n")
    rows = []
    for n, g in sel["gates"].items():
        cfg = sw[n]
        nun = sum(1 for d in cfg["designs"] if "skipped" in d)
        rows.append([n, cfg["lattice"], f(cfg["r_base"], 3), cfg["M"], cfg.get("extra_bands", 8), len(cfg["designs"]) - nun, nun,
                     *["PASS" if g[k] else "FAIL" for k in ("G1", "G2", "G3", "G4")]])
    o.append(table(["config", "lattice", "r/a", "M", "extra bands", "solved", "unresolved (excluded)", "G1", "G2", "G3", "G4"], rows))

    o.append("\n## T3. Criterion metrics\n")
    rows = []
    for n, m in sel["metrics"].items():
        rows.append([n, f(m["rel_gap"]), f(m["c2_max_sm_frac"]), f(m["c3_median_strip"], 4), f(m["c4_max_gbp10"], 4),
                     f(m["max_gbp5"], 4), f(m["max_gbp20"], 4), f"{m['c5_pixels_per_period']:.0f}", f(m["c6_fraction"], 2),
                     f"{m['ng_max_overall']:.0f}"])
    o.append(table(["config", "C1 gap", "C2 SM frac", "C3 strip", "C4 GBP(n≥10)", "GBP(n≥5)", "GBP(n≥20)", "C5 px/period",
                    "C6 frac", "max n_g (singular at band extrema)"], rows))

    o.append("\n## T4. Decision per pairing (pre-declared rule)\n")
    rows = []
    for k, d in sel["decisions"].items():
        w = ", ".join(f"{c}:{v['winner']}" for c, v in d["per_criterion"].items())
        rows.append([k, str(d["eligible"]), w, d["selected"], d["rule"]])
    o.append(table(["pairing", "eligible", "per-criterion winners", "selected", "rule"], rows))
    o.append(f"\nAll pairings select the same family: **{sel['unanimous']}** ({', '.join(sel['families_selected'])}).\n")

    o.append("\n## T5. Plain W1 and best-GBP design per configuration\n")
    rows = []
    for n, cfg in sw.items():
        if n == "provenance": continue
        ds = [d for d in cfg["designs"] if "skipped" not in d]
        plain = next((d for d in ds if d["f1"] == 1.0 and d["s"] == 0.0), None)
        best = max(ds, key=lambda d: d["gbp10"])
        for lab, d in (("plain W1", plain), ("best GBP(≥10)", best)):
            if d is None: continue
            fw = d["flat_window"]["10"]
            rows.append([n, lab, f"{d['f1']:g}/{d['s']:+.2f}", d["max_modes_in_gap"], f(d["sm_strict_even_frac"]),
                         f(d["sm_even_only_frac"]), (f"{fw['n_mean']:.1f}" if fw else "—"),
                         (f(fw['frac_bw']) if fw else "—"), (f(fw['gbp']) if fw else "—"),
                         f(d["window"]["strip_at_center"], 4) if d["window"] else "—",
                         f"{d['window']['edge_max']:.1e}" if d["window"] else "—",
                         (f"{d['window']['ng_at_center']:.2f} [{d['window']['ng_min']:.2f}–{d['window']['ng_max']:.0f}]"
                          if d["window"] else "—")])
    o.append(table(["config", "design", "f1/s", "max modes in gap", "SM-strict(even) frac", "SM even-only frac",
                    "flat-window n_g", "its frac. BW", "its GBP", "strip frac", "edge frac (window max)",
                    "n_g at window center [min–max]"], rows))

    import glob
    cfiles = sorted(glob.glob(str(RES / "convergence_*.json")))
    if cfiles:
        o.append("\n## T6. Numerical convergence (same analysis code)\n")
        rows = []
        for cf in cfiles:
            for r in pio.load_json(cf)["rows"]:
                def rel(a, b): return "—" if (a in (None, 0) or b is None) else f"{(b - a) / a * 100:+.2f}%"
                rows.append([r["config"], f"{r['f1']:g}/{r['s']:+.2f}", r["case"], r["bands_sufficient"],
                             rel(r["gap_w_lo_base"], r["gap_w_lo_case"]), rel(r["gap_w_hi_base"], r["gap_w_hi_case"]),
                             f(r["sm_strict_even_frac_base"]) + " → " + f(r["sm_strict_even_frac_case"]),
                             f(r["gbp10_base"], 4) + " → " + f(r["gbp10_case"], 4),
                             (f"{r['ng_at_center_base']:.2f} → {r['ng_at_center_case']:.2f}" if "ng_at_center_base" in r else "—"),
                             f"{r['vg_fd_rel_err_case']['median']:.3f}" if r["vg_fd_rel_err_case"]["median"] is not None else "—"])
        o.append(table(["config", "f1/s", "case vs base", "bands ok", "Δ gap_lo", "Δ gap_hi", "SM-even frac base → case",
                        "GBP(≥10) base → case", "n_g at window center", "median |v_g,MPB−FD|/v_g"], rows))
    o.append("\n## Environment\n\n" + "```\n" + json.dumps(env, indent=1) + "\n```\n")
    o.append("\nSweep provenance: `" + json.dumps(sel.get("provenance", {})) + "`\n")
    o.append("\n## T7. Supercell edge-energy fraction over all solved designs (window max)\n")
    rows = []
    for n in ("R", "R_M6", "H", "H_r030"):
        e = [d["window"]["edge_max"] for d in sw[n]["designs"] if "skipped" not in d and d.get("window")]
        rows.append([n, sw[n]["M"], f"{min(e):.1e}", f"{max(e):.1e}", sum(1 for x in e if x < 1e-3), len(e)])
    o.append(table(["config", "M", "min", "max", "designs with edge < 1e-3", "designs with window"], rows))
    OUT.write_text("\n".join(o))
    narr = ROOT / "docs" / "phase0p5_report_narrative.md"
    if narr.exists():
        (ROOT / "docs" / "phase0p5_report.md").write_text(narr.read_text() + "\n\n---\n# Appendix: generated tables\n\nAll values below are read from `results/phase0p5/*.json` by `experiments/phase0p5/06_report.py`.\n\n" + "\n".join(o))
    print("wrote", OUT)
