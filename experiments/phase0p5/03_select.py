"""Apply the pre-declared gates, criteria and decision rule (docs/phase0p5_criteria.md, Revisions 1-2)."""
import numpy as np
from _common import *  # noqa
from pcttd import io as pio

MARGIN = 0.10


def metrics(cfg, bulk):
    ds = [d for d in cfg["designs"] if "skipped" not in d]
    win = [d for d in ds if d.get("window")]
    lat = cfg["lattice"]
    rel_gap = next(g["rel"] for rw in bulk[lat]["rows"] if abs(rw["r"] - cfg["r_base"]) < 1e-9
                   for g in rw["gaps"] if g["n"] == cfg["n_gap"]) if any(
        abs(rw["r"] - cfg["r_base"]) < 1e-9 for rw in bulk[lat]["rows"]) else None
    if rel_gap is None:   # supplementary radius not in the bulk scan grid: use the sweep-resolution gap
        g = cfg["gap_at_sweep_res"]; rel_gap = (g["w_hi"] - g["w_lo"]) / (0.5 * (g["w_hi"] + g["w_lo"]))
    ng_max = [d["ng_max_strict_even"] for d in win]
    m = dict(
        n_designs=len(cfg["designs"]), n_valid=len(ds), n_with_window=len(win),
        rel_gap=rel_gap,
        c2_max_sm_frac=max([d["sm_strict_even_frac"] for d in ds], default=0.0),
        c3_median_strip=float(np.median([d["window"]["strip_at_center"] for d in win])) if win else None,
        c4_max_gbp10=max([d["gbp10"] for d in ds], default=0.0),
        max_gbp5=max([(d["flat_window"]["5"] or {}).get("gbp", 0.0) for d in ds], default=0.0),
        max_gbp20=max([(d["flat_window"]["20"] or {}).get("gbp", 0.0) for d in ds], default=0.0),
        ng_max_overall=max(ng_max, default=0.0),
        ng_max_min_over_windowed=min(ng_max, default=0.0),
        c5_pixels_per_period=(min(d["r_min"] for d in ds) and
                              ds[0]["H"] * (8.0 / min(d["r_min"] for d in ds)) ** 2) if ds else None,
        c5_res_req=8.0 / min(d["r_min"] for d in ds) if ds else None,
        c5_median_wall_s=float(np.median([d["wall_s"] for d in ds])) if ds else None,
        c6_fraction=float(np.mean([d["c6_pass"] for d in ds])) if ds else 0.0,
    )
    g2 = [d for d in ds if d["sm_strict_even_frac"] >= 0.01]
    g3 = [d for d in g2 if d["window"]["vg_sign_constant"] and d["window"]["edge_max"] < 1e-3]
    wi = [d["ng_max_strict_even"] for d in win if d["ng_max_strict_even"] > 0]
    gates = dict(
        G1=bool(m["rel_gap"] >= 0.05), G2=bool(g2), G3=bool(g3),
        G4=bool(wi and max(wi) >= 10 and max(wi) / min(wi) >= 2),
    )
    gates["detail"] = dict(G2_designs=[(d["f1"], d["s"]) for d in g2], G3_designs=[(d["f1"], d["s"]) for d in g3])
    return m, gates


def compare(a, b, key, higher_better=True):
    x, y = a[key], b[key]
    if x is None or y is None:
        return "tie"
    if higher_better:
        if x > y * (1 + MARGIN): return "A"
        if y > x * (1 + MARGIN): return "B"
    else:
        if x < y * (1 - MARGIN): return "A"
        if y < x * (1 - MARGIN): return "B"
    return "tie"


CRIT = [("C1", "rel_gap", True), ("C2", "c2_max_sm_frac", True), ("C3", "c3_median_strip", True),
        ("C4", "c4_max_gbp10", True), ("C5", "c5_pixels_per_period", False), ("C6", "c6_fraction", True)]


def decide(name_a, ma, ga, name_b, mb, gb):
    elig_a, elig_b = all(ga[k] for k in ("G1", "G2", "G3", "G4")), all(gb[k] for k in ("G1", "G2", "G3", "G4"))
    out = dict(eligible={name_a: elig_a, name_b: elig_b}, per_criterion={}, wins={name_a: 0, name_b: 0})
    for c, key, hb in CRIT:
        r = compare(ma, mb, key, hb)
        out["per_criterion"][c] = dict(metric=key, A=ma[key], B=mb[key], winner={"A": name_a, "B": name_b, "tie": "tie"}[r])
        if r == "A": out["wins"][name_a] += 1
        if r == "B": out["wins"][name_b] += 1
    if elig_a and not elig_b: out["selected"], out["rule"] = name_a, "only eligible candidate"
    elif elig_b and not elig_a: out["selected"], out["rule"] = name_b, "only eligible candidate"
    elif not elig_a and not elig_b: out["selected"], out["rule"] = None, "neither eligible: no selection, refer to review"
    else:
        wa, wb = out["wins"][name_a], out["wins"][name_b]
        if wa != wb:
            out["selected"], out["rule"] = (name_a if wa > wb else name_b), "most criteria won"
        else:
            for c in ("C4", "C2"):
                w = out["per_criterion"][c]["winner"]
                if w != "tie":
                    out["selected"], out["rule"] = w, f"tie on count; tie-break {c}"
                    break
            else:
                out["selected"], out["rule"] = None, "complete tie: refer to review"
    return out


if __name__ == "__main__":
    bulk = pio.load_json(RES / "bulk_gap_scan.json")
    sweep = pio.load_json(RES / "defect_sweep.json")
    M_, G_ = {}, {}
    for n, cfg in ((k, v) for k, v in sweep.items() if k != "provenance"):
        M_[n], G_[n] = metrics(cfg, bulk)
    dec = decide("H", M_["H"], G_["H"], "R", M_["R"], G_["R"])
    sens = decide("H_r030", M_["H_r030"], G_["H_r030"], "R", M_["R"], G_["R"])
    out = dict(metrics=M_, gates=G_, decision=dec, sensitivity_H_r030_vs_R=sens,
               decision_robust_to_H_radius=(dec["selected"] == ("H" if sens["selected"] == "H_r030" else sens["selected"])))
    pio.save_json(RES / "selection.json", out)
    import json; print(json.dumps(out, indent=1, default=float))
