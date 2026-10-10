"""Phase 5D: assemble results/phase5/phase5_handoff.{npz,json} from 5a_array.{npz,json} (selected arrays, definitions, provenance)."""
from _common import RES, np, CASES, THETA0
import json
from pcttd.io import save_json

z = np.load(RES / "5a_array.npz"); j = json.load(open(RES / "5a_array.json"))
keep = dict(f=z["f"], mask_common=np.ones(len(z["f"]), bool), theta_grid=z["theta_grid"], d=z["d"], theta0=z["theta0"], f0=z["f0"], dtau_target=z["dtau_target"], window=z["window"],
            pattern_freqs=z["pat_freqs"], patterns_res30=z["patterns_res30"], pattern_cases=z["pattern_cases"])
for res in (30, 40):
    for k in CASES + ["raw_unit", "trim_unit"]:
        for q in ("theta_peak", "gain", "hpbw", "strongest_db", "strongest_theta", "n_3db"):
            keep[f"res{res}_{k}_{q}"] = z[f"res{res}_{k}_{q}"]
        keep[f"res{res}_{k}_pointing_error"] = z[f"res{res}_{k}_theta_peak"] - THETA0
    for k in ("raw", "trim"):
        keep[f"res{res}_{k}_phase_error_vs_ideal"] = z[f"res{res}_{k}_eps"]
np.savez(RES / "phase5_handoff.npz", **keep)
meta = dict(source="results/phase5/5a_array.{npz,json}; Phase 4 inputs results/phase4/phase4_handoff.npz (unchanged)",
            frequency="mask_common 0.231-0.268 (domain-limited); array arrays are over this band only", cases=CASES, labels={"trim": "phase-intercept-trimmed PhC (never ideal TTD)", "raw": "raw PhC"},
            **{k: j[k] for k in ("setup", "definitions", "restrictions", "30", "40", "res_check", "amp_diag", "fixed_analytic_check_max_dev_deg")})
save_json(RES / "phase5_handoff.json", meta)
print(len(keep), "arrays")
