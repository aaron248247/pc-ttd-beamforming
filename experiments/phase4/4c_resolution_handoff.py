"""Phase 4C: res 30 vs res 40 comparison (criteria Revision 4 item 8) and the single machine-readable Phase 5 handoff file.
usage: python 4c_resolution_handoff.py [--handoff]    reads results/phase4/4b_res30.{npz,json}, 4b_res40.{npz,json}
-> results/phase4/4c_resolution.json ; with --handoff also results/phase4/phase4_handoff.npz + phase4_handoff.json
"""
from _common import RES, ROOT, np, F0, DTAU_TARGET, VALID, LENGTHS
import json, sys
from pcttd import channels as C
from pcttd.io import save_json, provenance

LABEL = "primary_poly_tol10"


def load(res):
    return dict(np.load(RES / f"4b_res{res}.npz")), json.load(open(RES / f"4b_res{res}.json"))


if __name__ == "__main__":
    a, ja = load(30); b, jb = load(40)
    f = a["f"]
    out = dict(rule="criteria Revision 4 item 8", checks={})
    ba, bb = ja["bands"][LABEL], jb["bands"][LABEL]
    sel = (f >= VALID[0] - 1e-9) & (f <= VALID[1] + 1e-9)
    # evaluate both resolutions on the SAME frequency window (the res-30 common band) so differences are not a band-edge artifact
    sl = slice(ba["i_low"], ba["i_high"] + 1) if ba else slice(None)
    d_inc = [float(np.mean(a["inc_poly"][m][sel]) - np.mean(b["inc_poly"][m][sel])) / DTAU_TARGET for m in range(3)]
    d_inc_c = [float(np.mean(a["inc"][m][sel]) - np.mean(b["inc"][m][sel])) / DTAU_TARGET for m in range(3)]
    out["band_res30"] = ba and [ba["f_low"], ba["f_high"]]; out["band_res40"] = bb and [bb["f_low"], bb["f_high"]]
    out["mean_increment_change_over_target_valid_domain_poly"] = d_inc
    out["mean_increment_change_over_target_valid_domain_central"] = d_inc_c
    out["amp_imbalance_db_max_valid"] = [float(a["imbalance_db"][sel].max()), float(b["imbalance_db"][sel].max())]
    out["trim_delta_res30"] = [float(x) for x in a["delta"]]; out["trim_delta_res40"] = [float(x) for x in b["delta"]]
    out["eps_trim_max_valid"] = [float(np.abs(a["eps_phi_trim"][:, sel]).max()), float(np.abs(b["eps_phi_trim"][:, sel]).max())]
    out["eps_raw_diff_max_valid"] = float(np.abs(a["eps_phi_raw"][:, sel] - b["eps_phi_raw"][:, sel]).max())
    # direct comparison of the complex transfer functions (relative to channel 0, i.e. the quantity the array uses)
    rel_a = a["H"] / a["H"][0]; rel_b = b["H"] / b["H"][0]
    out["max_abs_phase_diff_rel_valid"] = float(np.abs(np.angle(rel_a[:, sel] / rel_b[:, sel])).max())
    out["max_amp_rel_diff_valid"] = float(np.abs(np.abs(a["H"][:, sel]) / np.abs(b["H"][:, sel]) - 1).max())
    out["t_end_res30"] = ja["t_end"]; out["t_end_res40"] = jb["t_end"]
    out["tests"] = dict(
        mean_inc_change_le_2pct=bool(max(abs(x) for x in d_inc) <= 0.02),
        band_edges_move_le_0p002=bool(ba and bb and abs(ba["f_low"] - bb["f_low"]) <= 0.002 and abs(ba["f_high"] - bb["f_high"]) <= 0.002),
        imbalance_change_le_0p1dB=bool(abs(out["amp_imbalance_db_max_valid"][0] - out["amp_imbalance_db_max_valid"][1]) <= 0.1),
        eps_trim_change_le_0p1rad=bool(abs(out["eps_trim_max_valid"][0] - out["eps_trim_max_valid"][1]) <= 0.1))
    save_json(RES / "4c_resolution.json", out, with_provenance=True)
    print(json.dumps(out, indent=1)[:3000])
    if "--handoff" in sys.argv:
        prim = ja["bands"][LABEL]
        mask_common = a["mask_common"].astype(bool)
        arrs = dict(
            f=f, H_raw=a["H"], amp_raw=np.abs(a["H"]), phase_raw_wrapped=np.angle(a["H"]), phase_abs_unwrapped=a["phi_abs"], phase_rel_unwrapped=a["phi_rel"],
            H_trim=a["H_trim"], trim_delta=a["delta"], S11_raw=a["S11"], T=a["T"], R=a["R"], proj_res_in=a["res_in"], proj_res_out=a["res_out"], unitarity_defect=a["defect"],
            mask_common=mask_common, mask_validated_domain=sel, f0=F0, dtau_target=DTAU_TARGET, N=np.array(LENGTHS), amp_imbalance_db=a["imbalance_db"],
            H_raw_res40=b["H"], H_trim_res40=b["H_trim"])
        prov = dict(source="results/phase4/4b_res30.npz (production, res 30; res 40 confirmation arrays carry the suffix _res40)",
                    sign_convention=("MEEP exp(-i omega t), MEEP_SIGN = +1: a delay tau multiplies a monochromatic signal by exp(+i 2 pi f tau); tau_g = +d arg H / d omega; "
                                     "ideal TTD channel relative to channel 0: H_m/H_0 = exp(+i 2 pi f m dtau_target); array factor AF = sum_m H_m exp(-i k0 x_m sin(theta)) (pcttd.array)"),
                    H_definition=("H_m(f) = S21_m = a_f(out)/a_f(in), guided-mode amplitudes at x = -3 and x = N_m + 3 of the directly simulated D0 device with N_m periods; "
                                  "includes the full propagation phase over N_m + 6 periods; nothing de-embedded; channel 0 is the reference for relative phases"),
                    H_trim_definition="H_trim_m = H_raw_m exp(-i trim_delta_m); trim_delta_m = raw phase error at f0 (one frequency-independent constant per channel; amplitude and group delay untouched)",
                    labels="H_raw = raw PhC; H_trim = phase-intercept-trimmed PhC (never call either ideal TTD)",
                    N=LENGTHS, N0=LENGTHS[0], dN=LENGTHS[1] - LENGTHS[0], f0=F0, dtau_target=DTAU_TARGET, validated_domain=list(VALID),
                    common_band=[prim["f_low"], prim["f_high"]] if prim else None, common_band_definition=f"4b_res30.json bands.{LABEL} (see report for the final gate band)",
                    frequency_grid="0.220-0.275 step 5e-4 (111 points); use mask_common / mask_validated_domain", run_config=ja["cfg"], run_provenance=ja["run_provenance"],
                    analysis_provenance=provenance(), restrictions=["quantitative claims only inside mask_validated_domain (0.231-0.268, post hoc Phase 2 validation limit)",
                                                                    "delay granularity: one period ~ one carrier cycle", "pointwise delay/GDD from H is noise-limited"])
        np.savez_compressed(RES / "phase4_handoff.npz", **arrs, provenance_json=json.dumps(prov))
        save_json(RES / "phase4_handoff.json", prov, with_provenance=False)
        print("handoff written", list(arrs))
