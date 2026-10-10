# Phase 7 audit — repository, scientific chain, claim wording

Starting checkpoint: HEAD 3c10780, clean tree, `git diff --check` clean, 195 tests passed (0 failed, 0 skipped, 145 s), Phase 4/5/6 handoffs present. No discrepancy with the expected starting state. No new simulation was run in Phase 7; all new numbers are read from existing result files by `experiments/phase7/7a_final_results.py` (`results/final/final_results.json`, 149 values, each with source file and key). The headline values were also recomputed independently from the Phase 4 handoff (band-fit increments 3.9315/3.9308/3.9310 a/c, trimmed residual 0.02597 rad, amplitude imbalance 0.0227 dB) and agree with the Phase 4 report.

## 1. Scientific chain

| phase | question | method | result | restriction carried | feeds |
|---|---|---|---|---|---|
| 0.5 | which 2-D family | MPB screening with a pre-declared rule | triangular air holes, TE, W1; moderate-margin decision (R wins single-mode bandwidth; H wins on slow-light potential and FDTD cost); r/a = 0.30 is a post-selection engineering choice | flat n_g >= 10 windows are about 1 % wide: early warning that broadband slow light is hard | geometry for all later phases |
| 1 | is the FDTD extraction trustworthy | vacuum, dielectric, Bragg vs TMM, 2-D crystal vs MPB | second-order convergence; edges within 3.5e-3 of MPB | edge tolerance dominated by finite length/termination; no W1 yet | pipeline and estimator rules |
| 2 | does finite-W1 FDTD reproduce MPB | modal projection, 3-length fit | median n_g bias +0.3 %, max about 1 % in 0.231-0.268 | upper edge 0.268 chosen post hoc; tol 1e-7 required | validated domain |
| 3 | does first-row engineering pay | D0/D1/D2 | D1/D2: +7 %/+14 % delay, -43 %/-62 % bandwidth; D0 selected | three-point illustration; sign bug found, runs discarded, rerun | D0 for the family |
| 4 | four-channel family | N = 8-11 direct FDTD | increments within 0.04 % (band fit); imbalance 0.023 dB; raw intercept error up to 2.7 rad; trimmed residual 0.026 rad | first attempt discarded; run-time and estimator changes post hoc | actual H_m |
| 5 | element-level array | H_m as weights | trimmed PhC within 0.05 deg of ideal TTD; but pitch 1.96 lambda0 gives 0 dB grating lobes | coarse quantum (0.981 cycle) | motivates Phase 6 |
| 6 | mitigation | subarrays, Monte Carlo | Msub = 4 matches the quantum to 1.9 %; 0.189 deg; lobes mitigated | one angle; trim; replica lobes | final claim |

Overclaim check by phase: Phase 0.5 flags its own margin; Phase 1 and 2 are labeled verification and claim nothing about TTD; Phase 3 labels D1/D2 a three-point illustration; Phase 4 states that the favorable numbers reflect one flat-n_g mode at four lengths and says nothing about beams; Phase 5 declines to claim grating-lobe suppression; Phase 6 labels budgets as system-level. No phase overclaims beyond its numerical evidence. Two wording issues remain outside the phase reports (Section 4: README title and the term "dispersion-engineered").

## 2. Final claim audit

Candidate claim (from the task): "A dispersion-engineered W1 photonic-crystal waveguide can provide broadband true-time-delay behavior suitable for subarray-level phased-array beamforming. ... naturally matches ... avoiding the severe element-level grating-lobe problem ..."

| phrase | verdict | evidence | replacement |
|---|---|---|---|
| "dispersion-engineered" | **not justified, dropped** | the selected device D0 is the plain W1. The engineered first-row variants D1/D2 were evaluated and rejected (narrower band, 4-8x GDD). The only dispersion engineering in the project is a negative result | "W1 photonic-crystal delay line", "dispersion-characterized plain W1". Say explicitly that first-row engineering did not help |
| "broadband" | qualified | 14.8 % fractional (0.231-0.268), domain-limited lower bound, upper edge a post hoc validation limit; pointwise delay increments vary by several percent inside it | "14.8 % validated fractional bandwidth" every time |
| "true-time-delay" | correct only for the slope | delay increments equal the target to 0.04 % (band fit) but raw channels carry a constant phase intercept error up to 2.7 rad; trimmed channels have residual 0.026 rad that is mostly n_g curvature | "delay lines with TTD slope plus one constant phase trim per channel"; never call the trimmed family "ideal TTD". Keep "TTD" only with that qualifier |
| "suitable" | too strong | 2-D, lossless, no disorder, no couplers, no experiment | "supports the feasibility of" / "promising in simulation" |
| "naturally matches" | softened | the match is delay per period x f0 = n_g f0 = 0.981 vs required 1.000 at theta0 = 30 deg, Msub = 4. It is a numerical coincidence of one geometry, one carrier and one angle, not a design principle. A conventional line can be cut to any length | "happens to match to 1.9 %" |
| "avoiding the grating-lobe problem" | softened | equal-strength 0 dB lobes become ordinary sidelobes (-12.4 dB) plus replica lobes <= -16.6 dB; replicas remain | "mitigates"; "replaces ... by" |
| "beam-squint reduction" | supported | max pointing error 0.189 deg vs 2.69 deg fixed phase; R_squint 0.070; res 30/40 agree to 0.007 deg | keep, with the "even ideal subarray TTD has 0.156 deg" caveat |
| "requires constant phase trim" | supported | raw offset about -6.7 deg; Phase 4 trim used unchanged gives 0.73 deg | keep |

**Settled central claim.** In a two-dimensional, lossless FDTD/MPB model, finite sections of a plain W1 photonic-crystal waveguide provide delay increments of 3.93 a/c per lattice period (0.981 carrier cycle at band center, accurate to 0.04 % by band fit) over a validated 14.8 % fractional bandwidth. Using each of four channels (N = 8-11) for a four-element lambda0/2 subarray, which requires 1.000 cycle at 30 deg (a 1.9 % mismatch), and one constant phase trim per channel, the maximum pointing error over the band is 0.19 deg against 2.69 deg for fixed phase, and the Phase 5 equal-strength grating lobes become ordinary sidelobes (-12.4 dB) with residual subarray replica lobes (<= -16.6 dB). This supports the feasibility of subarray-level photonic-crystal delay lines for beam-squint reduction; it does not demonstrate a device. It is limited by the 2-D lossless model, the domain-limited band, the required trim, residual replica lobes, one steering angle, and the absence of loss, disorder, coupling and experiment.

## 3. Phase 5 interpretation and the strongest result

The Phase 5 outcome is interpreted as a **limitation imposed by the delay quantum**, not as a failure of the delay lines: element-level trimmed PhC tracked ideal TTD to 0.05 deg. The array failed to be a practical aperture because the quantum fixes d = 1.96 lambda0. This is fair to both phases and is how the manuscript presents it.

Which result is strongest? Quantitatively cleanest: the Phase 4/5 agreement between simulated channels and ideal TTD slope (0.04 %, 0.05 deg). Architecturally most useful: Phase 6 Msub = 4. The central claim rests on Phase 6, with Phase 4 as its numerical foundation. Phase 6 should not be oversold: the pointing improvement over fixed phase is the generic benefit of any TTD, and the PhC contributes only 0.055 deg beyond ideal subarray TTD (0.156 deg). What the study shows is that real simulated PhC channels reach ideal-TTD performance at this subarray size, not that PhC lines beat conventional ones.

## 4. Terminology and wording audit (findings and fixes)

* GRATING LOBE: Phase 5 element-level spatial-alias lobes (0 dB), Phase 6 subarray replica lobes (<= -16.6 dB), ordinary sidelobes (-12.4 dB) are kept distinct in all Phase 7 documents.
* "Robustness"/"tolerance": Phase 6 numbers are called system-level tolerance estimates; no fabrication claim.
* README title said "dispersion-engineered ... as true-time-delay elements" and the status table was stale (Phase 5 row "Phase 6 not started", "5-7 not started"). Fixed in this phase.
* `docs/research_plan.md` uses "dispersion-engineered" in the founding question. It is the approved Phase 0 text and is kept, with a note at its head.
* Phase 4 and Phase 5 report status lines ("Phase 5 has not been started"; Phase 5 grating lobes unresolved) are historical; superseding notes were added at their heads. Their bodies are unchanged.

## 5. Contradictions and inconsistencies found

1. README status stale and "dispersion-engineered" title: **fixed**.
2. README had no Phase 2 reproduction section; the Phase 2 command lines are inferred from script usage strings and the stored per-run configuration: **documented in `docs/reproducibility.md` (marked inferred)**.
3. Phase 6 report says "4 visible at f0" replica lobes; the analytic list in `6a_subarrays.json` has 3 (-30, 0, 90 deg). The metric counts all q != 0 with |sin| <= 1, which includes two grazing lobes at exactly +-90 deg (q = -3 and q = +1, sin = -1 and +1); the analytic list drops -90 deg through a floating-point boundary test. Phase 5 (decision D13) does not count +-90 deg endpoint skirts. **No numbers change**; the -16.6 dB headline is the maximum over the larger set, so excluding grazing lobes could only lower it. Noted in the manuscript limitations.
4. Phase 3 D0 local GDD at band center is +0.77 (local quadratic) vs -0.91 (band-average fit) and -0.79 (MPB): sign disagreement is documented in Phase 3 as FDTD local-GDD noise; the final documents use only the band-average and MPB values.
5. Trimmed pointing error looks larger in Phase 6 (0.189 deg) than Phase 5 (0.050 deg) although the same channels are used: different arrays (the Phase 6 fixed intra-subarray phases squint, 0.156 deg even for ideal subarray TTD). Explained in the final tables; not a contradiction.
6. Pointwise delay increments vary by up to about 6 % inside the band (cubic-fit curve, Figure 4b; up to 9.7 % pointwise in the Phase 4 report) while the band-fit increments are within 0.04 %. Both are kept and shown; the beam metrics are not sensitive because pointing depends on the phase slope.
7. Physical mapping tables round lambda0 differently (9.99 vs 10.0 mm at 30 GHz); cosmetic.
8. Unused code: `src/pcttd/edges.py` is exercised only by Phase 1 (the slope-of-ln-T estimator, validated on TMM and found not applicable to the 2-D gap). It is kept and labeled. No other dead modules were found (every other module is imported by at least two files).
9. No doc references a missing result file or script (checked by a reference scanner over all `docs/*.md` and the README, repeated as a test). Two kinds of file are intentionally untracked: raw plane-field dumps `results/phase{2,3,4}/raw*` (git-ignored, 134 MB; processed outputs are tracked) and `__pycache__`.
10. Tracked results from discarded methodology: `results/phase4/attempt1_*` (under-converged first attempt) are kept deliberately and are labeled in the Phase 4 report; the Phase 3 sign-bug runs were discarded and rerun (Phase 3 report Section 13: not used anywhere); the Phase 3 result files describe the reruns.

No unresolved contradiction affects a headline number.

## 6. Answers to the final audit questions

1. **Chain internally consistent?** Yes (Section 1); items in Section 5 are presentation-level.
2. **Final claim supported?** The settled claim in Section 2 is; the original wording was not (see table).
3. **"True-time delay" correct?** For the delay slope, yes; for the trimmed family only as "TTD slope plus constant trim".
4. **"Dispersion-engineered" justified?** No. Dropped.
5. **Validated bandwidth described accurately?** Yes, as 14.8 %, domain-limited, post hoc upper edge.
6. **Phase 4-6 restrictions preserved?** Yes; they are Table 7 of `final_tables.md` and the manuscript limitations.
7. **Numerical corrections disclosed?** Yes: Phase 2 post hoc boundary, Phase 3 sign bug, Phase 4 discarded attempt, run-time change, estimator change and cascade step-count erratum, Phase 6 first Monte Carlo sign issue and post hoc diagnostics.
8. **Phase 5 interpreted fairly?** Yes, as a quantum-imposed limitation.
9. **Phase 6 the strongest central result?** Strongest architectural result, yes; Phase 4 is the strongest quantitative result (Section 3).
10. **Reproducible?** Figures, tables, handoff and the Phase 5/6 analyses: yes in minutes from tracked results. Phases 2-4 FDTD: reproducible in principle (hours; commands documented, Phase 2 inferred) but raw dumps are not tracked.
11. **Figures publication-quality?** Nine multi-panel figures, colorblind-safe palette, PNG (300 dpi) and PDF; reviewed visually. Figure 1 is schematic.
12. **Reviewable as a draft paper?** Yes, as a computational feasibility study (`docs/manuscript_draft.md`), with literature placeholders.
13. **Executive summary understandable by a hiring manager?** Written for that audience (`docs/executive_summary.md`).
14. **Unresolved contradictions?** None affecting results.
15. **More simulation needed before arXiv?** No. See Section 7.

## 7. arXiv readiness: B — NEARLY READY, small gap remains

* **Manuscript-writing gap:** all `[REF]` placeholders; literature positioning (no novelty claim has been checked against prior art); phrasing the physical mapping; a copy edit. None needs simulation.
* **Scientific gap:** the project's own research plan made the conventional-waveguide baseline "core", and the only baseline is MPB-only (Phase 3: the H-family W1 needs about 2.1x fewer wavelengths of length than the R reference for the same delay at comparable bandwidth). No loss, coupling or subarray-level conventional comparison exists, so the manuscript cannot claim an advantage over conventional delay lines beyond group-index footprint; it claims feasibility only. This restricts the claim; it does not block a feasibility draft.
* **Numerical validation gap:** none blocking. Known and disclosed limits (domain-limited band, noise-limited pointwise derivatives, two-resolution convergence) are stated.

Category A would require the literature check; category C does not apply because no result is unvalidated.
