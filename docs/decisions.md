# Project decision log

Chronological, append-only. Phase 0.5 documents (`phase0p5_criteria.md` with Revisions 1-3, `phase0p5_report*.md`) are
the historical record of the geometry selection and are not rewritten. (A spelling-only normalization to American English was applied on 2026-10-03.)

## 2026-10-03 — Review of Phase 0.5

**D1. Geometry family (confirmed).** Triangular lattice, air holes in epsilon = 12, TE (Hz out of plane), W1 line defect.
This is the photonic-crystal architecture for the next phases.

**D2. Baseline radius: r/a = 0.30 — a post-selection engineering baseline choice within the selected H family.**
The pre-declared Phase 0.5 rule selected the H family and placed its gap-maximizing baseline at r/a = 0.40 (relative bulk gap
0.490, on the upper limit of the declared scan range). **That historical result is unchanged.** r/a = 0.30 is *not* the output
of the selection rule. It is an engineering choice made afterwards, for these reasons:
* it retains a substantial bulk band gap (relative width 0.280, Phase 0.5 table T3, config `H_r030`);
* all 25 screened defect geometries were numerically resolved at the screening resolution;
* it provides richer slow-light tunability in the Phase 0.5 screening sweep;
* it avoids relying on the extreme upper boundary of the radius scan;
* the purpose of the project is controllable true-time-delay behavior, not maximizing bulk gap width.

The r/a = 0.40 results remain part of the Phase 0.5 record and its sensitivity analysis. Config name `H_r030` in the Phase 0.5
result files is a compatibility-sensitive identifier and is kept.

**D3. Conventional-waveguide reference (carried forward).** The R air-channel waveguide (square-lattice rods, TM, W1)
is retained as a *reference* waveguide, not as a competitor for selection as the primary PhC architecture. It will anchor
later comparisons of delay per physical length, usable fractional bandwidth, transmission, reflection, delay ripple, GDD,
sensitivity, and eventually beam-squint performance. The project will quantify what the slow-light H structure gains and
sacrifices relative to this lower-dispersion reference; no system-level advantage of H is assumed.

**D4. Language.** American English spelling in all prose, comments, docstrings, documentation, figure labels, captions,
reports, and commit messages. Non-American spelling is preserved only in compatibility-sensitive identifiers, schema/result keys,
filenames, literal quotations, and similar machine-facing values.

**D5. Phase 1 scope.** FDTD verification only: (A) vacuum, (B) homogeneous dielectric, (C) 1-D Bragg stack vs transfer-matrix
method, (D) finite 2-D bulk photonic-crystal stop band vs MPB, with r/a = 0.30. No W1 / slow-light device simulation until the
Phase 1 checkpoint is reviewed.

## 2026-10-10 — Phase 3 (dispersion engineering / delay-bandwidth trade study)

**D6. Phase 3 outcome: PASS WITH RESTRICTIONS; Phase 4 candidate = D0 (plain W1).** Engineered first rows (D1, D2) raise n_g by only 7% and 14% inside the flat region and cost 43% and 62% of the +-10% usable fractional bandwidth
(FDTD, validated range 0.231-0.268: D0 14.8%, D1 8.4%, D2 5.6%) with 4-8 times higher GDD; transmission and reflection inside the bands are not the problem (T >= 0.987, R <= 1%). See `docs/phase3_report.md`.
Restrictions: FDTD-validated claims are limited to 0.231-0.268 (the Phase 2 post hoc boundary, which also bounds every primary band); FDTD pointwise GDD is noise-limited; R is MPB only.
**D7. Target delay for the length metric:** 25 carrier cycles at each band center (criteria Revision 3, 2026-10-09, before the Phase 3 FDTD results).
**D8. Bug fixed:** `w1fdtd.hole_list` discarded the sign of s (negative s became outward); first D1/D2 runs discarded and re-run; regression test added (`tests/test_w1fdtd_geometry.py`).

## 2026-10-10 — Phase 4 (four-channel TTD family)

**D9. Phase 4 outcome: PASS WITH RESTRICTIONS; Phase 5 not yet authorized.** Four directly simulated D0 channels, N0 = 8, dN = 1 (declared before any Phase 4 FDTD), target step 3.9323 a/c, f0 = 0.2495. Common usable band 0.231-0.268 (14.8%, domain-limited by the post hoc Phase 2 validation limit); band-fit adjacent steps within 0.2% of target, amplitude imbalance <= 0.023 dB, T >= 0.997, raw phase intercept errors up to 2.68 rad, trimmed (one constant per channel at f0) <= 0.026 rad, cascade agreement <= 0.006 rad, res 30 vs 40 agreement 0.006 rad. See `docs/phase4_report.md`; handoff `results/phase4/phase4_handoff.npz`.
Phase 3 history is unchanged: D0 was selected because D1/D2 traded too much bandwidth for modest extra delay.
**D10. Method change (criteria Revision 5):** `stop_when_dft_decayed` alone gave erratic, under-converged stop times; all Phase 4 runs use a common minimum run time of 3500 after the source (first attempt discarded, kept as `raw_attempt1`).
**D11. Estimator change (criteria Revision 5, post hoc, disclosed):** the band criterion on adjacent increments uses the local-quadratic estimator; the pre-declared pointwise central difference is noise-limited (5.5-7.5% rms at one-period spacing) and is reported alongside.
