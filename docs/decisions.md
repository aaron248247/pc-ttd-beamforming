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
