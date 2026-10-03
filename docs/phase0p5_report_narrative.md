# Phase 0.5 report — 2-D geometry selection with MPB

Status: Phase 0.5 complete; **awaiting scientific review before Phase 1.** Criteria and decision rule:
[phase0p5_criteria.md](phase0p5_criteria.md) (committed before any result; Revisions 1–3 are dated and disclosed there).
Raw data: `results/phase0p5/`. Figures: `figures/phase0p5/`. Tables in the appendix are generated from the saved JSON.

## 1. Decision

**Selected family: triangular-lattice air holes in ε = 12, TE (Hz), W1 defect (candidate H).**
The declared comparison (H vs R) selects H, and so do all three sensitivity pairings (appendix T4).

**This is a moderate-margin decision, not a landslide.** Read the reasons before accepting it:

* In the declared comparison R is ineligible *only* because it fails gate G3 (supercell edge energy 1.4×10⁻³ – 3.2×10⁻³ against a
  1×10⁻³ limit at M = 4). That is a supercell-width effect, not a property of the structure: at M = 6 every R design
  passes (≤ 1.1×10⁻⁴) and the physics metrics are unchanged to the digits reported (T3, "R_M6"). With R made eligible, H still
  wins on count (C1, C4, C5 vs C2; C3 and C6 tie), but see the next two points.
* H wins C1 only because the declared rule maximizes the gap and the maximum lies at/beyond the range limit r/a = 0.40
  (thin 0.2 a walls). At the conventional r/a = 0.30 (config H_r030) R wins C1, C2 *and* (R_M6 pairing) the count is 2–2,
  decided only by the pre-declared tie-break on C4 (slow-light GBP).
* R wins C2 (single-mode bandwidth) decisively: 0.385 vs 0.154 (H) and 0.138 (H_r030).

So the choice of H rests on **C4 (flat slow-light potential) and C5 (FDTD cost)**; R's strength is a wide, clean single-mode band.
Given the project's purpose (slow-light delay) that is the right tie-breaker, but it is a judgement the rule encodes, and I flag it.

## 2. Why H (what the data show physically)

* **Candidate R is, physically, an air-channel waveguide.** Removing one rod row gives one even TM band with a fractional bandwidth of 0.337
  (plain W1, T5; the gap is 0.391 wide) with a cutoff at k = 0, group index ≈ 2.03 at the window centre, and a large n_g only very close to cutoff
  and near the zone edge. There is no extended flat slow-light region: its best n_g ≥ 10 flat window has n_g ≈ 10 over
  0.3% fractional bandwidth. It is, in effect, a strong candidate for the **conventional-waveguide baseline** in the
  Phase 3/4 comparison, which is a useful by-product.
* **Candidate H has the controllable slow-light region.** Dispersion engineering with the first-row radius and shift moves
  the flat-n_g windows substantially (appendix T3/T5; figure f4). The best flat window found (config H_r030, f1 = 0.9, s = +0.10)
  has n_g ≈ 19.6 over 1.0% fractional bandwidth (GBP 0.195).
* **H is cheaper for FDTD** by the declared pixels-per-period proxy (T3, C5): 2.6× (H_r030) to 6.7× (H vs R_M6) fewer pixels per period, because its smallest feature is larger.
* **H is multimode over its wide gap** (up to 7 modes at r/a = 0.40 for the best-GBP design): the usable windows are
  separate even-mode single-mode intervals, each 10–15% wide, not one broad band. For FDTD this means the source spectrum
  must be confined to one window.

## 3. An early warning for the research question (not a conclusion)

In this 25-point, two-parameter screening sweep the flat slow-light windows with n_g ≥ 10 are **about 1% fractional
bandwidth or narrower in every configuration** (T5), and n_g ≈ 4 elsewhere in the H single-mode windows. A broadband array
band of several percent or more will not fit inside a flat n_g ≥ 10 window. This does not yet say TTD fails — flatness of
n_g is one measure; Phase 2/3 will use the actual transfer function, GDD and the array metric — but it makes risk 1 of the
Phase 0 plan (bandwidth too narrow) the leading threat, and the conventional-waveguide baseline correspondingly important.

## 4. Numerical verification

| Check | Result |
|---|---|
| Unit/regression tests | 30 pass (`results/phase0p5/pytest_output.txt`): 24 analytic/non-MPB (band analysis, geometry, regressions) and 6 MPB-vs-closed-form |
| MPB vs closed forms | n_g = √12 (homogeneous ε = 12) and = 1 (vacuum) to rtol 1e-4; degenerate cos/sin parity pair detected; strip/edge energy fraction on a uniform wave matches area fraction to the grid quantum |
| MPB v_g vs finite difference | median relative difference ≤ 0.4% across the converged cases (T6) |
| Resolution 24 → 48 | band-gap edges move ≤ 0.2%; n_g at window centre unchanged to the digits shown (T6) |
| Solver tolerance 1e-7 → 1e-9 | no change in any reported quantity (T6) |
| Supercell M → M+1 | gap and n_g unchanged; **window-threshold metrics move** (below) |
| Band sufficiency | enforced (`branches_sufficient`); selection aborts if any design fails |

**What is and is not numerically stable.** Gap edges and group index along a band are well converged. Metrics that depend
on a *threshold on n_g* are not: GBP(n_g ≥ 10) changes by up to −32% (R plain W1, res 48) and −14% (H_r030 best, M+1), and
the single-mode fraction by up to −14% (H best, M+1). Differences below ≈ 30% in GBP, and below ≈ 15% in SM fraction, are therefore
**not significant**. The decision-relevant C2 (≈ 2.5× gap between R and H) and the C4 H-vs-R margins are larger than that;
the C4 H(r=0.40)-vs-R margin (0.048 vs 0.031) is only ~50% and should be treated as moderate.

Supercell width: for H and H_r030 not all designs satisfy the edge-energy limit at M = 5 (T7: 11 of 14 and 12 of 25).
G3 needs only one passing design per candidate, so gating is unaffected, but those *specific* designs' windows are not converged
and Phase 1 will use a wider supercell/cladding for any H design it adopts.

## 5. Defects found and fixed during this phase (honest log)

1. **Geometry sign bug** (`copysign`): every inward row-shift was built as an outward shift; ±s sweep results were identical.
   Caught because the paired results agreed to 10 digits. Whole sweep re-run. Regression tests added (tests only covered s > 0 before).
2. **Branch-tracking bug**: a frequency cut dropped the even defect branch for some plain W1s, so the plain H_r030 W1 initially
   showed no single-mode window. Fixed; sufficiency check added; H/H_r030 re-solved with more bands.
3. **Pre-declared rule error**: "v_g > 0 on k ∈ [0, 0.5]" mis-classified the principal guided mode (folded zone). Corrected as
   Revision 1 with the physical argument; applied identically to both candidates.
4. **Script overshoot** of the declared radius range (r = 0.42 instead of the range limit 0.40).
5. Two test tolerances were set from the wrong quantity (a pixel-quantized energy fraction; a hand-derived bound) and were
   corrected from the physical/analytic reasoning, not loosened to pass; a vacuous assertion I drafted was replaced by exact expected outcomes.

## 6. What is safe to claim now

* In 2-D with ε = 12, a triangular air-hole lattice has a complete TE gap whose relative width grows monotonically with r/a to
  the range limit 0.40 (0.490), and r/a = 0.30 gives 0.280; a square rod lattice has a TM gap peaking at r/a ≈ 0.19 (0.391).
* The plain W1 of each lattice has an even-parity single-mode window inside the gap: fractional width 0.106 (H), 0.107 (H_r030), 0.337 (R).
* The R W1 is a nearly-vacuum-like air waveguide (n_g ≈ 2 over most of its band); H W1 single-mode windows have n_g ≈ 4,
  with larger n_g only in narrow bands near band extrema.
* Flat n_g ≥ 10 windows exist in all candidates but are ≤ ~1% wide in this sweep.

## 7. Not claimed / limitations

* Nothing yet about delay, transmission, reflection, GDD of a *finite* device, or squint. Those need FDTD (Phase 1+).
* Two defect parameters, 25 points, one baseline radius per candidate: not an optimization; the flat-window numbers are not upper bounds.
* 2-D, lossless, non-dispersive, no light line: any slow-light figure here is optimistic for a real slab.
* `max n_g` values in the appendix are singular at band extrema (they depend on the level grid) and are not used in any criterion except the G4 threshold.
* Bulk gaps are evaluated on a 14×14 interior grid plus the high-symmetry path; a true minimum between samples is possible but unlikely.
* Windows use a 1500-level grid and a 41-point k-sampling with Hermite interpolation using MPB group velocities.

## 8. Decision requested for Phase 1

1. Confirm **triangular-lattice air-hole, TE, W1** as the working geometry.
2. Choose the baseline radius: **r/a = 0.40** (the declared rule's output; 49% gap, walls only 0.2 a, 11/25 sweep designs unresolved at res 24)
   or **r/a = 0.30** (outside the rule; 28% gap, all 25 designs resolved, richer slow-light tuning: GBP 0.195 vs 0.048). I lean to 0.30 for
   Phase 1–3 because slow-light controllability is the project's point, but that departs from the pre-declared rule, so it is your call.
3. Whether to carry candidate R's air-channel waveguide forward as the Phase 3/4 **conventional-waveguide baseline**.
