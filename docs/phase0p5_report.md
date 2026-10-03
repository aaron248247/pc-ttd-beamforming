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
Given the project's purpose (slow-light delay) that is the right tie-breaker, but it is a judgment the rule encodes, and I flag it.

## 2. Why H (what the data show physically)

* **Candidate R is, physically, an air-channel waveguide.** Removing one rod row gives one even TM band with a fractional bandwidth of 0.337
  (plain W1, T5; the gap is 0.391 wide) with a cutoff at k = 0, group index ≈ 2.03 at the window center, and a large n_g only very close to cutoff
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
| Resolution 24 → 48 | band-gap edges move ≤ 0.2%; n_g at window center unchanged to the digits shown (T6) |
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


---
# Appendix: generated tables

All values below are read from `results/phase0p5/*.json` by `experiments/phase0p5/06_report.py`.

## T1. Bulk gap (declared rule, res 32)

| lattice | pol. | gap n|n+1 | r_base/a | w_lo | w_hi | rel. width |
|---|---|---|---|---|---|---|
| triangular | TE (Hz) | 1 | 0.400 | 0.246 | 0.405 | 0.490 |
| square | TM (Ez) | 1 | 0.190 | 0.289 | 0.430 | 0.391 |


Relative gap vs radius, triangular (n=1): r=0.200:0.065, r=0.240:0.145, r=0.280:0.233, r=0.320:0.326, r=0.360:0.417, r=0.385:0.467, r=0.395:0.483


## T2. Gates and sweep coverage

| config | lattice | r/a | M | extra bands | solved | unresolved (excluded) | G1 | G2 | G3 | G4 |
|---|---|---|---|---|---|---|---|---|---|---|
| H | triangular | 0.400 | 5 | 20 | 14 | 11 | PASS | PASS | PASS | PASS |
| R | square | 0.190 | 4 | 8 | 25 | 0 | PASS | PASS | FAIL | PASS |
| H_r030 | triangular | 0.300 | 5 | 20 | 25 | 0 | PASS | PASS | PASS | PASS |
| R_M6 | square | 0.190 | 6 | 8 | 25 | 0 | PASS | PASS | PASS | PASS |


## T3. Criterion metrics

| config | C1 gap | C2 SM frac | C3 strip | C4 GBP(n≥10) | GBP(n≥5) | GBP(n≥20) | C5 px/period | C6 frac | max n_g (singular at band extrema) |
|---|---|---|---|---|---|---|---|---|---|
| H | 0.490 | 0.154 | 0.9980 | 0.0483 | 0.1583 | 0.0280 | 5413 | 1.00 | 489 |
| R | 0.391 | 0.385 | 0.9646 | 0.0313 | 0.0589 | 0.0000 | 24931 | 1.00 | 393 |
| H_r030 | 0.280 | 0.138 | 0.9930 | 0.1945 | 0.1945 | 0.1373 | 9623 | 1.00 | 700 |
| R_M6 | 0.391 | 0.385 | 0.9640 | 0.0312 | 0.0583 | 0.0000 | 36011 | 1.00 | 181 |


## T4. Decision per pairing (pre-declared rule)

| pairing | eligible | per-criterion winners | selected | rule |
|---|---|---|---|---|
| declared: H vs R | {'H': True, 'R': False} | C1:H, C2:R, C3:tie, C4:H, C5:H, C6:tie | H | only eligible candidate |
| H vs R_M6 (R with converged supercell) | {'H': True, 'R_M6': True} | C1:H, C2:R_M6, C3:tie, C4:H, C5:H, C6:tie | H | most criteria won |
| H_r030 vs R (suppl. radius) | {'H_r030': True, 'R': False} | C1:R, C2:R, C3:tie, C4:H_r030, C5:H_r030, C6:tie | H_r030 | only eligible candidate |
| H_r030 vs R_M6 | {'H_r030': True, 'R_M6': True} | C1:R_M6, C2:R_M6, C3:tie, C4:H_r030, C5:H_r030, C6:tie | H_r030 | tie on count; tie-break C4 |


All pairings select the same family: **True** (H).


## T5. Plain W1 and best-GBP design per configuration

| config | design | f1/s | max modes in gap | SM-strict(even) frac | SM even-only frac | flat-window n_g | its frac. BW | its GBP | strip frac | edge frac (window max) | n_g at window center [min–max] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H | plain W1 | 1/+0.00 | 3 | 0.106 | 0.237 | — | — | — | 0.9995 | 2.1e-05 | 3.92 [3.92–21] |
| H | best GBP(≥10) | 0.8/-0.10 | 7 | 0.132 | 0.254 | 10.2 | 0.005 | 0.048 | 0.9967 | 2.1e-04 | 4.06 [3.97–489] |
| R | plain W1 | 1/+0.00 | 1 | 0.337 | 0.337 | 10.7 | 0.003 | 0.029 | 0.9727 | 1.8e-03 | 2.03 [1.79–56] |
| R | best GBP(≥10) | 0.9/+0.10 | 2 | 0.361 | 0.361 | 10.1 | 0.003 | 0.031 | 0.9695 | 2.9e-03 | 1.90 [1.58–98] |
| H_r030 | plain W1 | 1/+0.00 | 3 | 0.107 | 0.250 | 63.4 | 0.001 | 0.040 | 0.9978 | 1.1e-03 | 3.87 [3.87–4] |
| H_r030 | best GBP(≥10) | 0.9/+0.10 | 3 | 0.114 | 0.239 | 19.6 | 0.010 | 0.195 | 0.9969 | 3.1e-04 | 3.76 [3.76–4] |
| R_M6 | plain W1 | 1/+0.00 | 1 | 0.336 | 0.336 | 10.8 | 0.002 | 0.026 | 0.9720 | 4.6e-05 | 2.03 [1.79–56] |
| R_M6 | best GBP(≥10) | 0.9/+0.10 | 2 | 0.361 | 0.361 | 10.0 | 0.003 | 0.031 | 0.9686 | 8.7e-05 | 1.90 [1.58–70] |


## T6. Numerical convergence (same analysis code)

| config | f1/s | case vs base | bands ok | Δ gap_lo | Δ gap_hi | SM-even frac base → case | GBP(≥10) base → case | n_g at window center | median |v_g,MPB−FD|/v_g |
|---|---|---|---|---|---|---|---|---|---|
| H | 0.8/-0.10 | res48_tol1e-8 | True | -0.20% | +0.06% | 0.132 → 0.132 | 0.0483 → 0.0447 | 4.06 → 4.06 | 0.004 |
| H | 0.8/-0.10 | M+1 | True | -0.00% | -0.00% | 0.132 → 0.113 | 0.0483 → 0.0476 | 4.06 → 4.02 | 0.004 |
| H | 0.8/-0.10 | tol1e-9 | True | -0.00% | -0.00% | 0.132 → 0.132 | 0.0483 → 0.0483 | 4.06 → 4.06 | 0.004 |
| H | 1/+0.00 | res48_tol1e-8 | True | -0.20% | +0.06% | 0.106 → 0.107 | 0.0000 → 0.0000 | 3.92 → 3.92 | 0.003 |
| H | 1/+0.00 | M+1 | True | -0.00% | +0.00% | 0.106 → 0.106 | 0.0000 → 0.0000 | 3.92 → 3.92 | 0.003 |
| H | 1/+0.00 | tol1e-9 | True | -0.00% | -0.00% | 0.106 → 0.106 | 0.0000 → 0.0000 | 3.92 → 3.92 | 0.003 |
| H_r030 | 0.9/+0.10 | res48_tol1e-8 | True | -0.06% | -0.03% | 0.114 → 0.114 | 0.1945 → 0.1949 | 3.76 → 3.76 | 0.004 |
| H_r030 | 0.9/+0.10 | M+1 | True | +0.00% | +0.00% | 0.114 → 0.114 | 0.1945 → 0.1682 | 3.76 → 3.76 | 0.003 |
| H_r030 | 0.9/+0.10 | tol1e-9 | True | +0.00% | +0.00% | 0.114 → 0.114 | 0.1945 → 0.1945 | 3.76 → 3.76 | 0.003 |
| H_r030 | 1/+0.00 | res48_tol1e-8 | True | -0.06% | -0.03% | 0.107 → 0.108 | 0.0399 → 0.0398 | 3.87 → 3.87 | 0.003 |
| H_r030 | 1/+0.00 | M+1 | True | +0.00% | -0.00% | 0.107 → 0.109 | 0.0399 → 0.0393 | 3.87 → 3.87 | 0.003 |
| H_r030 | 1/+0.00 | tol1e-9 | True | -0.00% | -0.00% | 0.107 → 0.107 | 0.0399 → 0.0399 | 3.87 → 3.87 | 0.003 |
| R | 0.9/+0.10 | res48_tol1e-8 | True | -0.05% | +0.02% | 0.361 → 0.362 | 0.0313 → 0.0308 | 1.90 → 1.90 | 0.000 |
| R | 0.9/+0.10 | M+1 | True | +0.00% | +0.00% | 0.361 → 0.361 | 0.0313 → 0.0312 | 1.90 → 1.90 | 0.001 |
| R | 0.9/+0.10 | tol1e-9 | True | +0.00% | +0.00% | 0.361 → 0.361 | 0.0313 → 0.0313 | 1.90 → 1.90 | 0.000 |
| R | 1/+0.00 | res48_tol1e-8 | True | -0.05% | +0.02% | 0.337 → 0.337 | 0.0293 → 0.0199 | 2.03 → 2.03 | 0.000 |
| R | 1/+0.00 | M+1 | True | +0.00% | +0.00% | 0.337 → 0.336 | 0.0293 → 0.0262 | 2.03 → 2.03 | 0.000 |
| R | 1/+0.00 | tol1e-9 | True | -0.00% | +0.00% | 0.337 → 0.337 | 0.0293 → 0.0293 | 2.03 → 2.03 | 0.000 |


## Environment

```
{
 "python": "3.10.19 (main, Oct 21 2025, 16:43:05) [GCC 11.2.0]",
 "platform": "Linux-6.18.40.1-microsoft-standard-WSL2-x86_64-with-glibc2.35",
 "numpy": "1.24.3",
 "scipy": "1.11.1",
 "matplotlib": "3.10.6",
 "pandas": "2.0.3",
 "pytest": "7.4.0",
 "meep": "1.33.0-beta",
 "meep_with_mpi": true,
 "mpb_available": true,
 "git_head": "7783a200de84b17cd4c0913fc84f1a2f5033f73c"
}
```


Sweep provenance: `{"timestamp_utc": "2026-10-03T18:51:37+00:00", "git_head": "b06715e7cf8b227f79c62b50b8c50b8d38c1028b", "git_dirty": true, "meep": "1.33.0-beta", "numpy": "1.24.3", "scipy": "1.11.1"}`


## T7. Supercell edge-energy fraction over all solved designs (window max)

| config | M | min | max | designs with edge < 1e-3 | designs with window |
|---|---|---|---|---|---|
| R | 4 | 1.4e-03 | 3.2e-03 | 0 | 25 |
| R_M6 | 6 | 2.6e-05 | 1.1e-04 | 25 | 25 |
| H | 5 | 2.1e-05 | 6.2e-02 | 11 | 14 |
| H_r030 | 5 | 2.3e-04 | 1.8e-02 | 12 | 25 |
