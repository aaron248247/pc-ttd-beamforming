# Phase 0.5 — Geometry selection: pre-declared criteria

Status: **declared before any Phase 0.5 computation was run** (see git history:
this file is committed before any result exists).

All frequencies are normalized: `w = omega*a/(2*pi*c) = a/lambda`; wavevectors `k = k_phys*a/(2*pi)`.
Group index `n_g = c/v_g = dk/dw` (normalized units, so it is scale-free).

## Candidates (same host permittivity for fairness)

| | Candidate H (leading) | Candidate R (comparison) |
|---|---|---|
| Lattice | triangular | square |
| Scatterer | air holes in eps=12 | eps=12 rods in air |
| Polarization | TE (Hz out of plane) | TM (Ez out of plane) |
| Defect | W1 (one row removed along Gamma-K) | one row of rods removed along Gamma-X |
| Supercell | height 5*sqrt(3) a (9 hole rows + 1 vacant) | height 9 a (8 rod rows + 1 vacant) |

Polarization convention is MPB/MEEP: TE = (Ex, Ey, Hz); TM = (Ez, Hx, Hy).

## Step 1 — Bulk gap and baseline radius (rule, applied identically to both)

1. Radius scan: H: r/a in [0.20, 0.40]; R: r/a in [0.10, 0.35], step 0.02, refined to 0.005 around the peak.
2. Gap between bands n and n+1 is `max_k w_n < min_k w_{n+1}` over the full-BZ k-grid
   (high-symmetry path **and** an interior grid; the larger-width conclusion is the one that holds
   on the denser sample, i.e. gap uses all sampled k).
3. Considered gaps: n = 1, 2, 3.
4. **Candidate gap** = lowest n whose maximum relative width `(w_hi - w_lo)/w_mid` over the scan is >= 5%.
   **Baseline radius** = radius maximizing that relative width.

## Step 2 — Defect sweep (same for both)

Parameters: first-neighbour-row radius scale `f1 = r1/r in {0.8, 0.9, 1.0, 1.1, 1.2}` and
first-neighbour-row outward lateral shift `s in {-0.10, -0.05, 0, +0.05, +0.10}` (units of a;
positive = rows move away from the defect). 25 designs per candidate, `f1=1, s=0` is the plain W1.

At each of 41 k-points on [0, 0.5], all supercell bands with frequency inside the bulk gap are defect bands.
Per band: out-of-plane-field parity (Phi = Hz for TE, Ez for TM) about the defect centre line,
group velocity from MPB (Hellmann-Feynman), cross-checked against finite differences.

Definitions:
* `N_fwd(w)` = number of crossings of the level `w` by all in-gap bands over k in [0, 0.5]
  (equals the number of forward-propagating guided modes).
* **SM-strict window**: contiguous interval with `N_fwd = 1` and `v_g > 0`.
* **SM-even window**: same, counting only even-Phi bands (odd modes may coexist).
* **Flat slow-light window** (inside an SM-strict window): contiguous sample interval where every `n_g`
  is within +-10% of the interval-mean `n_g`, with mean `n_g >= n_min`.
  `GBP = mean(n_g) * (w_hi - w_lo)/w_center`. Reported for `n_min` in {5, 10, 20}.
* **Strip energy fraction**: fraction of `|D|^2/eps + |H|^2` inside `|y| <= 1.5 * row_pitch`
  (row pitch = sqrt(3)/2 a for H, a for R), evaluated for the band at the centre of its window.
  **Edge energy fraction**: same within 1 row pitch of the supercell boundary (supercell-convergence guard).

## Gates (a candidate must pass all to be eligible)

* **G1** bulk relative gap >= 5% at the baseline radius.
* **G2** at least one sweep design has an SM-strict window with fractional width `dw/w_c >= 1%`
  whose frequencies are all even-Phi (excitable by a symmetric source).
* **G3** that window has `v_g > 0` throughout and edge energy fraction < 1e-3 (supercell converged).
* **G4** some design has max `n_g` >= 10 inside an SM-strict window, and max `n_g` over the sweep varies by
  a factor >= 2 between the best and worst designs that have a window (the parameters actually tune slow light).

## Ranked criteria (compared only between eligible candidates)

A candidate "wins" a criterion if better by more than 10% (relative); otherwise tie.

| | Criterion | Metric | Better |
|---|---|---|---|
| C1 | Gap width | baseline relative gap | larger |
| C2 | Single-mode bandwidth | max over sweep of SM-strict fractional width (even band) | larger |
| C3 | Confinement | median over sweep of strip energy fraction | larger |
| C4 | Slow-light controllability | max over sweep of GBP (n_min = 10) | larger |
| C5 | Tractability | FDTD pixels per waveguide period = `H_cell * res_req^2`, `res_req = 8 px / r_min` (smallest radius in the sweep) | smaller |
| C6 | FDTD suitability | fraction of sweep designs with an even SM-strict window >= 1% wide that lies >= 5% of the gap width away from both gap edges | larger |

## Decision rule

1. If exactly one candidate is eligible it is selected.
2. If both are eligible: most criteria won; tie broken by C4, then C2.
3. If neither is eligible: no geometry is selected, the failure is reported, and the design is revisited with review.
   The criteria are not edited after results are seen. Any later change is recorded as a new, dated revision.

## Known limitations of the study design (declared in advance)

* Baseline radius is bulk-gap optimal, not defect optimal; radius is not co-swept with the defect parameters.
* Only two defect parameters are swept. Other known slow-light levers (second-row shifts, etc.) are not explored.
* 2-D plane-wave results have no light line; real slab devices will be worse.
* Screening resolution is lower than the final check resolution; the winner is re-checked at higher resolution.

---

## Revision 1 (dated 2026-10-03, after inspecting one raw dispersion — disclosed)

**What triggered it.** While debugging the first end-to-end MPB run I inspected the raw in-gap bands of one
plain triangular-lattice W1 (r = 0.30, M = 4, coarse resolution). The even guided mode spans nearly the whole gap but has
`dw/dk < 0` on the folded zone `k in [0, 0.5]`. My declared "SM-strict ... and `v_g > 0`" rule, and gate G3's
"`v_g > 0` throughout", therefore classified the principal guided mode as invalid. That rule was a physics error, not a
property of either geometry.

**Why the original rule was wrong.** Bloch wavevectors are defined modulo `2*pi/a`. A band crossing at reduced
`k in [0, 0.5]` with negative slope is a positive-power-flow mode whose true Bloch wavevector lies in `(0.5, 1)`; its
time-reversed partner `(-k, -v_g)` always exists. So every crossing of a level over `[0, 0.5]` is exactly one forward-propagating
mode, regardless of the sign of `dw/dk` there. A backward-power-flow mode is not a distinct physical category in a periodic
waveguide at the single-band level.

**Changes (applied identically to H and R; no threshold or weight changed):**
* `N_fwd(w)` = number of crossings (unchanged), but the extra condition `v_g > 0` is **removed** from SM-strict / SM-even.
* `n_g = c/|v_g|` (absolute value).
* G3 now reads: window has a single crossing at every level (so `v_g` cannot change sign inside it; recorded as
  `vg_sign_constant`) **and** edge energy fraction < 1e-3.
* Everything else (G1, G2, G4, C1–C6 definitions, 10% margin, decision rule) is unchanged.

**Consistency note.** At the time of this revision I had seen: the bulk gap of each candidate at its first-guess radius
(H: r = 0.30, R: r = 0.20) and the raw H dispersion above. I had not seen any R defect dispersion nor any ranked metric.

---

## Revision 2 (dated 2026-10-03; made after the bulk-gap scan and **before any defect sweep was run**)

**Bulk results seen at this point (so the record is complete).** Triangular TE gap n=1 relative width increases
monotonically with r over the whole declared range and is still increasing at the limit r/a = 0.40
(0.490 at 0.40; the scan script's refinement had overshot to 0.42, 0.509, outside the declared range — an
implementation bug, not a rule change). Square TM gap n=1 peaks in the interior near r/a = 0.19 (0.391).

**Changes.**
1. *Clip clarification (bug fix, not a rule change).* The refinement around the peak is restricted to the declared scan
   range. For H this gives `r_base = 0.40` (the range limit; the unconstrained maximum lies at larger r).
2. *Resolution-validity exclusion (new, applies to both lattices, declared before the defect sweep).* A sweep design
   is **unresolved** at the screening resolution if its smallest edge-to-edge wall `min_edge_gap < 3/res`
   (= 0.125 a at res 24). Unresolved designs are recorded but not solved, and are excluded from every ranked
   metric and gate. This is a numerical-validity condition; it does not depend on any band-structure result.
3. *Supplementary sensitivity configuration "H_r030" (not decision-bearing).* Because the gap-optimal H radius sits on a range
   limit with thin (0.2 a) walls, the identical sweep and analysis are also run for triangular holes at the conventional
   r/a = 0.30, and the decision outcome is reported both with and without it. The decision rule (above) is evaluated on
   the declared configurations H (r_base from rule) and R. If H_r030 would change the outcome, that is reported as a
   finding, "decision sensitive to baseline radius", and the choice is referred to review.

---

## Revision 3 (dated 2026-10-03; made after the first full defect sweep — disclosed, no rule or threshold changed)

Defect results had been seen when these were made. None of them alters a gate threshold, criterion definition,
margin, or the decision rule. They are code-correctness fixes and *added* sensitivity analyses, listed so the record is complete.

**Code-correctness fixes found during the phase (all with regression tests):**
1. `geometry.py`: `math.copysign(s, j)` returned `|s|` with the sign of `j`, so every inward shift (`s < 0`) was built as the
   outward shift. The whole first sweep was discarded and re-run after the fix.
2. `bands.class_branches`: a frequency cut (`wcut`) plus minimum-count truncation silently dropped the even defect
   branch at small k for some designs. Replaced by an explicit sufficiency check (`branches_sufficient`); designs that
   fail the check are rejected by the selection script. H and H_r030 were re-solved with 20 extra bands (R passes with 8).
3. The scan script's refinement overshot the declared radius range (Revision 2, item 1).

**Added sensitivity analyses (reported next to, never instead of, the declared comparison):**
* `R_M6`: the R configuration with a wider supercell (M = 6), because R at the declared M = 4 failed G3 solely through the
  edge-energy criterion. R_M6 is eligible and its metrics match R to the digits reported.
* `H_r030`: as defined in Revision 2.
* Four pairings are evaluated (H|R, H|R_M6, H_r030|R, H_r030|R_M6). The **declared** decision is H vs R.
