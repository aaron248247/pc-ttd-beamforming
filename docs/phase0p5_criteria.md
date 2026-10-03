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
