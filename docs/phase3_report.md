# Phase 3 report: dispersion engineering and delay-bandwidth trade study

Date: 2026-10-10. Baseline: Phase 2 checkpoint 38d4c82. Units: a = c = 1, f = a/lambda, delays in a/c. All numbers below are generated into
[phase3_report_tables.md](phase3_report_tables.md) (tables T1-T9) from `results/phase3/*.json`; figures are in `figures/phase3/`.
Processed results, scripts and tests are committed; raw plane-field dumps (`results/phase3/raw/`, 18 files) are git-ignored.

## 1. Objective

Does modifying the nearest rows of holes around the W1 defect produce a practically better delay-bandwidth tradeoff than the plain W1?
The goal is a quantified tradeoff (delay per length, usable bandwidth, ripple, dispersion, transmission, reflection, length, numerical robustness) for a
small predeclared set of designs, not maximum group index.

## 2. Starting Phase 2 validation basis

Plain W1 (triangular lattice of air holes, eps = 12, r/a = 0.30, TE, Hz-even MPB band 8), custom MPB-mode projection (not `get_eigenmode_coefficients`).
FDTD vs MPB n_g in the Phase 2 validated range 0.231-0.268 (14.8% fractional bandwidth): median bias +0.3 to +0.5%, RMS 0.4-0.7%, max <= about 1.4%; over the
full pre-declared interval 0.2312-0.2742 the maximum is 2.9%. Quantitative FDTD claims in this report rely on that pipeline and on decay_tol <= 1e-7 (tol = 1e-6 gave 5-14%
rms n_g noise and is not used). **The upper boundary 0.268 was chosen post hoc in Phase 2, after the data were seen, and is recorded as such in
`docs/phase2_5_criteria.md` Revision 2.** The restriction to 0.231-0.268 is a Phase-2-derived validation limit; it is **not** a Phase 3 predeclared criterion.
The source/monitor-distance and PML sensitivity checks of Phase 2 were done only at tol = 1e-6 and are not precision convergence evidence; they were not repeated.

## 3. Designs

D0 plain (f1 = 1.0, s = 0); D1 (f1 = 1.2, s = -0.05); D2 (f1 = 1.2, s = -0.10), the three predeclared Phase 0.5 H_r030 sweep points; f1 scales the radius of the rows |j| = 1,
s shifts them outward (negative = toward the defect), identically in MPB (`geometry.triangular_w1`) and FDTD (`w1fdtd.hole_list`; a regression test now enforces this, see 13).
Device: N = 4, 8, 12 periods of the engineered rows between plain-W1 leads (reference planes 3 periods outside the section); the interface mismatch is physical and is not removed.
Three designs are a **tradeoff illustration, not a Pareto frontier**.

## 4. Numerical methods

* MPB: 3a (res 32, nk = 101, M = 6), 3c (res 48 check), 3d (analysis). FDTD: `3b_fdtd.py` (Phase 2 runner; production res 30 for all nine device runs; res 40 confirmation for all nine; decay_tol = 1e-7; the runner refuses larger tolerances unless marked exploratory; none were). Frequencies 0.220-0.275, step 5e-4.
  Configuration, run times, stop times and provenance (git head, code_dirty flag, library versions, timestamps) are in T1 and in each `3e_*.json`/raw file.
* Modal S21, S11 by Lorentz-reciprocity projection onto the plain-W1 forward/backward modes (`guided.project`); de-embedding by the exact wave-cascade (T-matrix) method and by a linear fit of the S21 phase over the three lengths; the 2 pi branch is fixed with the MPB Bloch phase of the *device-section* band.
* Delay estimators: central difference of the unwrapped phase slope (primary, as in Phase 2) and a local quadratic (sensitivity).
* Criteria Revision 3 (dated 2026-10-09, before any Phase 3 FDTD result) fixed the 25-carrier-cycle target delay and the handling of the validated range.

## 5. ttdmetrics definitions

`src/pcttd/ttdmetrics.py` was recreated from scratch (the earlier zero-byte stub was not recovered); 10 unit tests (`tests/test_ttdmetrics.py`) cover fractional bandwidth
(f_high - f_low)/f_center, GDD = d tau / d omega (checked on linear and quadratic delays), median, ripple (max-min)/median, max relative deviation, contiguous runs,
usable-band search, GBP = fractional bandwidth x median n_g, required length, common band, NaN/mask handling and tolerance monotonicity.
Usable band = widest contiguous sample interval on which T >= 0.5 for every N, the device-section band exists, tau is finite, and every sample lies within +-tol of the interval's own median
(all windows tested, since the median changes with the window); ties are broken by smaller ripple. Primary tol = 10%; 5% and 25% are sensitivity only.
Target delay = 25 carrier cycles at the band center (Revision 3); length in periods and in free-space wavelengths at the center.

## 6. MPB results (MPB-only conclusions)

(T2, T8; `figures/phase3/3a_bands_ng.png`.) All three designs keep the same Hz-even band 8: it is monotone decreasing over the retained branch, its lower edge rises from 0.2134 (D0) to 0.2295 (D1) and 0.2368 (D2),
and n_g rises toward that edge faster for D2 than D1. Resolution/k-grid: 51-point spline derivative vs Hellmann-Feynman n_g differs <= 6e-4 (0.231-0.268); res 48 vs 32 n_g differs <= 7e-4 inside the MPB +-10% bands
(up to 13% for D2 at 0.237, i.e. at the band edge, and 2.6% for D1 at 0.231; those edge regions are not used). D1 and D2 have narrow gaps to even band 9 (0.0017, 0.0003 at the closest k); the band tracking is also corroborated by FDTD (Section 7).
Confinement (strip energy) was not recomputed for D1/D2.
MPB-only +-10% n_g-flatness bands in the validated range: D0 14.8% (n_g 3.91), D1 8.2% (4.19), D2 5.4% (4.47); unrestricted: D0 18.0%, D1 10.1%, D2 7.2%.

## 7. FDTD results (FDTD-validated conclusions, 0.231-0.268 unless stated)

(T3-T7; figures 3b-3d.) Primary usable band (+-10%, T >= 0.5 for N = 4, 8, 12, validated range):

| design | band | fractional BW | median n_g | ripple | GDD per period, band average (MPB center) | GBP | periods for 25 cycles (wavelengths) |
|---|---|---|---|---|---|---|---|
| D0 | 0.2310-0.2680 | 14.8% | 3.93 | 13.2% | -0.91 (-0.79) | 0.583 | 25.5 (6.36) |
| D1 | 0.2465-0.2680 | 8.4% | 4.20 | 14.6% | -3.86 (-3.40) | 0.351 | 23.1 (5.95) |
| D2 | 0.2535-0.2680 | 5.6% | 4.47 | 15.0% | -7.07 (-6.63) | 0.248 | 21.5 (5.59) |

(GDD units (a/c)^2 per period.) T_min >= 0.987 and R_max <= 0.010 inside every band; |1 - R - T| <= 0.008; projection residual median 0.025-0.030, max 0.055.
FDTD n_g vs MPB inside the bands: median +0.4% / +0.2% / +0.1%, max 3.3% / 2.7% / 1.9% (D0/D1/D2).
**All three primary bands end exactly at 0.268, so the validated-range limit, not the delay-flatness rule, bounds them on the high side.** Without the restriction (not validated above 0.268, T4)
the same rule gives 18.0% / 9.9% / 7.6% (D0/D1/D2); the ranking is unchanged. The 5% / 25% sensitivity bands are in T4 (validated: D0 12.9% / 14.8%, D1 5.8% / 10.8%, D2 3.8% / 7.6%); the 25% band of D0 saturates at the range edges.
The "(max - min)/median" ripple of 13-15% is consistent with the +-10% rule (which allows up to 20%) and is partly FDTD estimator noise: the MPB ripples over the same bands are 9.4%, 13.1%, 14.7%.

## 8. Usable-band construction

Delay flatness is the binding rule for D1 and D2 (their delay falls steeply from the band edge: n_g about 12.5 at 0.232 for D1 and 13 at 0.240 for D2), so their flat region starts well inside the transmitting band.
D0 is flat from the validated lower edge. Transmission is never the limiting condition inside the designs' bands: T >= 0.5 for all N is reached just above the device-section band edge (about 0.229 for D1 and 0.237 for D2; N = 4 turns on last), and below it the device is evanescent (T <= 5e-3 for N = 12, R -> 1).

## 9. Delay-bandwidth tradeoff

D1 and D2 buy 6.9% and 13.7% more delay per period (median n_g vs D0) and shorten the 25-cycle line by 9% and 16% (23.1, 21.5 vs 25.5 periods; 5.95, 5.59 vs 6.36 wavelengths), at the cost of 43% and 62% of the +-10% fractional bandwidth (8.4%, 5.6% vs 14.8%),
and of 40% and 57% of the group-delay-bandwidth product. D0 has the best GBP (0.58 vs 0.35, 0.25). Within this three-point set no engineered design improves delay-bandwidth performance; the plain W1 is better in bandwidth, GBP, GDD and tolerance to the band edge,
the engineered rows only trade bandwidth for a small length saving. (`figures/phase3/3g_tradeoff.png`; this is a three-point illustration, not a frontier, and a larger parameter search could change the picture.)

## 10. GDD / dispersion

GDD per period (a/c)^2: band-average slope of tau vs omega, FDTD -0.91 (D0), -3.86 (D1), -7.07 (D2); MPB values at the band center -0.79, -3.40, -6.63, i.e. the engineered designs are 4-8 times more dispersive; the sign is negative (delay decreases with frequency) for all.
The FDTD *local* GDD (second derivative of the measured phase) is noise-limited (|GDD| up to 9-16 inside the bands, and center values that differ from the band average by several units, T3); it is reported but should not be used quantitatively. MPB GDD is the reliable pointwise curve; the band-average FDTD slope agrees with it within about 10%.
The maximum |GDD| of the MPB curves over the bands is in `3d_mpb_analysis.json` / T8.

## 11. Interface mismatch

Inside the usable bands R <= 0.26% (D0), 0.73% (D1), 1.0% (D2), with Fabry-Perot-type oscillations of R between about 1e-4 and 1e-2 (`3b_transmission_reflection.png`); unitarity holds to <= 0.8%.
The larger interface reflection is concentrated at the engineered band edge (e.g. R = 1.6% at 0.232 for D1 and 2.3% at 0.240 for D2 at N = 12, rising to unity below the edge). Increased delay therefore does not cost transmission inside the flat band, but it does cost flatness and bandwidth because the useful region sits above the steep edge.
The cascade and the S21-phase estimators agree: median delay over the bands from the length pairs (4,8), (8,12), (4,12) and the three-length fit differ by <= 0.4% (T5); pointwise differences between pair estimators up to 7% (D0) reflect the short lever arm (dN = 4) of the noisy phase, not a bulk-delay change (the 4-12 pair differs from the fit by <= 1.1%).
No resonant phase slope was misread as bulk delay: the fit residual of the phase vs N is <= 0.005 rad over the bands.

## 12. Conventional R reference (MPB only; incomplete reference comparison)

R: square lattice of rods, TM, r/a = 0.19, one rod row removed, MPB res 24, nk = 61, M = 6 (Phase 0.5 geometry; no FDTD, no transmission or interface information). +-10% n_g band 0.3698-0.4193, fractional BW 12.5% (lower bound: the retained branch was cut at 0.42 by the 3a selection window),
median n_g 1.82, ripple 11.7%, GDD at the center -0.61, GBP 0.229, 34.7 periods (13.7 wavelengths) for the 25-cycle target (T9). Like-for-like against the D0 MPB band (14.8%, n_g 3.91, GDD -0.79, GBP 0.58, 25.7 periods, 6.4 wavelengths from MPB):
the H-family W1 needs about 2.1 times fewer wavelengths of length for the same delay at comparable fractional bandwidth. This is an MPB-only statement about the guided band; it says nothing about R's coupling, loss or disorder, nor does R share D0's FDTD validation.

## 13. Convergence and numerical limitations

* **Incident (found and fixed during Phase 3).** The first D1/D2 FDTD runs used a wrong geometry: `w1fdtd.hole_list` used `copysign(s, j)`, which discards the sign of s, so negative s shifted the rows *outward* instead of toward the defect (Phase 2 used s = 0 and is unaffected).
  The symptom was a D2 device that transmitted (T = 1) below its MPB band edge with the D0 delay and a constant 2.6 rad fit residual. The runs were discarded (not used anywhere), the code fixed, a regression test `tests/test_w1fdtd_geometry.py` added that compares the FDTD hole layout with the MPB geometry for four (f1, s) pairs, and D1/D2 re-run at res 30 and 40. D0 (s = 0, f1 = 1) is independent of the sign convention, so its runs were kept.
* Resolution: res 40 vs res 30 (T7): bands identical for D0 (14.8%), 8.0% vs 8.4% (D1) and 5.0% vs 5.6% (D2) (band-edge shifts of 0.001-0.0015 caused by ripple noise); median n_g within 0.25%; pointwise central-difference n_g differs up to 3-6% (estimator noise). The usable bands are not converged better than about +-0.5 percentage points of fractional bandwidth.
* Production tolerance 1e-7 for all runs; stop times are not monotone in N for D1 (t_end 1228, 2717, 1838) and D2 (5583, 3261, 2546) near the engineered band edge (long ringing); res 30 and res 40 agree despite that.
* The upper-edge anomaly of Phase 2 (above 0.268) is not investigated further; no Phase 3 claim uses the FDTD data above 0.268 as validated.
* Transmission > 1 beyond 1% appears only at f < 0.23 for D0 (max 1.024 at N = 4), outside the validated range.
* Not done: independent FDTD pulse check for D1/D2; J = 8 cell test; source/monitor/PML precision convergence; confinement for D1/D2; any robustness/disorder study.

## 14. Phase 3 gate: PASS WITH RESTRICTIONS

1. Modes identified: yes (continuous band 8, FDTD Bloch phase tracks MPB, transmission turns on at the MPB band edge).
2. FDTD extraction trustworthy: yes within 0.231-0.268 (tol 1e-7; res 30/40 agreement; unitarity <= 0.8%; projection residual <= 5.5%; delay vs MPB median <= 0.4%, max <= 3.3%).
3. Contiguous usable band under the declared T >= 0.5 and +-10% rules: yes for all three designs.
4. Delay stability across length: yes (band median agrees to 0.4% across pairs and with the fit).
5. Interface effects: understood for the usable bands (R <= 1%) and near the edges; the restriction is that they are measured for engineered-to-plain transitions of one geometry, not for channel families.
6. Justified Phase 4 candidate: D0 (below).
7. R represented honestly: MPB only, incomplete, lower-bound bandwidth.
8. Caveats documented.
Restrictions: FDTD claims are limited to 0.231-0.268 (post hoc Phase 2 boundary) and all primary bands are bounded by it; FDTD local GDD is noise-limited; resolution convergence of band edges is about half a percentage point of bandwidth; no pulse check or J = 8 test for D1/D2.
No criterion was weakened; the only additions are dated Revision 3 entries made before FDTD results.

## 15. Scientific interpretation

For this r/a = 0.30 W1, modifying the first rows as in D1/D2 raises n_g by only 7-14% inside the flat region, moves the band edge up, and makes the flat region narrower and more dispersive, with no transmission penalty inside the band;
the plain W1 already provides an approximately flat n_g of about 3.9 over about 15% bandwidth (and about 18% unrestricted per MPB/FDTD). Higher n_g (to ~13) exists only within a few percent of the engineered band edge where the delay is far from flat.
Relative to R (MPB only) the H-family W1 shortens the line by about 2x at comparable fractional bandwidth. These are statements for three points of a two-parameter family at one baseline radius.

## 16. Limitations

Three designs; one baseline radius; 2-D, lossless, infinite lateral continuation within the finite cell (J = 6); no disorder/fabrication or loss; the +-10% band is also bounded by the validated range;
the 25-cycle target is a stated normalization, not a system requirement; R compared by MPB only; GDD from FDTD is noisy.

## 17. What Phase 4 is allowed to consume

* Design **D0 (plain W1)**: four channels of N0 + m dN periods with the length set to be validated directly; per-period delay 3.93 a/c, flat (+-10% of median) over 0.231-0.268 (fractional bandwidth 14.8%, center 0.2495), T >= 0.99, R <= 0.3%.
* Use the FDTD-validated range 0.231-0.268 for all quantitative claims; treat 0.268-0.275 as unvalidated.
* Use decay_tol <= 1e-7 for production FDTD; the new `ttdmetrics` for bands and ripple; the cascade/length-fit machinery with the MPB branch reference.
* D1/D2 numbers may be used as a labeled tradeoff illustration. R may be cited as MPB-only.
* Phase 4 has not been started; it needs separate authorization.

### Phase 4 candidate selection

D0 is selected. It has the widest usable band (14.8% vs 8.4%, 5.6%), the best GBP, the lowest GDD (about 4-8 times lower), the flattest in-band delay and the lowest reflection, and its delay increment per period (3.93 a/c, about 0.98 carrier cycles at 0.2495) is the same one the Phase 4 family would use.
D1 and D2 save 9% and 16% of the length for the 25-cycle target, but at 43% and 62% bandwidth loss and much higher dispersion and band-edge sensitivity: they are not preferred for a four-channel family whose common band is the intersection of the channel bands.
