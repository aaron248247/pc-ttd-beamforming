# Phase 6 report — delay-granularity mitigation and robustness

Criteria: `docs/phase2_5_criteria.md` Revision 7 (commit 65f1ab3, declared before any Phase 6 result). Tables: `docs/phase6_report_tables.md` (Tables A–K). Handoff: `results/phase6/phase6_handoff.{npz,json}`. Figures: `figures/phase6/6_01 … 6_17`.
Starting checkpoint 959e4c9 (Phase 5, 173 tests passing). Phase 4 and Phase 5 result files are unchanged (sha256-checked in `tests/test_phase6.py`).

**Gate: PASS WITH RESTRICTIONS.** Phase 7 has not been started.

## 1. What was computed

* Hierarchical array: four PhC channels = four subarrays of Msub λ0/2-spaced elements, W(q,p) = G_q(f)·exp(+i2πf0 p d_e sinθ0), element (q,p) at (qMsub+p)λ0/2, full element-level far field (exact factorization, not point-source subarrays; unit-tested against `array.array_factor` and `beam.analyze_case`). Library: `src/pcttd/subarray.py`, `src/pcttd/robustness.py`.
* Msub = 1 (hybrid, Case C), 2, 4, 8; cases fixed phase / ideal element-level TTD / ideal subarray TTD / raw PhC / trimmed PhC / Phase 4 `H_trim` as is. Trimmed = constant phase per channel at f0 chosen for the **required** step (`trim_phase`); it changes the intercept only. "Phase 4 H_trim as is" is trimmed for 3.9323 a, not for the required step.
* Part B on the best architecture, system-level only (no new FDTD): deterministic perturbation sweeps with first-order analytic predictions, Monte Carlo (seed 20261010, 1000 samples per tolerance point, uniform errors on all channels), tolerance budget, res 40 check, period-quantization scan.
* Minor deviations from Revision 7, disclosed: (i) production grid is 40001 points in sin θ for pointing and lobes (finer than the 20001 declared for lobes; same method); (ii) post hoc additions: a measured *replica-lobe* level at the analytic subarray grating positions, an *efficiency loss at fixed input power* column, a *calibrated absolute-length* delay variant, and the Msub = 1…12 quantization scan; (iii) peak-power loss is reported as a positive loss so P95 is the bad tail (the first Monte Carlo run reported it with the opposite sign and was discarded; nothing else changed).

## 2. Delay requirement and match (Table A)

| Msub | elements | required step | PhC/required (ρ) | mismatch | f0·(τ_PhC − τ_req) |
|---|---|---|---|---|---|
| 1 | 4 | 0.25 cycle | 3.924 | +292 % | +0.731 cycle |
| 2 | 8 | 0.50 cycle | 1.962 | +96.2 % | +0.481 cycle |
| 4 | 16 | 1.00 cycle | **0.9811** | **−1.89 %** | **−0.0189 cycle** |
| 8 | 32 | 2.00 cycles | 0.4906 | −50.9 % | −1.019 cycle |

Case B (ideal reference only): required step λ0/4 = 1.0020 a = 0.25 cycle; **τ_PhC/τ_half = 3.9244**. Msub = 4 aligns the PhC one-period quantum with the required subarray step to 1.9 %; the other sizes do not (a PhC with one-period increments cannot supply them).

## 3. Beam results (res 30, mask_common; Tables B–D, Figures 6_03–6_10)

Analytic law for a retrimmed channel family, checked numerically: sinθ − sinθ0 = (f0/f − 1) sinθ0 (1 − ρ), i.e. the subarray-level squint is (1 − ρ) times the fixed-phase squint. A constant phase cannot change group delay, so a mismatched delay step leaves a squint that trimming cannot remove.

| Msub | trimmed PhC max abs error | R_squint | strongest other lobe | verdict |
|---|---|---|---|---|
| 1 (hybrid, element-level + phase) | 7.67° | 2.86 (worse than fixed phase, opposite sign) | −4.8 dB (ambiguity) | fails; analytic prediction 7.5° |
| 2 | 2.37° | 0.88 | −12.2 dB | no squint suppression |
| **4** | **0.189°** | **0.070** | **−12.4 dB** | **meets tracking, squint, lobe criteria** |
| 8 | 1.42° | 0.53 | −12.4 dB | squint only halved |

* **Best architecture: Msub = 4, phase-intercept-trimmed (re-trimmed) PhC subarray TTD** (the only qualifying case; selected by the declared rule). Tracking tolerance tol_pt = 0.05·HPBW_ideal = 0.367°; max error 0.189° (RMS 0.110°; low/f0/high = +0.157/0/−0.189°); offset-removed R_squint = 0.070 (stretch target ≤ 0.20 met).
* Hypothesis test: the Msub = 4 alignment **is confirmed numerically** (|ρ − 1| = 1.9 % ≤ 5 %, and the retrimmed case meets the criteria). The Phase 5 "too coarse" delay quantum is nearly ideal for four-element subarrays. The 1.9 % mismatch is **not rounded away**: it contributes a slope residual of ±0.05° over the band (analytic ±0.05°; numerical difference to ideal subarray TTD 0.055°). If only the Phase 4 trim is used (not re-trimmed to 4.008 a), the intercept error costs −0.59° at f0 (analytic −0.62°) and 0.73° max over the band (tol_pt 0.37° exceeded). Other sizes are poor: with k = 1 period per step, Msub = 3 and 5 give R_squint ≈ 0.25 (Table K), but nothing matches as well as 4.
* Even ideal inter-subarray TTD has 0.156° max error here: the fixed intra-subarray phases squint and pull the lobe. The PhC adds only 0.055° on top.
* Peak power: min 0.983 (−0.07 dB) for trimmed Msub = 4. HPBW at f0: 7.35° (Phase 5: 7.69°); aperture 8 λ0 (Phase 5: 7.85 λ0 pitch-spanned).
* Raw PhC (Msub = 4): constant offset ≈ −6.7° with strongest lobe −6.7 dB, but offset-removed spread only 0.12 of fixed phase: **raw squint is suppressed, absolute pointing is not**; constant phase trimming is therefore still required.
* Amplitude: |H| = 1 diagnostics change nothing visible (trimmed: max error and peak power equal to within 1e-4 / 0.002); amplitude mismatch does not matter at the Phase 4 level (≤ 1 dB imbalance).

## 4. Grating lobes (Table B/C/D, Figures 6_05–6_07)

* λ0/2 element spacing: zero visible element-level spatial-alias replicas over the band (ideal element-level TTD: strongest other lobe −13.15 dB, the uniform-array sidelobe). Verified, not assumed.
* The Phase 5 equal-strength grating lobes (0 dB) are gone: the trimmed Msub = 4 array has **one** lobe within 3 dB at every frequency, strongest other lobe −12.38 dB (ordinary sidelobe; ideal subarray TTD −12.39 dB) and strongest *subarray replica* (quantization) lobe **−16.6 dB** (worst over the band, at the band edges; nulled at f0 by the element factor).
* Not removed completely: the subarray pitch is 2 λ0, so replica lobes exist (4 visible at f0); they are suppressed by the λ0/2 element factor, which squints with frequency, so they grow from a null at f0 to about −17 dB at the band edges. The mitigation is from ≈ 0 dB to ≤ −16.6 dB, not elimination.

## 5. Robustness (Part B; Tables E–I, Figures 6_12–6_17)

System-level, on the Phase 4 `H_q(f)`; nominal max error 0.189°. Uniform independent errors on all channels. Combined L2 test (2°, 0.5 %, 0.25 dB): P95 max error **0.352°** ≤ 1.0° (passed). Monte Carlo is reproducible (tested).

* Phase-trim error ±a: P95 max error 0.385° at 5°, 0.60° at 10°. First-order prediction (Δsinθ = g/(2πf Msub d_e), g = LS slope over channels) agrees with the deterministic response to ~5 % (it slightly overestimates; the element factor stiffens the peak). Pointing shifts only; peak power moves < 0.05 dB.
* Delay error (e of q·τ_PhC): uncalibrated error produces a **frequency-independent** pointing shift (analytic and numeric: 0.31° for e = 1 % deterministic ramp); P95 max error 0.45° at 1 %, 0.74° at 2 %, 1.58° at 5 %. Calibrated (intercept re-trimmed, pure slope error) is far milder: P95 0.357° even at 10 % (shift ∝ (f−f0)/f; analytic ≈ numeric).
* Absolute-length reading (error relative to the whole (N+6)τ_PhC ≈ 55–67 a/c device delay): uncalibrated, P95 0.5° is reached at only 0.14 %, i.e. about 0.02 carrier cycle of absolute delay; calibrated at f0: 2.2 %. A length error is only benign if the constant phase is re-trimmed after fabrication.
* Amplitude error ±A: P95 max error 0.236° at 6 dB; fixed-input-power efficiency loss < 0.03 dB at 1 dB; strongest lobe degrades from −12.4 to −7.0 dB P95 at 6 dB. Amplitude mismatch is not a pointing limiter in this range.
* **System-level tolerance estimates** (P95 of max error including nominal, limit 0.5° / 1.0°; not fabrication specifications): phase ±7.6° / ±20°; delay (uncalibrated, of q·τ_PhC) ±1.2 % / ±2.9 %; delay (calibrated, of q·τ_PhC) > ±10 % (beyond the sweep) for both; delay (uncalibrated, absolute) ±0.14 % / ±0.36 %; delay (calibrated, absolute) ±2.2 % / ±5.3 %; amplitude > ±6 dB for both (beyond the sweep). Median-based limits are in Table I.
* Combined P95 max error: L1 0.26°, L2 0.35°, L3 0.56°, L4 0.92°.

## 6. Resolution check (Table J, Figure 6_17)

Msub = 4 trimmed PhC: max |θ30 − θ40| = 0.007° (limit tol_pt/3 = 0.122°); change of max error 0.0013°; change of R_squint 0.0005; change of peak power 0.0019; change of strongest lobe 0.041 dB (limit 0.5). All criteria pass; best architecture and gate unchanged. (Raw Msub = 1 differs by 0.29° but is not the best architecture.)

## 7. Gate questions

1. **Best use of the PhC delay quantum:** Msub = 4 (λ0/2 elements, one PhC channel per four elements).
2. **Msub = 4 alignment:** yes — required 1.000 cycle vs 0.981 (−1.89 %).
3. **Phase 5 equal-strength grating lobes removed:** replaced by one design lobe and ordinary sidelobes (−12.4 dB), with residual subarray replica lobes ≤ −16.6 dB; mitigated, not eliminated.
4. **Broadband squint suppression preserved:** yes; R_squint = 0.070, max error 0.189° vs 2.69° for fixed phase.
5. **Pointing error from 0.981 vs 1.000:** slope residual ±0.05° (band), 0.055° vs ideal subarray TTD; −0.59° at f0 if the Phase 4 trim is used unchanged.
6. **Constant phase trimming still required:** yes (raw offset 6.7°; Phase 4 trim as is 0.73°).
7. **Amplitude mismatch:** no (≤ 0.001° in pointing at 1 dB; ≤ 0.24° P95 at 6 dB).
8. **Phase-trim sensitivity:** P95 max error 0.39° at ±5°, budget ±7.6° (0.5°) and ±20° (1.0°).
9. **Delay sensitivity:** ±1.2 %/±2.9 % of q·τ_PhC uncalibrated; strongly dependent on whether the intercept is re-trimmed; absolute uncalibrated length error must stay below ≈ 0.14 %.
10. **Tolerance ranges:** Section 5.
11. **res 30 vs res 40:** stable (all declared criteria pass).
12. **Interpretation:** the PhC channels are best interpreted as **subarray-level TTD with a constant-phase trim per channel** (a hybrid at the subarray level: PhC supplies the delay slope, a constant supplies the intercept, fixed phase handles the intra-subarray steering). Element-level use (Msub = 1) fails: a constant phase cannot supply the missing 0.75-cycle of group delay and the squint is worse than fixed phase.
13. **Phase 6 improves the Phase 5 interpretation:** yes materially — the coarse quantum is a match, not a defect, for four-element subarrays — with the restrictions below.

## 8. Physical mapping (illustrative, normalized units retained)

| carrier | a | λ0 | element pitch | subarray pitch (Msub = 4) | aperture (16 el.) | delay step (PhC / required) |
|---|---|---|---|---|---|---|
| 30 GHz | 2.49 mm | 10.0 mm | 5.0 mm | 20 mm | 80 mm | 32.7 ps / 33.4 ps |
| 94 GHz | 0.80 mm | 3.19 mm | 1.6 mm | 6.4 mm | 25.5 mm | 10.4 ps / 10.6 ps |
| 1550 nm optical | 387 nm | 1.55 µm | 775 nm | 3.1 µm | 12.4 µm | 5.07 fs / 5.17 fs |

This mapping holds when the PhC delays the same carrier that is radiated: an optical phased array at 1550 nm, or a scaled microwave/mm-wave PhC delay structure feeding antennas at its own operating frequency. **It does not describe RF-photonics:** if the optical device delays an RF-modulated optical carrier that feeds an RF subarray, the required delay is RF cycles (33.4 ps at 30 GHz, about 6600 optical periods), the 5 fs period quantum is negligible, and the optical carrier-phase intercept is irrelevant to the RF phase. The "one carrier cycle per period" argument is specific to the carrier-delay mapping. The 2-D dielectric structure is not claimed to be a microwave antenna feed network.

## 9. Restrictions (Phase 4 and 5 restrictions travel unchanged)

Domain-limited to mask_common 0.231–0.268; 2-D, lossless, uniform, no disorder, coupling, bends, thermal or length-dependent loss; four one-period channels of one flat-n_g D0 mode; raw PhC has constant intercept errors up to ≈ 2.7 rad and needs one constant trim per channel (trimmed ≠ ideal TTD); common minimum run time 3500 and estimator choices from Phase 4; subarray replica lobes remain (≤ −16.6 dB); the fixed intra-subarray phase gives a residual ±0.16° squint even for ideal TTD; Monte Carlo uses uniform independent errors on the measured transfer functions, with no fabrication model, so budgets are system-level estimates; Msub = 8 with 2-period steps and fractional periods were not simulated; uncalibrated absolute-length errors are severe.

## 10. Unresolved issues

* No fabrication-disorder FDTD or coupling between adjacent channels; a length error shifts both the intercept and the slope, and the realistic calibrated/uncalibrated mix is unknown.
* Msub = 4 is selected for this delay quantum and 30° steering; other steering angles change the required step (Msub·λ0/4·sinθ0/0.5) and the alignment.
* Channel-to-channel bandwidth (common band) is only 14.8 %; the usable band of the tolerance results is that same band.
* Phase 7 documentation/finalization is not done.
