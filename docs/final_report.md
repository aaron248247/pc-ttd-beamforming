# Final project report — photonic-crystal delay lines for subarray-level beam-squint reduction

Final repository state for Phases 0–7 (computational study; no experiment). Numbers come from `results/final/final_results.json` (`docs/final_tables.md`); figures are in `figures/final/`. The claim audit and the list of contradictions found are in `docs/phase7_audit.md`; the manuscript draft is `docs/manuscript_draft.md`.

## 1. Abstract

A two-dimensional, lossless FDTD/MPB study of a plain W1 photonic-crystal (PhC) waveguide (triangular air holes, ε = 12, r/a = 0.30, TE) as a delay medium for phased-array beamforming. After a verification ladder, the plain W1 gives a group index of about 3.9 flat to ±10 % over a validated 0.231–0.268 c/a band (14.8 % fractional bandwidth, domain-limited). Four directly simulated channels (8–11 periods) differ by 3.931 a/c per period (within 0.04 % of target by band fit; 0.981 carrier cycle at band center), with a large constant phase-intercept error (up to 2.7 rad) that one constant trim per channel reduces to 0.026 rad. As element-level delays of a 30°-steered four-element array the trimmed channels track ideal TTD to 0.05°, but the delay quantum forces a 1.96 λ0 pitch and 0 dB grating lobes. The same quantum matches the delay of four-element λ0/2 subarrays at 30° to 1.9 %; with one channel and one trim per subarray, the maximum pointing error over the band is 0.19° against 2.69° for fixed phase (R_squint 0.070), the strongest sidelobe is −12.4 dB and replica lobes are at most −16.6 dB. This is a worked near-match example (exact match at 29.4°; useful window about 26–34° under the stricter post hoc criteria, from the later angle sweep), not a broadly scannable architecture, and it supports feasibility only within the model's limits (Section 15).

## 2. Motivation

Fixed phase shifters squint with frequency; true-time delay does not. PhC slow-light waveguides promise short delay lines, but bandwidth, dispersion, delay quantization and the actual (non-ideal) channel responses decide whether they help a beam. The project asked, without assuming a favorable answer, when and why PhC delays improve wideband steering (the validated band is 14.8 %).

## 3. Research question

What delay, bandwidth and dispersion does a PhC line-defect waveguide give, does a family of such channels behave as a delay family, and what happens to the beam when the simulated channels are used as array weights?

## 4. Photonic-crystal platform

Triangular lattice of air holes in ε = 12, TE, W1 defect (Figure 1a). The family was chosen in Phase 0.5 by MPB screening; the decision was moderate-margin (the square-lattice rod waveguide R has a much wider single-mode band; the hole lattice wins on flat slow-light potential and FDTD cost). The baseline radius r/a = 0.30 is a post-selection engineering choice (the rule's output was 0.40; both are recorded). The rod waveguide R was kept as an MPB-only conventional reference.

## 5. Numerical methods

MPB for bands and n_g; MEEP FDTD for finite devices; the guided mode is projected from DFT plane fields onto MPB modes with a Lorentz-reciprocity overlap (MEEP's own mode-coefficient routine did not reproduce the supercell bands and was abandoned). Delay is the derivative of the unwrapped phase; sign convention established numerically. Resolution 30 (production) and 40 (confirmation); DFT decay tolerance 1e-7 (1e-6 gave 5–14 % rms noise). Beam patterns are analytic array factors with the simulated H_m(f) as weights.

## 6. Verification (Figure 2)

Vacuum and dielectric delays equal L/c and nL/c up to the closed-form Yee correction, with second-order convergence. An eight-period Bragg stack agrees with the transfer-matrix method (complex-transmission error falls with order about 1.9). A 2-D PhC stop band agrees with the parity-resolved MPB gap within 3.5e-3. Two derivative-estimator pitfalls and several run-termination pitfalls were found and fixed; they are listed in the Phase 1 report.

## 7. W1 delay characterization

FDTD n_g agrees with MPB in 0.231–0.268 with a median bias of +0.3 % and a maximum of about 1 % (res 40). Delay is linear in length (phase residual ≤ 0.0045 rad over three lengths). A time-domain pulse agrees with the phase-derived delay to 0.4–0.9 %. The upper edge 0.268 was chosen after the data were seen; above it, the projection residual, the unitarity defect and the n_g bias grow together, with the cause only partly identified.

## 8. Delay-bandwidth trade study (Figure 3)

| design | flat-delay band | fractional BW | median delay per period | GDD per period |
|---|---|---|---|---|
| D0 plain W1 | 0.231–0.268 | 14.8 % | 3.93 a/c | −0.91 (a/c)² |
| D1 (first row ×1.2, shift −0.05) | 0.2465–0.268 | 8.4 % | 4.20 | −3.86 |
| D2 (×1.2, −0.10) | 0.2535–0.268 | 5.6 % | 4.47 | −7.07 |

Engineering the first rows bought 7–14 % more delay at the price of 43–62 % of the bandwidth and 4–8 times the dispersion. D0 was selected. **The delay line used is therefore a plain W1, not a dispersion-engineered one**; engineering was tried and rejected. Three designs are an illustration, not a frontier, and every band is bounded by the validated range on the high side.

## 9. Four-channel family (Figure 4)

N = 8, 9, 10, 11, each a separate FDTD run. Band-fit increments 3.9315/3.9308/3.9310 a/c (target 3.9323), amplitude imbalance ≤ 0.023 dB, transmission 0.997–1.0025. Raw channels carry a frequency-independent intercept error of 0, −1.20, −2.41, +2.67 rad (a consequence of the Bloch phase, not an analysis error); one constant trim per channel leaves ≤ 0.026 rad, mostly the curvature of n_g(f). Pointwise increments deviate by several percent inside the band; only the band-fit values are within 0.04 %. All channels are one flat-n_g mode at four lengths, so the favorable numbers are not four independent engineering successes. The first FDTD attempt was discarded (under-converged); a common minimum run time and a different increment estimator were adopted after seeing data and are disclosed.

## 10. Element-level array (Figure 5)

Four isotropic elements steered to 30°, actual H_m as weights. Fixed phase squints +2.69° to −2.26°; ideal TTD holds 30°; raw PhC suppresses squint but points 6.3° low; trimmed PhC stays within 0.050° of ideal TTD (R_squint 0.019). The pitch is forced to 1.96 λ0 by the one-cycle quantum, and every case has 0 dB grating lobes. The conclusion was "scientifically promising for squint suppression, but not a usable aperture".

## 11. Delay-granularity limitation

A half-wave array at 30° needs 0.25 cycle of delay per element; the PhC supplies 0.981 per period. A constant phase cannot supply the missing group delay (a trimmed Msub = 1 hybrid squints more than fixed phase: 7.7° error). Sub-period delays and fractional-period lengths were not simulated.

## 12. Subarray-level architecture (Figures 6–8)

One PhC channel per Msub λ0/2-spaced elements. Required step Msub/4 cycles: 0.25, 0.50, 1.00, 2.00 for Msub = 1, 2, 4, 8; PhC step 0.981. Only Msub = 4 matches (−1.89 %); the others mismatch by +292 %, +96 %, −51 %. For Msub = 4 with a re-trimmed constant phase per channel: max pointing error 0.189° (fixed phase 2.69°; ideal subarray TTD 0.156°), R_squint 0.070, strongest other lobe −12.4 dB, strongest subarray replica −16.6 dB, minimum peak power −0.075 dB, HPBW 7.35° at f0. The match is a coincidence of n_g·f0 ≈ 0.98 at this carrier and angle; replica lobes remain; res 30 and 40 agree to 0.007°. Used unchanged, the Phase 4 trim gives 0.73° error (above the 0.37° tolerance), so the trim must be re-chosen for the architecture; with no trim the beam points 6.7° off.

## 13. Robustness (Figure 9)

System-level only: mathematical perturbations on the measured channel transfer functions (Monte Carlo seed 20261010, 1000 samples per point). Combined test (±2°, ±0.5 %, ±0.25 dB): P95 maximum error 0.35°. P95 = 0.5° at ±7.6° trim phase error, ±1.2 % uncalibrated delay error, and beyond ±6 dB amplitude error. Uncalibrated absolute length errors are severe (about 0.14 % of the device delay). These are system-level tolerance estimates, not fabrication tolerances.

## 14. Physical interpretation

Delay per period in carrier cycles is n_g·f0. A slower mode increases both footprint savings and the delay quantum. At element level the quantum implies a sparse array and grating lobes; at subarray level the same quantum is needed at 30°. The scaling mapping (a = 2.49 mm at 30 GHz, 0.80 mm at 94 GHz, 387 nm at 1550 nm) applies when the PhC delays the radiated carrier; it does not describe RF-photonic beamforming, and the 2-D structure is not claimed to be a microwave feed network. Footprint: one carrier cycle needs 0.254 λ0 in the W1, 1.0 λ0 in air and 0.289 λ0 in an ideal homogeneous ε = 12 medium (the W1 is 11.9 % shorter; the ε = 12 row is an analytic benchmark, not a modeled waveguide), so a large footprint advantage over a high-index propagation reference is not demonstrated. The earlier MPB-only rod-waveguide comparison (about 2.1×) has been withdrawn from the manuscript as unlike-for-like.

## 15. Limitations

2-D; lossless; no loss, disorder, channel coupling, bends, thermal effects or realistic couplers; validated band only (upper edge chosen post hoc); post-data Phase 4 convergence correction and estimator change; constant phase trim required; residual replica lobes; system-level robustness only; four channels of one mode; the angle sweep reuses the same channels and shows a narrow useful window; Msub = 8 with two-period steps not simulated; no simulated conventional line; no experiment; no novelty claim (prior art: PhC TTD lines, subarray TTD and the required-delay relation are published; see `docs/closest_prior_art.md`). Full list: `docs/final_tables.md`, Table 7.

## 16. Conclusions

Within a 2-D lossless model, a plain W1 PhC delay line gives about one carrier cycle of delay per period with 0.04 % band-fit accuracy over a validated 14.8 % band; it fails as an element-level delay (grating lobes) and works as a four-element-subarray delay (1.9 % match, 0.19° pointing error, −12.4 dB sidelobes, ≤ −16.6 dB replicas) with one constant phase trim per channel. The result supports feasibility near a matched angle; it does not establish an advantage over conventional delay lines or a realizable device.

## 17. Reproducibility

See `docs/reproducibility.md`. Final figures, tables and handoff regenerate from tracked results in under two minutes (`experiments/phase7/7a…7d`); the test suite runs in about two and a half minutes.

---

### Historical corrections kept (numerical quality control)

Phase 0.5 geometry sign bug and branch-tracking bug (rerun); Phase 2 post hoc validated-band boundary and abandoned built-in mode solver; Phase 3 sign bug in the FDTD hole-shift code (D1/D2 runs discarded and rerun, regression test added); Phase 4 first FDTD attempt discarded for under-convergence, run-time and estimator changes after seeing data, cascade step-count erratum; Phase 6 first Monte Carlo sign convention discarded and rerun, post hoc diagnostic additions. Details in `docs/decisions.md` and the phase reports.
