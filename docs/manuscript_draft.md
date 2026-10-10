# Subarray-Level Delay Lines from W1 Photonic-Crystal Waveguides for Beam-Squint Reduction: A Computational Study

*Draft for review. Computational study only: no experiment, no hardware, no fabrication model. Numbers in this draft are generated from `results/final/final_results.json` (see `docs/final_tables.md`); `[REF]` marks places where a literature citation is needed and none has been supplied.*

Candidate titles considered: (1) "Photonic-Crystal True-Time Delay for Broadband Subarray Beamforming" (rejected: "broadband" is 14.8 %, "true-time delay" needs a trim qualifier); (2) "W1 Photonic-Crystal Delay Lines for Broadband Phased-Array Subarray Steering" (rejected: same); (3) "Broadband Beam-Squint Suppression with Photonic-Crystal Subarray Delay Lines" (rejected: overstates both); (4) the title above, chosen as the most conservative that the evidence supports.

## Abstract

Beam squint in wideband phased arrays is removed by true-time delay (TTD), and photonic-crystal (PhC) line-defect waveguides are a candidate delay medium because of their slow-light group index. We study this in a two-dimensional, lossless model combining plane-wave-expansion band structures and finite-difference time-domain (FDTD) simulation of a plain W1 waveguide (triangular air-hole lattice, ε = 12, r/a = 0.30, TE). After verifying the FDTD delay extraction against analytic and band-structure references, we find a group index of about 3.9, flat to ±10 % over a validated band of 0.231–0.268 c/a (14.8 % fractional bandwidth, a domain-limited lower bound). Four directly simulated sections of 8–11 periods give delay increments of 3.931 a/c per period (within 0.04 % of the target by band fit, 0.981 carrier cycle at band center), amplitude imbalance below 0.023 dB, and a large constant phase-intercept error (up to 2.7 rad) that one constant trim per channel reduces to 0.026 rad. Used as element-level delays of a four-element, 30°-steered array, the delay quantum of about one carrier cycle forces an element pitch of 1.96 λ0 and equal-strength (0 dB) grating lobes, although the trimmed channels track ideal TTD to 0.05°. The same quantum matches the delay required by four-element λ0/2 subarrays at 30° (1.000 cycle) to 1.9 %. With one PhC channel per subarray and one constant phase trim per channel, the maximum pointing error over the band is 0.19° against 2.69° for fixed phase (R_squint = 0.070), the strongest ordinary sidelobe is −12.4 dB, and the residual subarray replica lobes are at most −16.6 dB. System-level perturbation estimates give a P95 pointing error of 0.35° for a combined phase, delay and amplitude perturbation. The result supports the feasibility of subarray-level PhC delay lines for beam-squint reduction. It is limited by the 2-D lossless model, the validated band, the required phase trim, residual replica lobes, a single steering angle, and the absence of loss, disorder, coupling and experiment.

## 1. Introduction

Fixed phase shifters steer a beam correctly at one frequency only; the beam angle drifts with frequency (beam squint) [REF]. True-time delay removes squint but requires delay lines whose length grows with the aperture and the delay [REF]. Slow-light PhC waveguides reduce the physical length per unit delay [REF], and PhC delay lines have been proposed for beamforming [REF]. Whether a PhC delay line is useful for broadband beamforming depends on the delay–bandwidth–dispersion trade, on how delay is quantized across channels, and on what the beam pattern does with the actual (not ideal) channel transfer functions.

This study asks, in a controlled 2-D computational setting, what delay, bandwidth and dispersion a W1 PhC waveguide provides, whether a family of four channels built from it behaves as a delay family, and what happens to the beam when those simulated channels are used as array weights. The answer was not assumed to be favorable. The study found one unfavorable result (element-level use fails because of grating lobes) and one favorable architectural result (a four-element subarray matches the delay quantum).

**Contributions of this study** (not claimed as novel relative to the literature, which has not been surveyed here): (i) a verified FDTD/MPB workflow for extracting complex transfer functions and delays of finite W1 sections, with documented numerical-quality events; (ii) a quantified delay–bandwidth–dispersion comparison of a plain W1 with two first-row-engineered variants; (iii) a four-channel family built from directly simulated complex transfer functions; (iv) identification of the delay-quantum limitation at the element level; (v) the observation that the one-period quantum nearly matches the delay required by a four-element λ0/2 subarray at 30°; (vi) system-level sensitivity estimates for that architecture. **Standard theory used:** the array factor, beam squint from fixed phase, the Floquet–Bloch band theory of PhC waveguides, and the equivalence of a frequency-linear phase to a delay.

## 2. Model and Methods

**Platform.** Triangular lattice of air holes in ε = 12, r/a = 0.30, TE polarization (Hz out of plane), W1 line defect (one row removed), 2-D, lossless, non-dispersive. Normalized units a = c = 1, f = a/λ. The geometry family was selected from candidates by MPB screening with a pre-declared rule (a moderate-margin decision; r/a = 0.30 is an engineering choice made after the rule's output of 0.40, both recorded). Figure 1 shows the structure and architecture.

**Solvers.** MPB for band structures and group index; MEEP FDTD (version 1.33.0-beta) for finite devices. The guided mode is projected from the recorded fields onto the MPB forward and backward modes by a Lorentz-reciprocity overlap, giving S21 = a_f(out)/a_f(in) between reference planes three periods outside the device. Delay is τ_g = d arg H / dω; the sign convention (MEEP exp(−iωt)) was established numerically. Production resolution is 30 pixels per a with a 40 pixel confirmation; DFT decay tolerance 1e-7.

**Array model.** Beam patterns are computed analytically from the array factor with the simulated complex H_m(f) as weights (isotropic elements, no mutual coupling), over 18001–20001 angle samples (Phases 5–6). Beam metrics: design-lobe angle, pointing error, R_squint (maximum pointing error divided by that of fixed phase), half-power beamwidth, peak power, strongest other lobe.

**Terms.** *Phase-intercept trim* is a frequency-independent phase correction chosen at f0 = 0.2495; it does not change amplitude or group delay. *True-time delay* is used for a phase slope corresponding to a physical delay; trimmed channels are described as TTD-slope channels with a constant trim and are never called ideal TTD.

## 3. Numerical Verification

The FDTD extraction was verified in a ladder (Figure 2). In vacuum and homogeneous dielectrics the extracted delay equals L/c and nL/c up to the closed-form Yee-grid correction, with error converging at second order (vacuum mean error 2.2e-3 a/c at 40 px/a for L = 10). An eight-period quarter-wave Bragg stack agrees with the transfer-matrix method, with the median relative error in the complex transmission falling from 2.6e-2 to 1.4e-3 between 80 and 320 px/a (order about 1.9). A finite 2-D PhC stop band agrees with the parity-resolved MPB gap to within 3.5e-3 for every variant studied, dominated by finite-length and termination effects.

For the W1, the FDTD group index agrees with the MPB value in 0.231–0.268 with a median bias of +0.3 % and a maximum of 1.0 % (res 40, tolerance 1e-7). The default decay tolerance (1e-6) gave 5–14 % rms noise and is not used. The upper edge 0.268 of the validated band was chosen after inspecting the data (Section 10). A time-domain pulse check agrees with the phase-derived delay to 0.4–0.9 %.

## 4. Photonic-Crystal Delay Line

The guided Hz-even band (MPB band 8) has n_g,min = 3.87 at f = 0.2564. In the validated band the plain W1 (D0) has a median delay of 3.932 a/c per period, flat to ±10 % over 0.231–0.268 (14.8 % fractional bandwidth), with transmission above 0.99 and reflection below 0.3 % in band. Both band ends coincide with validated-range limits, so 14.8 % is a lower bound for flatness within the model; the same rule without the restriction gives 18.0 % but that range is not validated. Two first-row-engineered variants (D1: first-row radius ×1.2, shift −0.05; D2: ×1.2, −0.10) were simulated at N = 4, 8 and 12 (Figure 3): they raise the median delay by 7 % and 14 % but reduce the flat-delay bandwidth to 8.4 % and 5.6 % (−43 %, −62 %) and increase the group-delay dispersion by factors of 4 and 8 (band-average −3.86 and −7.07 against −0.91 (a/c)² per period). Higher group indices exist only within a few percent of the engineered band edge, where delay is far from flat. D0, the plain W1, was therefore used downstream. **The delay line is not dispersion-engineered in this study; first-row engineering was tested and rejected.** Three designs are an illustration, not a Pareto frontier.

Relative to a square-lattice rod waveguide reference (MPB only, no FDTD), the plain W1 needs about 2.1 times fewer free-space wavelengths of length for the same delay at comparable fractional bandwidth. No loss, coupling or FDTD baseline comparison was performed.

## 5. Multi-Channel Delay Family

Four sections of the same plain W1 with N = 8, 9, 10, 11 periods were simulated directly (Figure 4); channels were not built by cascading. Over the validated band, the band-fit adjacent delay increments are 3.9315, 3.9308 and 3.9310 a/c against the 3.9323 target (errors −0.02 % to −0.04 %). The pointwise increments are not that flat: a cubic fit of the phase gives increments between about 3.88 and 4.18 a/c across the band, following the rise of n_g toward the lower edge. Transmission is 0.997–1.0025 (the excess over unity is numerical) and the amplitude imbalance among channels is at most 0.023 dB. A cascade prediction from shorter devices reproduces longer ones to 0.006 rad and 0.5 %.

The raw channels are not the target phase: the per-period phase advance differs from 2πf0·τ_target by −1.2035 rad, giving frequency-independent intercept errors of 0, −1.20, −2.41 and +2.67 rad (m = 0–3). One constant trim per channel, chosen at f0, reduces the phase error to at most 0.026 rad (1.5°) over the band; what remains is mostly the curvature of n_g(f). Trimming leaves amplitude and group delay unchanged, so the trimmed family is a delay family with a fixed phase added; it is not ideal TTD. Because all channels are one uniform medium at different lengths, these favorable numbers reflect one flat-n_g mode, not four independently engineered elements.

## 6. Beamforming Results: Element-Level Array

Four isotropic elements steered to 30° with the Phase 4 transfer functions as weights (Figure 5). The delay quantum is 0.981 cycle per period while a half-wave array at 30° needs 0.25 cycle, so the element pitch is set by the quantum: d = τ/sin30° = 7.86 a = 1.96 λ0. Fixed phase squints by +2.69° to −2.26° over the band (4.9° peak to peak, about 0.6 beamwidths). Ideal TTD holds 30°. The raw PhC channels suppress the squint (offset-removed spread 0.2 of fixed phase) but point 6.3° low because of the intercept error. The trimmed channels stay within 0.050° of ideal TTD (R_squint 0.019). But at 1.96 λ0 pitch, every case has grating lobes as strong as the design lobe (0 dB), with three to five equal lobes across the band: the result suppresses squint of the design lobe and does not make a usable aperture. This is a limitation imposed by the delay quantum; the delay lines themselves track ideal TTD.

## 7. Subarray Architecture

Each PhC channel feeds a subarray of Msub λ0/2-spaced elements with fixed intra-subarray steering phases; the PhC provides the inter-subarray delay, so the required step is Msub·(λ0/2)·sin30°/c = Msub/4 cycles (Figure 6). The one-period quantum (0.981 cycle) matches this to 1.9 % only for Msub = 4 (1.000 cycle required, mismatch −1.89 %); Msub = 1, 2 and 8 mismatch by +292 %, +96 % and −51 %. This match is a coincidence of n_g·f0 ≈ 0.98 for this geometry, carrier and steering angle, not a design principle; at other angles the required step changes, and a conventional line can be cut to any length. Msub = 8 with two-period steps was not simulated.

Results for Msub = 4 (16 elements, 4 PhC channels, subarray pitch 2λ0), with a constant trim re-chosen for the required 4.008 a step (Figures 7–8):

* maximum pointing error 0.189° (rms 0.110°) over 0.231–0.268 against 2.69° for fixed phase; R_squint = 0.070; tracking tolerance 0.367° (5 % of the half-power beamwidth) met. Even ideal subarray TTD has 0.156° because the fixed intra-subarray phases squint; the PhC adds 0.055°.
* the equal-strength 0 dB grating lobes are replaced by one design lobe and ordinary sidelobes (strongest other lobe −12.4 dB; ideal subarray TTD −12.4 dB). Subarray replica lobes at the 2λ0 pitch remain, at most −16.6 dB (growing from a null at f0), suppressed by the squinting λ0/2 element factor. This is mitigation, not elimination.
* minimum peak power 0.983 (−0.075 dB); half-power beamwidth 7.95°, 7.35°, 6.83° at the band low, center and high; resolutions 30 and 40 agree to 0.007° and 0.04 dB.
* the raw PhC (without trim) points 6.7° off (strongest lobe −6.7 dB); the Phase 4 trim unchanged (set for 3.9323 a rather than 4.008 a) gives 0.73° error. A constant phase trim re-chosen for the architecture is required.
* element-level use with a constant trim (Msub = 1) fails: a constant phase cannot supply the missing 0.75 cycle of group delay, and the squint is worse than fixed phase.

## 8. Robustness and Tolerances

Part of the analysis perturbs the measured transfer functions mathematically (uniform independent errors per channel; Monte Carlo seed 20261010, 1000 samples per point; Figure 9). These are **system-level tolerance estimates, not fabrication tolerances; no fabrication-disorder FDTD was run.** The P95 of the maximum pointing error is 0.35° for a combined perturbation of ±2° phase, ±0.5 % delay, ±0.25 dB amplitude (nominal 0.19°). A P95 error of 0.5° is reached at ±7.6° of trim phase error, ±1.2 % uncalibrated delay error (relative to q·τ_PhC), and not within ±6 dB of amplitude error; if the intercept is re-trimmed after a delay change, ±10 % delay error stays below 0.5°. In absolute-length terms the uncalibrated limit is much tighter (about 0.14 % of the total device delay, or about 0.02 carrier cycle): a length error is benign only if the constant phase is re-trimmed afterward. Amplitude mismatch is not a pointing limiter in this range.

## 9. Discussion

*Physical interpretation.* The delay of a uniform W1 section is τ = n_g·N·a/c, so the delay per period in carrier cycles is n_g·f0 (0.981 here). The quantum is coarse because one period already delays about one cycle. Slow light helps footprint (here n_g ≈ 3.9 against about 2 for the rod-waveguide reference) but also makes the minimum delay step large. At element level that forces a sparse array; at subarray level (one channel per four λ0/2 elements) the same step is what is needed at 30°. The constant phase trim is required because the Bloch phase of the W1 mode gives each channel a frequency-independent intercept error; the trim supplies it and the PhC supplies the slope.

*What the PhC contributes.* The beam-squint reduction relative to fixed phase is the generic benefit of any TTD. What the study shows is that the simulated PhC channels reproduce ideal subarray TTD pointing to 0.055° over a 14.8 % validated band. It does not show a performance advantage over a conventional delay line; the only comparison with a conventional guide is MPB-only.

*Carrier mapping.* The normalized results can be scaled (for example a = 2.49 mm at 30 GHz, 0.80 mm at 94 GHz, 387 nm at 1550 nm) when the PhC delays the radiated carrier. They do not describe RF-photonic beamforming, where the optical device delays an RF-modulated carrier and the required delay is RF cycles; the one-cycle-per-period argument is specific to the carrier-delay mapping. The 2-D structure is not claimed to be a microwave antenna feed network.

## 10. Limitations

* 2-D simulation (TE, no out-of-plane radiation); lossless, non-dispersive dielectric; no propagation loss; no fabrication disorder; no channel coupling; no bends; no thermal effects; no realistic input/output couplers; isotropic elements without mutual coupling. These bound every conclusion.
* Finite validated band only: 0.231–0.268 (14.8 %). The upper edge 0.268 was chosen after inspecting the Phase 2 data (the projection residual, unitarity defect and n_g bias grow above about 0.265–0.268; the cause was only partly identified), so the band is domain-limited, not a measured physical band edge.
* Post-data corrections in the channel-family phase: the first FDTD attempt was discarded for under-convergence (erratic stop times; one channel had a 1.9 % unitarity defect), a common minimum run time was added after seeing the data, and the delay-increment estimator was changed after seeing that the pre-declared pointwise central difference was noise-limited. Delay claims therefore rest on band fits and phase slopes. A cascade step-count error in the original criteria text was found and corrected.
* A constant phase trim per channel (raw intercept errors up to 2.7 rad) is required, and how it is realized physically is not addressed.
* Residual subarray replica lobes (strongest −16.6 dB); the count of "visible" replicas depends on whether grazing lobes at ±90° are counted. The fixed intra-subarray phase leaves a residual ±0.16° squint even for ideal TTD.
* Only four PhC channels of one flat-n_g mode at four lengths; only one steering angle (30°) and one match (Msub = 4, one-period step); Msub = 8 with two-period steps was not simulated.
* Robustness numbers are system-level estimates on the measured transfer functions, not fabrication tolerances.
* Two-resolution (30 and 40) convergence, not a full series. Pointwise group-delay derivatives and local GDD are noise-limited.
* No conventional-delay-line comparison beyond an MPB-only reference; no literature survey or novelty check.
* No experiment.

**Numerical-quality events retained for transparency.** The Phase 2 validated-band boundary was chosen post hoc. A sign bug in the hole-shift code (negative shifts built as outward shifts) was found when a design transmitted below its band edge; the runs were discarded, the code fixed with a regression test, and the runs repeated. A similar sign bug had been found and fixed in the MPB geometry earlier. The first channel-family FDTD attempt was discarded. A first Monte Carlo run reported peak-power loss with the wrong sign and was discarded and rerun. Several diagnostics (replica-lobe level, efficiency loss, calibrated absolute-length delay) were added after seeing results and are labeled post hoc. These are examples of numerical quality control, not hidden steps.

## 11. Conclusion

In a 2-D lossless model, a plain W1 photonic-crystal waveguide provides delay increments of about one carrier cycle per period, accurate to 0.04 % in band fit, over a validated 14.8 % fractional bandwidth. Used at element level, this quantum forces a 1.96 λ0 pitch and 0 dB grating lobes. Used for four-element λ0/2 subarrays at 30°, it matches the required delay to 1.9 %; with one constant phase trim per channel the maximum pointing error over the band is 0.19° against 2.69° for fixed phase, with sidelobes at −12.4 dB and replica lobes at most −16.6 dB. This supports the feasibility of subarray-level PhC delay lines for beam-squint reduction. Establishing an advantage over conventional delay lines, realistic loss, coupling and disorder, and the realization of the phase trim are the next steps.

## References

[REF] — to be supplied: beam squint and true-time-delay beamforming; PhC slow-light waveguides and delay lines; PhC delay lines for beamforming; MPB; MEEP; subarray TTD architectures and quantization lobes.

## Figures

1. `figures/final/fig01_concept` — W1 channel and subarray architecture (2-D FDTD domain vs analytic beamforming model).
2. `fig02_verification` — FDTD verification (convergence, Bragg vs TMM, W1 n_g vs MPB).
3. `fig03_tradeoff` — D0/D1/D2 delay–bandwidth–dispersion.
4. `fig04_channels` — four-channel transfer functions.
5. `fig05_phase5_limit` — element-level array and grating lobes.
6. `fig06_architecture` — delay match, quantization, sidelobe level versus Msub.
7. `fig07_patterns` — Msub = 4 patterns at band low, center, high.
8. `fig08_pointing` — pointing versus frequency.
9. `fig09_robustness` — system-level perturbation sensitivity.
