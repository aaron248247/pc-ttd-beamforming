# Delay Quantum of W1 Photonic-Crystal Waveguides in Subarray Beamforming: A Two-Dimensional Computational Study

*Computational study only: no experiment, hardware, or fabrication model. Numerical values are generated from `results/final/final_results.json`, `results/final/angle_sweep.json`, and `results/final/conventional_reference.json` and summarized in `docs/final_tables.md`. Citations marked † still require bibliographic confirmation; `[CITATION NEEDED]` marks an unresolved source.*

## Abstract

Beam squint in wideband phased arrays is reduced by true-time delay, and photonic-crystal (PhC) slow-light waveguides have been used as delay lines. The simulated complex transfer functions of a specific PhC delay family are used directly as array weights. The model is two-dimensional and lossless: a plain W1 line-defect waveguide (triangular air-hole lattice, ε = 12, TE) treated with plane-wave band-structure and FDTD simulation. One lattice period delays the carrier by about one cycle (0.981 cycle at band center). This step is too coarse for element-level steering at λ0/2 spacing: using the four simulated channels element by element forces an effective pitch of 1.96 λ0 and produces grating lobes at about 0 dB. The same delay step is close to what four-element λ0/2 subarrays require near 30° (1.000 cycle required, −1.9 % mismatch; exact match at 29.4°). In that worked example, with one constant phase trim per channel, the maximum pointing error over the validated 14.8 % band (0.231–0.268 c/a) is 0.19° against 2.69° for fixed phase. The useful angular window is narrow: about 26–34° under the stricter post hoc criteria. Compared with ideal propagation in a homogeneous ε = 12 medium, the selected W1 shortens a one-cycle delay by only 12 %. The study is a computational case study of delay granularity and hierarchy selection, not a demonstration of a generally superior or broadly scannable PhC beamformer. It is limited by the 2-D lossless model, the narrow validated band, the required phase trim, and the absence of loss, disorder, coupling and experiment.

## 1. Introduction

Beam squint arises when fixed phase shifters steer a wideband array: the beam direction drifts with frequency. Replacing phase shifts by true-time delay (TTD) removes the frequency dependence of the steering phase in principle, and TTD for wideband phased arrays, including photonic and microwave-photonic implementations, is an established subject [Capmany & Novak 2007; Yao 2009; Mailloux 2005‡]. Integrated optical beamforming networks have used ring-resonator delay lines, for example a chip whose four delay lines are configured to feed a 16-element array with continuously tunable delay [Burla et al. 2013, abstract-level reading only].

PhC waveguides are an established route to compact delay. Line-defect (W1) slow-light waveguides and their dispersion were characterized experimentally two decades ago [Notomi et al. 2001; Vlasov et al. 2005], their bandwidth and dispersion engineering is well studied [Frandsen et al. 2006; Li et al. 2008; Mori et al. 2007; Povinelli et al. 2005; Krauss 2007; Baba 2008†], and the general delay–bandwidth limits of slow light are known [Tucker et al. 2005]. PhC waveguide delay lines have also been demonstrated, not merely proposed, as true-time-delay elements for phased arrays: a silicon-nanomembrane PhC waveguide array with four wavelength-tunable channels [Lin et al. 2012] and one-dimensional fishbone PhC delay lines [Chung et al. 2018]. Those works delay the RF envelope of an optical carrier and tune the delay continuously with wavelength.

Subarray-level TTD combined with element-level phase shifters is likewise a standard architecture. Its quantization (replica) lobes are textbook material [Mailloux 2005‡], and recent hybrid delay–phase work states the required inter-subarray delay and sizes the subarray from the delay range available. In particular Nguyen and Kim give the delay of the m-th TTD as ((2m−1)N − 1)ψ/(4 f_c), so adjacent TTDs differ by N sinθ/2 carrier cycles at half-wave spacing, and a criterion linking subarray size, maximum delay, angle and carrier [Nguyen & Kim, arXiv:2212.07484, Eqs. 24–29 of v7, to be checked against the typeset version; see also Dai et al. 2022].

None of this is new here. The question addressed is narrower and computational: given one specific, physically simulated W1 transfer-function family, what happens when those actual complex channel responses are used directly in a phased-array model, and how does the finite delay step constrain element-level versus subarray-level use? In the prior PhC work the delay is tuned continuously, so the minimum step of a fixed-length section does not arise. We did not identify, within the literature searched, a study combining verified W1 transfer functions with this particular full-wave-to-array workflow. The search was a web search without access to IEEE Xplore, Scopus or Web of Science, was strongest for titles and abstracts, and could miss conference papers, theses and antenna-literature work; absence from it is not evidence of absence. The full texts of Lin et al., Chung et al. and Nguyen and Kim were read; Burla et al. was read at abstract level and Mailloux was not read.

**Scope of the study.** This work evaluates: (i) a verified 2-D W1 delay model (FDTD checked against analytic and band-structure references); (ii) extraction of the complex channel transfer functions of four sections; (iii) direct propagation of those responses into array-level beam models; (iv) the element-level result, which is negative because the delay step is coarse; (v) a worked subarray-level case near a matched angle; (vi) the dependence of that case on steering angle; (vii) system-level robustness estimates; and (viii) a length comparison with ideal propagation references. The theory used is standard: the array factor, beam squint from fixed phase, Floquet–Bloch band theory, and the equivalence of a frequency-linear phase to a delay.

## 2. Model and Methods

**Platform.** Triangular lattice of air holes in ε = 12, r/a = 0.30, TE polarization (H_z out of plane), W1 line defect (one row removed); 2-D, lossless, non-dispersive. Normalized units a = c = 1, f = a/λ. The lattice family was selected from candidates by MPB screening with a pre-declared rule; the baseline radius r/a = 0.30 is an engineering choice made after the rule returned 0.40 (both recorded in the repository). Figure 1 shows the structure and the two array architectures.

**Solvers.** Band structures and group index: MPB [Johnson & Joannopoulos 2001]. Finite devices: MEEP FDTD, version 1.33.0-beta [Oskooi et al. 2010]. The guided mode is projected from recorded fields onto the MPB forward and backward modes by a Lorentz-reciprocity overlap, giving S21 = a_f(out)/a_f(in) between reference planes three periods outside the device. Group delay is τ_g = d arg H/dω. We use the MEEP convention exp(−iωt), so a delay τ corresponds to exp(+iωτ); this sign was checked numerically. Production resolution is 30 pixels per a, with 40 pixels as confirmation; DFT decay tolerance 1e-7.

**Array model.** Beam patterns are computed analytically from the array factor with the simulated complex H_m(f) as weights (isotropic elements, no mutual coupling). Metrics: design-lobe angle, pointing error, R_squint (maximum pointing error divided by that of fixed phase), half-power beamwidth (HPBW), peak power, strongest other lobe.

**Slope and intercept.** A transfer function has a group-delay slope and a constant phase intercept. A constant phase trim changes the intercept and does not change the group delay. The trimmed architecture is therefore a "TTD slope plus constant phase trim", chosen at f0 = 0.2495; it is never called ideal TTD.

**Standard array relations used.** For half-wave elements, adjacent subarrays of Msub elements steered to θ0 need a delay of q_req = Msub sinθ0/2 carrier cycles. A W1 section of k periods supplies q_phc = k f0 τ_PhC, where τ_PhC is the delay per period. Requiring q_req = q_phc gives the match condition Msub sinθ_match/2 = k f0 τ_PhC. This is the equality of required and available delay and is not claimed as new.

## 3. Numerical Verification

The FDTD extraction was verified in a ladder (Figure 2). In vacuum and homogeneous dielectrics the extracted delay equals L/c and nL/c up to the closed-form Yee-grid correction, with error converging at second order (vacuum mean error 2.2e-3 a/c at 40 px/a for L = 10). An eight-period quarter-wave Bragg stack agrees with the transfer-matrix method; the median relative error in the complex transmission falls from 2.6e-2 to 1.4e-3 between 80 and 320 px/a (order about 1.9). A finite 2-D PhC stop band agrees with the parity-resolved MPB gap to within 3.5e-3 for every variant studied, dominated by finite-length and termination effects.

For the W1, the FDTD group index agrees with the MPB value over 0.231–0.268 with a median bias of +0.3 % and a maximum of 1.0 % (resolution 40, tolerance 1e-7). The default decay tolerance (1e-6) gave 5–14 % rms noise and was not used. The upper edge 0.268 of this validated band was set after inspecting the data (Section 12), so the band is domain-limited rather than a measured physical band edge. A time-domain pulse check agrees with the phase-derived delay to 0.4–0.9 %.

## 4. The W1 Delay Line

The guided H_z-even band (MPB band 8) has n_g,min = 3.87 at f = 0.2564. Over the validated band the plain W1 has a median delay of 3.932 a/c per period and an n_g within ±10 % of its median across 0.231–0.268 (14.8 % fractional bandwidth), with transmission above 0.99 and reflection below 0.3 %. The ±10 % criterion is a choice, and a group index near 4 is moderate for a slow-light W1 (compare n_g of order 30 in flat-band designs [Li et al. 2008]). Because both band ends coincide with validated-range limits, 14.8 % is a lower bound on flatness within the model; the same rule without the restriction gives 18.0 %, but that range is not validated.

First-row engineering was tested and rejected. Two variants (first-row radius ×1.2 with hole shifts −0.05 and −0.10) raise the median delay by 7 % and 14 % but reduce the flat-delay bandwidth to 8.4 % and 5.6 % and raise the group-delay dispersion by factors of 4 and 8 (Supplementary Figure S1). This is consistent with the known delay–bandwidth–dispersion trade [Frandsen et al. 2006; Li et al. 2008; Tucker et al. 2005]; three designs are an illustration, not a Pareto frontier. The plain W1 was therefore used downstream, and the delay line is not dispersion-engineered.

## 5. Multi-Channel Delay Family

Four sections of the same plain W1 with N = 8, 9, 10, 11 periods were simulated directly (Figure 3); channels were not built by cascading. The band-fit delay increment between adjacent channels is 3.9315, 3.9308 and 3.9310 a/c, within 0.02–0.04 % of the 3.9323 a/c target. This is a statement about the band fit, not about flatness: pointwise estimates of the adjacent increment vary by several percent across the validated band, roughly 3.88–4.18 a/c, following the rise of n_g toward the lower band edge. Transmission is 0.997–1.0025 (the excess over unity is numerical) and the amplitude imbalance among channels is at most 0.023 dB. A cascade prediction from shorter devices reproduces longer ones to 0.006 rad and 0.5 %.

The raw channels do not have the target phase. The per-period phase advance differs from 2πf0 τ_target by −1.2035 rad, so the frequency-independent intercept errors are 0, −1.20, −2.41 and +2.67 rad for the four channels. One constant trim per channel, chosen at f0, reduces the phase error to at most 0.026 rad (1.5°) over the band; the remainder is mostly the curvature of n_g(f). Trimming leaves amplitude and group delay unchanged. Because all four channels are one uniform medium at different lengths, these favorable numbers reflect one flat-n_g mode and not four independently engineered elements.

## 6. Element-Level Use: a Negative Result

We first use the four trimmed channels as element-level delays of a four-element array steered to 30° (Figure 4). One W1 period delays the carrier by 0.981 cycle. An array with λ0/2 element spacing steered to 30° needs only 0.25 cycle per element. The available step is therefore about four times too coarse, and using the channels element by element forces the element pitch to d = τ/sin30° = 7.86 a = 1.96 λ0. At that pitch every case has grating lobes as strong as the design lobe (about 0 dB), three to five equal lobes across the band. The result is not a usable aperture, even though the trimmed channels track ideal TTD to within 0.050° (R_squint 0.019) and fixed phase squints by +2.69° to −2.26° over the band (4.9° peak to peak, about 0.6 beamwidths). The untrimmed channels suppress the squint spread but point 6.3° low because of their intercept errors.

The problem is the mismatch between the available delay step and the element-level spatial sampling, not the PhC in isolation: the same trimmed channels follow ideal TTD in delay. Quantization-lobe mechanisms of this kind are standard array theory [Mailloux 2005‡]. We retain the result because it states plainly what a physically valid delay line can fail to do. A constant phase trim cannot repair it: element-level use with a constant trim (Msub = 1) cannot supply the missing 0.75 cycle of group delay, and the squint is worse than with fixed phase. The prior PhC delay-line work [Lin et al. 2012; Chung et al. 2018] avoids the issue by tuning the delay continuously with wavelength; a fixed-length section has no such freedom.

## 7. Subarray-Level Use: a Worked Near-Match Example

Each W1 channel can instead feed a subarray of Msub λ0/2-spaced elements with fixed intra-subarray steering phases, the PhC providing the inter-subarray delay (Figure 5). This is the standard subarray-TTD hierarchy; the present case is a worked example of applying the standard relation of Section 2 to the available W1 delay step.

At θ0 = 30° with Msub = 4 the required adjacent-subarray delay is 1.000 cycle. The W1 one-period delay is 0.981 cycle, a mismatch of −1.89 %. The exact analytic match for Msub = 4, k = 1 is sinθ = 0.4906, θ = 29.377°, so the 30° case is near, but not at, the exact match. For other subarray sizes the mismatch at 30° is +292 %, +96 % and −51 % for Msub = 1, 2 and 8 (one-period steps). Msub = 8 with two-period steps was not simulated. The near-match depends on n_g f0 ≈ 0.98 for this geometry and carrier together with a subarray of four; it is not presented as a design principle.

Results for Msub = 4 (16 elements, four PhC channels, subarray pitch 2 λ0), with a constant trim re-chosen for the required 4.008 a step (Figures 6 and 7):

* The maximum pointing error over 0.231–0.268 is 0.189° (rms 0.110°), against 2.69° for fixed phase (R_squint = 0.070). The tracking tolerance of 0.367° (5 % of the half-power beamwidth) is met. Ideal subarray TTD itself has 0.156° because the fixed intra-subarray phases squint; the PhC adds 0.055°.
* The 0 dB grating lobes are replaced by one design lobe and ordinary sidelobes (strongest other lobe −12.4 dB; ideal subarray TTD −12.4 dB). Subarray replica lobes at the 2 λ0 pitch remain at −16.6 dB or lower (growing from a null at f0), suppressed by the squinting λ0/2 element factor. This mitigates the grating-lobe problem; it does not eliminate it.
* The minimum peak power is 0.983 (−0.075 dB); the HPBW is 7.95°, 7.35° and 6.83° at the band low, center and high. Resolutions 30 and 40 agree to 0.007° and 0.04 dB.
* The raw PhC channels without trim point 6.7° off. The trim set for the 3.9323 a step (rather than 4.008 a) gives 0.73° error. A constant phase trim re-chosen for the architecture is required.

## 8. Dependence on Steering Angle

The sweep reuses the same four simulated channels at every angle and resets only the constant phase trim for each θ0. No new FDTD was performed. The sweep therefore tests array-level reuse of one validated channel family, not independent PhC physics at every angle. It runs from 5° to 70° in 1° steps (Figure 8). At 30° it reproduces the values of Section 7.

**Match table.** Table 1 lists the exact matching angles from sinθ_match = 2 k f0 τ_PhC/Msub = 1.96222 k/Msub. Only k = 1 channels (single-period increments) were simulated; entries with k > 1 are arithmetic extrapolations and assume that the Bloch-phase intercept of a longer section remains trimmable by one constant, which was not verified. Entries marked – have sinθ > 1.

| Msub | subarray pitch | k = 1 | k = 2 | k = 3 | k = 4 | k = 5 | k = 6 | k = 7 | k = 8 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.5 λ0 | none (sinθ = 1.96) | – | – | – | – | – | – | – |
| 2 | 1 λ0 | 78.85° | – | – | – | – | – | – | – |
| 4 | 2 λ0 | **29.38°** | 78.85° | – | – | – | – | – | – |
| 8 | 4 λ0 | 14.20° | 29.38° | 47.38° | 78.85° | – | – | – | – |
| 16 | 8 λ0 | 7.04° | 14.20° | 21.59° | 29.38° | 37.82° | 47.38° | 59.15° | 78.85° |

*Table 1. Match angles for Msub subarray elements and k PhC periods per adjacent-subarray increment (f0 τ_PhC = 0.98111). Only k = 1 channels were simulated.*

**Useful window.** The thresholds below were chosen after the 30° result was seen and are descriptive, not tests passed. Over 0.231–0.268:

| criterion (post hoc) | angular window |
|---|---|
| R_squint ≤ 0.20 | 24–36° |
| max difference from ideal subarray TTD ≤ 0.05 HPBW | 26–34° |
| max pointing error ≤ 0.5° | 22–34° |
| max pointing error ≤ 1° (for context) | 16–39° |
| all three stricter criteria | 26–34° |

The window edges are known to the 1° grid. The window is asymmetric about the exact match at 29.4° because the intra-subarray squint adds to the delay-mismatch term on one side and subtracts on the other. Outside it the error grows: the maximum pointing error is 0.049° at 28°, 0.49° at 34°, 1.04° at 40°, 2.4° at 52° and 7.8° at 70°, and toward small angles 0.63° at 20°, 1.08° at 14° and 1.76° at 5°. R_squint is below 1 over 14–70°, but with errors of several degrees at the large-angle end this is not meaningful scanning in terms of the beamwidth (about 7°). The strongest other lobe stays at −12.4 dB up to 52° and exceeds −10 dB from 53°; replica lobes reach 0 dB from 64°; the ideal-subarray reference also degrades at wide angles, so this is array geometry and not a PhC effect.

A first-order analytic estimate of the mismatch contribution, sinθ_pk − sinθ0 ≈ (f0/f − 1) sinθ0 (1 − ρ) with ρ = q_phc/q_req, agrees with the numerical peak-angle trajectory to within about 0.15° up to 40° (0.43° at 70°). It is a guide to the mismatch term, not a substitute for the full array evaluation near the match. Re-trimming the phase for each angle also adds a calibration step that the 30° robustness figures do not cover.

The reasonable conclusion is that the Msub = 4 channel family has a narrow useful angular window, about ±4° around 30° under the stricter criteria. The sweep does not show a broadly scannable architecture. Wider coverage would need additional delay increments or other (Msub, k) combinations from Table 1, which were not simulated.

## 9. Robustness and Tolerances

Part of the analysis perturbs the measured transfer functions mathematically (uniform independent errors per channel; Monte Carlo seed 20261010, 1000 samples per point; Supplementary Figure S2). These are system-level tolerance estimates, not fabrication tolerances; no fabrication-disorder FDTD was run, and slow-light disorder scattering is a separate, known problem [Hughes et al. 2005†; Mazoyer et al. 2009]. They were derived at 30° only. The 95th percentile of the maximum pointing error is 0.35° for a combined perturbation of ±2° phase, ±0.5 % delay and ±0.25 dB amplitude (nominal 0.19°). A P95 error of 0.5° is reached at ±7.6° of trim phase error and ±1.2 % uncalibrated delay error (relative to q·τ_PhC), and not within ±6 dB of amplitude error; if the intercept is re-trimmed after a delay change, ±10 % delay error stays below 0.5°. In absolute-length terms the uncalibrated limit is tight (about 0.14 % of the total device delay, or about 0.02 carrier cycle); a length error is benign only if the constant phase is re-trimmed afterward. Amplitude mismatch is not a pointing limiter in this range.

## 10. Length Compared with Propagation References

For one carrier cycle of delay at f0 = 0.2495 (λ0 = 4.008 a), a non-dispersive medium of group index n_g needs L = λ0/n_g. This is an analytic comparison of length at fixed delay only.

| medium | n_g | L/λ0 | L/a | length relative to W1 |
|---|---|---|---|---|
| air | 1 | 1.0000 | 4.008 | 3.93× |
| ideal homogeneous ε = 12 | √12 = 3.464 | 0.2887 | 1.157 | 1.135× |
| selected plain W1 | 3.932 | 0.2543 | 1.019 | 1.00 |

*Table 2. Length for one carrier cycle at f0. The ε = 12 row is an ideal homogeneous-medium propagation benchmark, not a modeled practical waveguide.*

The W1 is about 3.93 times shorter than air, which holds for any slow-wave medium, and about 11.9 % shorter than ideal ε = 12 propagation. A real dielectric waveguide has effective and group indices set by its geometry, and no conventional line was simulated. Width, cladding, loss, couplers and bends are not compared, and a conventional line can be cut to any length whereas the simulated W1 channels exist only in whole periods. At the modest selected group index, a large footprint advantage over a high-index dielectric propagation reference is not demonstrated. The experimentally demonstrated PhC delay lines cited above work at group indices of order 20–30, where the length advantage is larger [Lin et al. 2012; Chung et al. 2018].

## 11. Discussion

The question is what the study teaches, beyond a worked example.

1. **A valid delay line can be a poor element-level delay.** The simulated W1 sections are accurate delay elements: the band-fit increment is within 0.04 % of the target, and the trimmed channels track ideal TTD at element level to 0.05°. Yet their minimum step of 0.981 carrier cycle is four times the 0.25 cycle needed for half-wave steering at 30°, so they cannot be used element by element.
2. **Delay granularity and array geometry must be considered together.** The delay in carrier cycles per period is n_g f0, a property of the medium and carrier. Whether that step is usable is decided by the array: it must equal Msub sinθ0/2 for the subarray size and angle at hand (Section 2). This relation is standard; the computed channels provide a concrete case of the step against it.
3. **The same delay family can become usable at a different hierarchy level.** Four-element λ0/2 subarrays turn the same physical delay step from unusable to nearly matched near 30°. The sparse subarray pitch (2 λ0) leaves replica lobes at −16.6 dB or lower, so the improvement is mitigation of the grating-lobe problem and not its elimination.
4. **The match is local in angle.** The useful window for the Msub = 4 channels is about 26–34° under the stricter criteria. Table 1 shows that other matches need other (Msub, k) combinations, mostly with k > 1 channels that were not simulated.
5. **A modest group index does not show a footprint advantage.** The selected W1 is 12 % shorter than ideal ε = 12 propagation for the same delay. The footprint argument for PhC delay lines belongs to much higher group indices, which come with the bandwidth, dispersion and disorder penalties noted above [Tucker et al. 2005; Krauss 2007; Mazoyer et al. 2009].
6. **Full-wave-to-system modeling exposes constraints that ideal TTD models can hide.** An ideal-delay model contains neither the whole-period quantum nor the Bloch-phase intercept (up to 2.7 rad here) that a constant trim has to remove.

*Carrier mapping.* Normalized results can be scaled when the PhC delays the radiated carrier (for example a = 2.49 mm at 30 GHz, 0.80 mm at 94 GHz, 387 nm at 1550 nm). They do not describe RF-photonic beamforming, where the optical device delays an RF-modulated optical carrier and the relevant delay is in RF cycles; the PhC delay lines of Lin et al. and Chung et al. are of that type, and their delay increments are continuously wavelength-tunable. The one-cycle-per-period argument is specific to the carrier-delay mapping. The 2-D structure is not claimed to be a microwave feed network.

## 12. Limitations

* 2-D simulation (TE, no out-of-plane radiation); lossless, non-dispersive dielectric; no propagation loss; no fabrication disorder; no bends; no input/output couplers; no channel coupling; no thermal effects; isotropic elements with no real element pattern and no mutual coupling; no experiment.
* Only four channels, all one flat-n_g mode at four lengths, and only one simulated delay increment (one period). Entries of Table 1 with k > 1 and all Msub = 8 and 16 hierarchies are arithmetic only.
* Validated band 0.231–0.268 c/a (14.8 %). The upper edge was chosen after inspecting the data: the projection residual, unitarity defect and n_g bias grow above about 0.265–0.268, and the cause was only partly identified. 14.8 % is a domain-limited lower bound, and the delay is not flat to better than several percent pointwise within it.
* One constant phase trim per channel (raw intercept errors up to 2.7 rad) is required; how it would be realized physically is not addressed, and it must be reset for each steering angle.
* The angle sweep reuses the same channel family at every angle. The useful window is narrow. Its thresholds are post hoc; the 1° grid limits the edge resolution; the resolution-40 check was done only at 30°.
* The robustness numbers are system-level Monte Carlo estimates on the measured transfer functions, derived at 30° only; they are not fabrication tolerances.
* Residual subarray replica lobes (strongest −16.6 dB); the count of visible replicas depends on whether grazing lobes at ±90° are counted. The fixed intra-subarray phase leaves a residual squint of about ±0.16° even for ideal TTD.
* No advantage over conventional delay lines is demonstrated; the ε = 12 reference in Table 2 is an ideal analytic medium and no conventional line was simulated.
* Two-resolution (30 and 40) convergence only, not a full series. Pointwise group-delay derivatives and local group-delay dispersion are noise-limited.
* The literature search was limited (Section 1); Burla et al. was read at abstract level only, and Mailloux was not read.

**Disclosures about analysis choices.** The validated band boundary was set after inspecting the data. The first attempt at the four-channel simulation was discarded as under-converged (erratic stop times; one channel with a 1.9 % unitarity defect); the reported channels come from a repeated run with a common minimum run time, which was added after seeing the data. The pre-declared pointwise central-difference delay-increment estimator was replaced by band fits and phase slopes after it proved noise-limited, and all delay claims rest on the latter. Several diagnostics (replica-lobe level, efficiency loss, calibrated absolute-length delay) and the angular-window thresholds were defined after seeing results and are labeled post hoc. Earlier code defects found and fixed during development (with regression tests) and the development chronology are recorded in the repository documentation and are not material to the results reported.

## 13. Conclusion

In the model studied here, the W1 channels work at the subarray level only near a matched operating angle, while the same physical delay step is poorly suited to element-level half-wave steering. In a 2-D lossless model, a plain W1 gives about one carrier cycle of delay per period. Used element by element this forces a 1.96 λ0 pitch and 0 dB grating lobes. Used for four-element λ0/2 subarrays near 30°, with one constant phase trim per channel, the maximum pointing error over the validated 14.8 % band is 0.19° against 2.69° for fixed phase, with a useful window of about 26–34°. The calculation therefore isolates a delay-granularity and hierarchy effect; it does not show a generally superior or broadly scannable PhC beamformer, nor a large footprint advantage over a high-index propagation reference.

## References

Entries marked † have details (volume, pages or author list) still to be confirmed against the publisher record. ‡ marks a source whose content was not read.

1. M. Notomi, K. Yamada, A. Shinya, J. Takahashi, C. Takahashi, I. Yokohama, "Extremely large group-velocity dispersion of line-defect waveguides in photonic crystal slabs," Phys. Rev. Lett. 87, 253902 (2001).
2. Y. A. Vlasov, M. O'Boyle, H. F. Hamann, S. J. McNab, "Active control of slow light on a chip with photonic crystal waveguides," Nature 438, 65–69 (2005).
3. L. H. Frandsen, A. V. Lavrinenko, J. Fage-Pedersen, P. I. Borel, "Photonic crystal waveguides with semi-slow light and tailored dispersion properties," Opt. Express 14, 9444–9450 (2006).
4. J. Li, T. P. White, L. O'Faolain, A. Gomez-Iglesias, T. F. Krauss, "Systematic design of flat band slow light in photonic crystal waveguides," Opt. Express 16, 6227–6232 (2008).
5. D. Mori, S. Kubo, H. Sasaki, T. Baba, "Experimental demonstration of wideband dispersion-compensated slow light by a chirped photonic crystal directional coupler," Opt. Express 15, 5264–5270 (2007).
6. M. L. Povinelli, S. G. Johnson, J. D. Joannopoulos, "Slow-light, band-edge waveguides for tunable time delays," Opt. Express 13, 7145–7159 (2005).
7. T. F. Krauss, "Slow light in photonic crystal waveguides," J. Phys. D: Appl. Phys. 40, 2666–2670 (2007).
8. T. Baba, "Slow light in photonic crystals," Nat. Photonics 2, 465–473 (2008), doi:10.1038/nphoton.2008.146 †.
9. R. S. Tucker, P.-C. Ku, C. J. Chang-Hasnain, "Slow-light optical buffers: capabilities and fundamental limitations," J. Lightwave Technol. 23, 4046–4066 (2005).
10. E. Hughes, L. Ramunno, J. F. Young, J. E. Sipe, "Extrinsic optical scattering loss in photonic crystal waveguides: role of fabrication disorder and photon group velocity," Phys. Rev. Lett. 94, 033903 (2005) †.
11. S. Mazoyer, J. P. Hugonin, P. Lalanne, "Disorder-induced multiple scattering in photonic-crystal waveguides," Phys. Rev. Lett. 103, 063903 (2009).
12. C.-Y. Lin, H. Subbaraman, A. Hosseini, A. X. Wang, L. Zhu, R. T. Chen, "Silicon nanomembrane based photonic crystal waveguide array for wavelength-tunable true-time-delay lines," Appl. Phys. Lett. 101, 051101 (2012), doi:10.1063/1.4739003.
13. C.-J. Chung, X. Xu, G. Wang, Z. Pan, R. T. Chen, "On-chip optical true time delay lines featuring one-dimensional fishbone photonic crystal waveguide," Appl. Phys. Lett. 112, 071104 (2018), doi:10.1063/1.5006188.
14. M. Burla, D. Marpaung, L. Zhuang, A. Leinse, M. Hoekman, R. Heideman, C. Roeloffzen, "Integrated photonic Ku-band beamformer chip with continuous amplitude and delay control," IEEE Photon. Technol. Lett. 25(12), 1145–1148 (2013), doi:10.1109/LPT.2013.2257723 (abstract-level reading only).
15. J. Capmany, D. Novak, "Microwave photonics combines two worlds," Nat. Photonics 1, 319–330 (2007), doi:10.1038/nphoton.2007.89.
16. J. Yao, "Microwave photonics," J. Lightwave Technol. 27, 314–335 (2009).
17. R. J. Mailloux, *Phased Array Antenna Handbook*, 2nd ed., Artech House (2005) ‡ (cited for the general subarray-TTD and quantization-lobe topics only; edition details and the form of the relations were not verified).
18. D. Q. Nguyen, T. Kim, "Joint delay-phase precoding under true-time delay constraints in wideband sub-THz hybrid massive MIMO systems," arXiv:2212.07484 (v7, 2024); journal venue not confirmed.
19. L. Dai, J. Tan, H. V. Poor, "Delay-phase precoding for wideband THz massive MIMO," IEEE Trans. Wireless Commun. (2022), doi:10.1109/TWC.2022.3157315, arXiv:2102.05211 (background; not read closely).
20. A. F. Oskooi, D. Roundy, M. Ibanescu, P. Bermel, J. D. Joannopoulos, S. G. Johnson, "Meep: A flexible free-software package for electromagnetic simulations by the FDTD method," Comput. Phys. Commun. 181, 687–702 (2010).
21. S. G. Johnson, J. D. Joannopoulos, "Block-iterative frequency-domain methods for Maxwell's equations in a planewave basis," Opt. Express 8, 173–190 (2001), doi:10.1364/OE.8.000173.

`[CITATION NEEDED]`: a verified primary source for beam squint under fixed phase shifters (cited here only through the general references [15, 16, 17]); the Yee FDTD paper and standard PhC/FDTD textbooks (not re-verified); a verified source for closed-form replica-lobe levels of delay-quantized subarrays.

## Figures

Main text (file names under `figures/final/`):

1. Figure 1 — `fig01_concept`: W1 channel and the element-level and subarray architectures (2-D FDTD domain versus analytic beamforming model).
2. Figure 2 — `fig02_verification`: FDTD verification (convergence, Bragg versus transfer-matrix, W1 n_g versus MPB).
3. Figure 3 — `fig04_channels`: four-channel transfer functions and delay family.
4. Figure 4 — `fig05_phase5_limit`: element-level array and grating lobes (negative result).
5. Figure 5 — `fig06_architecture`: required-versus-available delay, quantization and sidelobe level versus Msub.
6. Figure 6 — `fig07_patterns`: Msub = 4 patterns at band low, center and high.
7. Figure 7 — `fig08_pointing`: pointing versus frequency.
8. Figure 8 — `fig10_angle_sweep`: dependence on steering angle.

Supplementary:

* Figure S1 — `fig03_tradeoff`: delay–bandwidth–dispersion comparison of the plain and first-row-engineered W1.
* Figure S2 — `fig09_robustness`: system-level perturbation sensitivity.
* Figure S3 — `fig11_match_map`: match-angle map for (Msub, k) (graphical form of Table 1).
