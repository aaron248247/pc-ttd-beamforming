# Phase 4 report: four-channel true-time-delay family (plain W1, D0)

> **Status note (Phase 7).** The status line below is historical. Phases 5 and 6 were completed afterwards; Phases 5 and 6 contain the later array-level work. The Phase 4 conclusions are unchanged.


Status: **gate PASS WITH RESTRICTIONS** (2026-10-10). Phase 5 has not been started and is not authorized until this report is reviewed.
Tables: [phase4_report_tables.md](phase4_report_tables.md) (generated). Figures: `figures/phase4/`. Definitions: `docs/phase2_5_criteria.md`, Revisions 4 and 5.

## 1. Objective

Determine whether four directly simulated finite sections of the Phase 3 candidate D0 (plain W1, r/a = 0.30, eps = 12, TE, Hz-even guided band 8), with systematically increasing length,
form a usable monotonic multi-channel delay family, and deliver their **actual complex transfer functions H_m(f)** for the Phase 5 array model.

## 2. Phase 3 starting point

Checkpoint 457825d (142 tests passing): D0 chosen because D1/D2 gain only 7%/14% in delay per period but lose 43%/62% of the +-10% usable bandwidth and are 4-8 times more dispersive. D0: validated usable band 0.231-0.268 (14.8%),
median per-period delay 3.93 a/c. Carried-forward restrictions (all respected): FDTD claims only inside 0.231 <= f <= 0.268 (a post hoc Phase 2 validation limit, not a physics band edge); decay tolerance <= 1e-7;
the negative-s geometry fix and its regression test are untouched; the custom MPB-derived guided-mode projection is used (no `get_eigenmode_coefficients`); pointwise FDTD derivatives (hence GDD) are noise-limited.

## 3. Channel-family definition

One structure, four lengths N_m = N0 + m dN with **N0 = 8, dN = 1**: N = 8, 9, 10, 11. Each channel is a separate full FDTD run (no cascading to build channels). Reference planes x = -3 and x = N_m + 3, source d_src = 10,
lead = plain W1 (identical to the D0 section: the device is uniform with no interface). Production res 30; res 40 confirmation of **all four** channels. Frequency grid 0.220-0.275, step 5e-4.

## 4. N0 and dN selection

Declared on 2026-10-10 in criteria Revision 4 (commit 7a24234) **before any Phase 4 FDTD result existed**; the repository did not define them. dN = 1 is the smallest integer-period step; N0 = 8 lies inside the already validated
length range and adds no unnecessary length. Target step dtau_target = 1 x 3.9323 a/c (the Phase 3 D0 median per-period delay over the validated band, res 30; res 40 gives 3.9308); trimming frequency f0 = 0.2495 (center of the
Phase 3 band). Neither was refit to Phase 4 data. One period gives f0 x 3.9323 = 0.981 carrier cycle, a coarse delay quantum that is characterized, not hidden; no sub-period interpolation is performed.

## 5. Numerical method

Same as Phase 3 (MEEP, DFT planes (Ey, Hz), modal projection with the MPB lead mode and Lorentz-reciprocity overlap), `decay_tol = 1e-7`. **A method change was required (criteria Revision 5):**
the first attempt, with `stop_when_dft_decayed` alone, stopped at erratic, non-monotone times (res 30: t_end = 1965, 801, 2816, 2462 for N = 8-11; res 40: 636-1043 for N = 9-11) and N = 9 was visibly under-converged
(unitarity defect 1.9%, T up to 1.018, cascade phase mismatch above 1 rad). That attempt was discarded (kept in `results/phase4/raw_attempt1/` and `attempt1_*`; its H_m differ from the final ones by <= 0.01 rad, so the phase itself was robust, but the amplitude/unitarity diagnostics were not) and every channel
was rerun with a **common minimum run time of 3500 after the source** (stop times 3524-4586 at res 30, 3582-4965 at res 40). With it the unitarity defect is <= 0.34% (res 30) / 0.40% (res 40).
This is a numerical-method change motivated by a convergence diagnostic, not a relaxation of any acceptance criterion.

## 6. Complex transfer functions

H_m(f) = S21_m = a_f(out)/a_f(in) between x = -3 and x = N_m + 3. It contains the **entire** propagation phase over N_m + 6 periods; nothing is de-embedded and the reference-plane offset is common to all channels, so it cancels in H_m/H_0.
Convention: MEEP exp(-i omega t), MEEP_SIGN = +1, a delay tau gives exp(+i omega tau), tau_g = +d arg H / d omega (verified: all channel delays are positive and increase with N at every sample in the validated domain). Phase branch fixed by the MPB Bloch phase
(consistent with a plain unwrap to 4e-15 rad). Saved per channel: complex S21, S11, T, R, projection residuals (in/out), backward amplitude at the output plane, unitarity defect, unwrapped phase, delay, stop time, configuration and provenance
(`results/phase4/4b_res{30,40}.{npz,json}`).

## 7. Transmission and reflection (Fig. 4a-4c)

In the validated domain (res 30) T = 0.997-1.0025 (the excess over 1 is numerical, <= 0.25%), R <= 0.25%, projection residual <= 0.054 (0.048 at res 40), |unitarity defect| <= 0.0034. There is no systematic dependence on length:
|H_m(f0)| = 0.9997, 1.0002, 0.9997, 1.0002 for N = 8-11. The structure is lossless in the model and has no interfaces, so loss-driven amplitude growth with length cannot appear here (see Limitations).

## 8. Relative phase (Fig. 4d)

arg[H_m/H_0] is smooth and increases with frequency. The raw channels do **not** realize the target phase: the per-period phase advance at f0 differs from 2 pi f0 dtau_target by -1.2035 rad (-0.192 cycle, -69 degrees), a physical consequence of the Bloch phase
of the W1 mode, not an error of analysis. The relative-phase intercept errors at f0 are therefore delta_m = 0, -1.2035, -2.4075, +2.6725 rad (= -3.6107 + 2 pi) for m = 0-3 (res 40: -1.2046, -2.4085, +2.6688).

## 9. Incremental delay (Figs. 4e-4g)

All delay statements below use the local-quadratic (+-3 samples) estimator or least-squares phase slopes; see Section 16 for why the pre-declared central difference is not used as the primary estimator.
Over the common band, least-squares band-fit adjacent increments are 3.9315, 3.9308, 3.9310 a/c (errors -0.02%, -0.04%, -0.03% against the 3.9323 target; res 40: 3.9270, 3.9257, 3.9261, -0.14% to -0.17%).
At f0 the local-quadratic increments are 3.907, 3.870, 3.935 (-0.65%, -1.60%, +0.08%). Pointwise inside the band the increments stay within +-9.7% of the target (res 40: +-9.0%), within +-3.9% only for f >= 0.238
(11.9% band). The upward excursion at the low end follows the MPB n_g rise toward the band edge (4.19 at 0.231). Relative-to-channel-0 delays Dtau_m follow m dtau_target to the same accuracy (max 8.5%, 7.1%, 8.0% pointwise for m = 1, 2, 3; band-fit 3.932, 7.862, 11.793 a/c against 3.932, 7.865, 11.797).
The family is monotonic: absolute delay increases with N at every sample of the validated domain.

## 10. Raw phase error (Fig. 4h)

epsilon_phi_m = arg(H_m/H_0) - 2 pi f m dtau_target (branch with epsilon(f0) in (-pi, pi]). Maximum over the common band: 1.212, 2.425, 2.676 rad for m = 1, 2, 3. It is almost entirely a frequency-independent intercept
(m times -1.2035 rad); its variation across the band equals the trimmed error below.

## 11. Phase-intercept trimming (Fig. 4i)

Raw PhC: H_m as simulated. Phase-intercept-trimmed PhC: H_m exp(-i delta_m), one constant per channel chosen at f0 = 0.2495 only; amplitude and group delay are unchanged by construction (unit-tested), and the error is zero at f0 to numerical precision.
Residual trimmed phase error over the common band: max 0.0086, 0.0172, 0.0260 rad (0.5, 1.0, 1.5 degrees) for m = 1, 2, 3 (res 40: 0.0085, 0.0171, 0.0245), i.e. about 15 times below the pi/8 reporting threshold. Removing the best-fit line in f leaves 0.0069, 0.0143, 0.0210 rad,
so the residual is mostly curvature of n_g(f), not a delay-slope error. Trimming therefore corrects an intercept that is large (up to 153 degrees), and what remains is physical dispersion, not hidden. The trim is a constant phase that a real implementation would have to supply
(about 0.19 cycle per period of length here); this study does not show how.

## 12. Amplitude mismatch (Fig. 4j)

Amplitude (not power) imbalance 20 log10(max|H_m|/min|H_m|): max 0.023 dB in the band, 0.0046 dB at f0; |H_m|/|H_0| in [0.9974, 1.0026]. No systematic growth with length. This is consistent with the numerical noise level (T within +-0.25% of 1).

## 13. Common usable band (Fig. 4k)

Primary band (criteria Revisions 4 and 5): **0.2310-0.2680, center 0.2495, fractional bandwidth 14.8%**, identical at res 30 and res 40. It is **domain-limited**: every other criterion (T >= 0.5, |defect| <= 0.02, residual <= 0.10, imbalance <= 1 dB) is satisfied over the entire
analysis grid or beyond the domain, and the increment criterion alone (local quadratic) would extend it to 0.2305-0.2725. Because 0.268 was itself a post hoc Phase 2 limit, this band is a validation-domain statement, not a measured physical limit. Sensitivity: +-5% increments give 0.2380-0.2680 (11.9%; res 40: 0.2365-0.2680).
The pre-declared central-difference rule gives a fragmented band, widest run 0.2425-0.2550 (5.0%) at res 30 and 0.2405-0.2680 (10.8%) at res 40: the result changes with resolution because it follows noise, see Section 16. The Phase 5 center frequency, the band center, equals f0, so no secondary trim is needed.

## 14. Cascade validation (Fig. 4l)

T-matrix prediction from directly simulated shorter devices vs the directly simulated longer device: (8, 9) -> 11 (2 steps), (9, 10) -> 11 (1 step), (8, 9) -> 10 (1 step). Max phase difference 0.0058, 0.0033, 0.0043 rad and max relative |S21| difference 0.0036, 0.0041, 0.0020 over the band (res 40: <= 0.0049 rad, <= 0.0047);
all accepted (<= 0.1 rad, <= 3%). The agreement holds over the whole validated domain, with no frequency region where the cascade becomes unreliable. The direct channels were not replaced by cascaded ones. Caveat: D0 is uniform (no interface),
so this verifies the machinery and the direct runs' mutual consistency, not interface physics. (An initial cascade comparison showing 1.5 rad was my own step-count error, criteria Revision 5 item 3, plus the under-converged N = 9.)

## 15. Resolution validation

Res 30 vs res 40, all four channels, the four pre-declared checks (Revision 4 item 8): band-mean increments change +0.13%, +0.12%, +0.15% of target (<= 2%); band edges identical; imbalance 0.023 vs 0.022 dB; maximum trimmed phase error 0.0260 vs 0.0245 rad. All pass.
Max phase difference of H_m/H_0 between resolutions in the validated domain 0.0058 rad, max relative |H| difference 0.0021. Conclusions on increments, band, amplitude mismatch and phase error do not change.

## 16. Phase 4 gate

1. All four channels directly simulated and numerically trustworthy: **yes**, after the common-run-time correction (Section 5); res 40 agrees to 0.006 rad / 0.2%.
2. Monotonic delay with N: **yes**, at every sample of the validated domain.
3. Adjacent steps near the target: **yes**: band-fit -0.02% to -0.04% (res 40: -0.14% to -0.17%); pointwise within +-9.7% (local quadratic); +-3.9% for f >= 0.238.
4. Nontrivial common band: **yes**, 0.231-0.268 (14.8%), domain-limited.
5. Amplitude differences small enough: **yes**, <= 0.023 dB.
6. Raw phase errors characterized: **yes**; large and almost purely an intercept (up to 2.68 rad).
7. Trimming improves broadband consistency without hiding delay errors: **yes**, residual <= 0.026 rad, mostly n_g curvature, and the delay slope is unchanged by trimming.
8. Direct vs cascade: **agrees** (<= 0.006 rad, <= 0.5%).
9. H_m(f) trustworthy for Phase 5: **yes within the stated restrictions**.

**Outcome: PASS WITH RESTRICTIONS.** Restrictions that travel with the handoff:
(a) quantitative use only inside the validated domain 0.231-0.268 (post hoc limit); (b) the band is domain-limited, not a measured band edge;
(c) the pre-declared pointwise central-difference increment criterion is noise-limited (5.5-7.5% rms against the smooth estimator at one-period spacing) and does not give a stable band; the primary estimator was changed to the local quadratic
**after** seeing data (Revision 5), with the pre-declared result reported next to it; delay claims rest on phase slopes and phase errors; (d) the raw family needs a per-channel constant phase (trim) to match the target phase at f0 (intercept errors up to 2.68 rad); the trimmed family must not be called ideal TTD;
(e) the delay quantum is 0.98 carrier cycle per period; (f) the model is lossless and perfectly uniform; (g) the first FDTD attempt was discarded for under-convergence (disclosed above).

## 17. Scientific interpretation

Because the four channels are the same uniform medium at different lengths, H_m/H_0 is essentially exp(i m theta(f)) with theta the per-period phase advance: the whole family is fixed by one flat-n_g mode, which is why the increments are equal to <1%, amplitudes match to 0.02 dB and the trimmed phase error is below 1.5 degrees across 14.8%.
The favorable numbers therefore reflect Phase 3's flat-n_g result, not independent channel engineering, and they say nothing yet about beam squint. The unfavorable facts are the coarse quantum (about one carrier cycle per period) and the large intercept error that a
real design would have to supply; neither is solved by Phase 4. No claim about beamforming is made.

## 18. Limitations

Lossless 2-D TE model; no fabrication disorder, no input/output coupling structure, no loss or scattering growth with length, no thermal or bend effects; single geometry (D0); frequency range limited to the post hoc validated domain;
four lengths of one period spacing, so a mid-band spacing/loss trend with larger N cannot be assessed; pointwise delay (and GDD) from H are noise-limited; the cascade check does not test interfaces; res 30 and 40 are not a full convergence series (two points), though they agree to 0.006 rad.

## 19. What Phase 5 is allowed to consume

`results/phase4/phase4_handoff.npz` (+ `.json` sidecar): common frequency grid; raw complex H_m (4 x 111) with amplitude, wrapped, absolute and relative unwrapped phase; trimmed H_m and trim constants; S11, T, R, residuals, unitarity defect;
`mask_common` (0.231-0.268), `mask_validated_domain`; f0, dtau_target, N_m; res 40 confirmation H arrays; sign convention, definitions and provenance in the JSON string. Phase 5 must use these H_m directly (raw and trimmed, always labeled), restrict quantitative conclusions to `mask_common`, and
state the delay-granularity and trimming restrictions above. Phase 5 requires explicit authorization after review of this report.
