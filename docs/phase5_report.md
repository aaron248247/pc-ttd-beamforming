# Phase 5 report — four-element array with the actual Phase 4 transfer functions

> **Update after Phase 6 (added in Phase 7).** The grating-lobe limitation reported here (equal-strength lobes at d = 1.96 lambda0) is mitigated at subarray level in Phase 6 (Msub = 4: strongest other lobe -12.4 dB, replica lobes <= -16.6 dB). The Phase 5 numbers and conclusions below are unchanged and still describe the element-level array.


**Gate result: PASS WITH RESTRICTIONS.** Criteria: `docs/phase2_5_criteria.md`, Revision 6 (declared and committed, 6553a1d, before any PhC array result was computed). Tables: `docs/phase5_report_tables.md` (generated). Figures: `figures/phase5/`. Results: `results/phase5/5a_array.{npz,json}` and `results/phase5/phase5_handoff.{npz,json}`.
Code: `src/pcttd/array.py` (unchanged), `src/pcttd/beam.py` (new), `experiments/phase5/5a_array.py`, `5b_figures.py`, `5c_tables.py`, `5d_handoff.py`; tests `tests/test_array.py` (extended), `tests/test_phase5.py`.

## 1. What was computed

The directly simulated complex transfer functions H_m(f) of the four D0 plain-W1 channels (N = 8, 9, 10, 11) from `results/phase4/phase4_handoff.npz` were used as array weights, with no fitted or substituted delays. Four cases, one four-element uniform linear array of isotropic elements, steered to theta0 = 30 deg:

| case | weights |
|---|---|
| fixed phase | exp(+i 2 pi f0 m dtau_target), frequency independent |
| ideal TTD | exp(+i 2 pi f m dtau_target) |
| raw PhC | handoff `H_raw` (amplitude and phase as simulated) |
| phase-intercept-trimmed PhC | handoff `H_trim` (one constant phase per channel, chosen at f0; amplitude, group delay, dispersion untouched). It is **not** ideal TTD. |

Convention: MEEP exp(-i omega t); AF(theta, f) = sum_m H_m exp(-i 2 pi f x_m sin theta), x_m = m d. Normalized power |AF|^2/16.

**Array parameters.** dtau_target = 3.9323 a/c, f0 = 0.2495 (lambda0 = 4.008 a). The pitch is set by the delay quantum, not by lambda0/2: **d = dtau_target / sin(30 deg) = 7.8646 a = 1.962 lambda0** (1.82-2.11 lambda over the band). One lattice period is 0.981 carrier cycle of delay, while a half-wave array at 30 deg needs 0.25 cycle. The array is therefore **sparse by construction and has grating lobes**. Design-lobe tracking window |sin(theta) - 0.5| <= 1/(2 f_max d) = 0.2372 (theta 15.2-47.5 deg). All primary numbers are over mask_common, 0.231-0.268 (75 samples).

## 2. Gate questions

1. **Fixed phase reproduces the analytic squint: yes.** Peak of the array factor vs sin(theta) = (f0/f) sin(theta0): maximum deviation 6e-12 deg. Pointing runs from +2.69 deg (0.231) through 0 (f0) to -2.26 deg (0.268), a 4.94 deg peak-to-peak squint (about 0.6 HPBW).
2. **Ideal TTD stays at theta0: yes**, to 4e-12 deg, peak power exactly 1 (unit-tested).
3. **Raw PhC reduces frequency-dependent squint: yes.** theta_peak spans 23.19-24.14 deg (0.95 deg peak-to-peak), against 4.94 deg for fixed phase. Offset-removed spread ratio 0.199.
4. **Raw PhC has an absolute pointing offset: yes, large.** Its lobe sits at 23.73 deg at f0, 6.27 deg (0.8 HPBW) below the target, -6.81 to -5.86 deg over the band. Its R_squint as defined is 2.53 because that metric includes the offset. The offset follows directly from the Phase 4 intercept errors: they are linear in the element index (0, -1.204, -2.408, -3.612 rad, the last equal to +2.67 rad wrapped), i.e. a slope of -1.204 rad per element, which the linear model turns into a shift of up to 7.0 deg against 6.81 deg numerically. Because the intercept is a *frequency-independent* phase, it behaves like a fixed phase shifter and itself squints (about 1/f); the raw residual spread is consistent with that effect (the linear model gives a shift that falls roughly as 1/f). **Low squint and correct absolute pointing are different properties: raw PhC has good delay slope and a wrong carrier phase.**
5. **The single trim places the lobe at theta0 and keeps the delay behavior: yes.** theta_peak = 29.950-30.006 deg. Trimming changes the phase intercept only (verified by test: |H| and the phase slope are identical to raw to 1e-8).
6. **Trimmed PhC vs ideal TTD.** Max |theta_trim - theta_ideal| = 0.050 deg (RMS 0.015 deg; 0.65% of the ideal HPBW; tolerance 0.384 deg), R_squint = 0.0186 (below the 0.20 stretch target, which is not a criterion). The Phase 4 residual phase error (max 0.026 rad) predicts 0.0501 deg with the phase-slope model (difference 1e-5 deg), so a 0.026 rad residual is indeed small enough that the pointing is indistinguishable from ideal TTD within the common band at this resolution of the analysis.
7. **Amplitude mismatch: negligible.** Replacing |H_m| by 1 (same phases) changes theta_peak by <= 2e-6 deg and peak power by <= 0.0022 (0.0095 dB), for a Phase 4 imbalance of <= 0.023 dB. The raw peak power ranges 0.9979-1.0002 (|H| can slightly exceed 1 at the band edge because of the residual unitarity defect, <= 0.34%).
8. **Stable under res 30 vs res 40: yes.** All predeclared checks met for both PhC cases (Table F): theta_peak differs by <= 0.010 deg, max-error and trimmed-vs-ideal values by <= 0.004 deg, R_squint by <= 0.0013, peak power by <= 0.002, grating-lobe level by <= 3e-6 dB. Trimmed R_squint 0.0186 vs 0.0177; max trimmed-vs-ideal 0.050 vs 0.048 deg.
9. **Grating lobes are severe: yes, by construction.** Every visible principal lobe has the same level as the design lobe (strongest other lobe 0.000 dB re the design peak in all four cases, all frequencies, because isotropic elements and equal weights give equal grating lobes). For ideal TTD the visible replicas are at -31.3 and -0.55 deg at f0 (-36.9, -2.9 deg at 0.231; -67.4, -26.7, +1.5, +77.0 deg at 0.268), so 3-5 equal lobes appear across the band. Near-endfire radiation also stays within about 0.1-0.7 dB of the peak at +-90 deg at f0 (grating skirts). The design lobe remains trackable over the whole band in all four cases (the farthest excursion, raw PhC, uses 45% of the window). **The PhC TTD suppresses squint of the m = 0 lobe; it does not suppress grating lobes, and this array is not usable as a standalone conventional aperture.**
10. **Scientifically promising for squint suppression: yes, with the restrictions below.** In this lossless 2-D model the actual four simulated channels reproduce ideal TTD pointing to 0.05 deg once the constant phase intercept is trimmed, over 14.8% fractional bandwidth. The result is a mixed one: the one-period delay quantum (0.981 carrier cycle) forces d of about 2 lambda0.

## 3. Key numbers (res 30, common band, degrees unless stated)

| case | max abs pointing error | RMS | error at 0.231 / f0 / 0.268 | R_squint | max vs ideal TTD | HPBW at 0.231 / f0 / 0.268 |
|---|---|---|---|---|---|---|
| fixed phase | 2.687 | 1.446 | +2.687 / 0 / -2.258 | 1.000 | 2.687 | 8.55 / 7.69 / 7.00 |
| ideal TTD | 4e-12 | 1e-12 | 0 / 0 / 0 | 0 | 0 | 8.31 / 7.69 / 7.16 |
| raw PhC | 6.809 | 6.299 | -6.809 / -6.273 / -5.865 | 2.534 (offset-removed 0.199) | 6.809 | 7.82 / 7.27 / 6.79 |
| phase-intercept-trimmed PhC | 0.0501 | 0.0150 | -0.0501 / 0 / -0.0150 | 0.0186 | 0.0501 | 8.30 / 7.69 / 7.16 |

HPBW is local to the design lobe (the beam broadens at 30 deg by 1/cos(theta); the raw lobe at about 23.7 deg is narrower). See the generated tables for gain, strongest-lobe, resolution and amplitude tables (A-H).

## 4. Distinctions kept explicit

* *True time delay*: frequency-dependent phase slope for a constant delay. *Phase intercept*: constant phase at f0. *Beam squint*: frequency dependence of the design-lobe angle. *Absolute pointing error*: lobe not at theta0 even if squint is small. *Grating lobe*: extra principal maximum because d > lambda/2.
* The array simultaneously has near-ideal TTD slope, low squint, a wrong absolute angle before trimming, and severe grating lobes.
* Design-lobe metrics (windowed) and full-pattern metrics (-90..+90 deg) are computed separately; the full patterns are shown in `figures/phase5/5d_patterns_full.png`.

## 5. Sparse-array discussion and routes not taken here (not implemented, not Phase 6 work)

The one-period quantum forces d = 1.96 lambda0 at 30 deg. Ways around it that remain to be studied: sub-period delay engineering (smaller delay increment per added period, or partial-period length steps), a lower delay per lattice period, subarray-level TTD (one PhC channel feeding a subarray of lambda/2 elements, so the TTD pitch is the subarray pitch and the grating lobes are suppressed by the subarray element pattern), hybrid phase + TTD, and an RF-photonic architecture where one photonic delay channel feeds a subarray. None is evaluated here.

## 6. Physical mapping (Table H)

f_normalized = a/lambda0 = 0.2495, lambda0 = 4.008 a. The one-period delay and channel increment are 3.9323 a/c; d = 7.8646 a; the N = 8-11 sections are 8-11 a long.

| scaling | a | one-period delay = dtau | array pitch d | N = 8-11 section |
|---|---|---|---|---|
| microwave, carrier 30 GHz (lambda0 = 9.99 mm) | 2.49 mm | 32.7 ps | 19.6 mm | 19.9-27.4 mm |
| mm-wave, carrier 94 GHz (lambda0 = 3.19 mm) | 0.796 mm | 10.4 ps | 6.26 mm | 6.4-8.8 mm |
| optical phased array, carrier 1550 nm | 387 nm | 5.07 fs | 3.04 um | 3.1-4.3 um |

Photonic carrier vs RF beamforming frequency: in an RF-photonic system the beamforming phase is set by the delay of the modulation envelope, not the optical-carrier phase. At 1550 nm one period delays the envelope by 5.07 fs, so a 30 GHz half-wave RF array (d = 5 mm, 30 deg) would need 8.34 ps per element step, about 1643 periods (635 um) per step if the group delay per period persisted for such lengths (not demonstrated). The Phase 5 one-period matching is therefore an optical-carrier phased-array scaling (beamformed wave = PhC carrier), not an RF-photonic design. The 2-D dielectric model is not a manufacturable RF-photonic beamforming system.

## 7. Restrictions (all Phase 4 restrictions travel with every Phase 5 conclusion)

1. Quantitative use is limited to mask_common, 0.231-0.268.
2. The upper edge 0.268 is a post hoc Phase 2 validation limit: the 14.8% band is domain-limited, not a demonstrated physical band edge.
3. Raw PhC channels have large frequency-independent phase-intercept errors (0, -1.2035, -2.4075, +2.6725 rad on the chosen branch) and steer to the wrong angle (about 6.3 deg low at f0).
4. The trimmed family is *phase-intercept-trimmed PhC*, one constant correction per channel at f0; never ideal TTD. It requires the trimming to be physically realizable.
5. Trimming changes the phase intercept only, not amplitude, group delay or dispersion.
6. The delay quantum is coarse (0.981 carrier cycle per period), forcing d about 2 lambda0 and grating lobes.
7. The model is 2-D, lossless, uniform, without fabrication disorder, realistic input/output coupling, thermal effects, bends or length-dependent loss; isotropic elements, no mutual coupling.
8. The favorable agreement reflects one flat-n_g D0 mode used at four lengths, not four independently engineered delay elements.
9. Phase 4 used a post-data numerical-method correction (common minimum run time 3500 after the source).
10. The predeclared pointwise central-difference increment estimator was noise-limited; the accepted interpretation uses local-quadratic delay estimates, phase slopes and phase error.
Additionally, for Phase 5: results come from res-30 production (res-40 confirmation agrees); no claim of fabrication readiness, practical antenna readiness, loss performance, 3-D validity, experimental demonstration or grating-lobe suppression is made.

## 8. Unresolved issues

* The tolerance tol_pt (0.05 HPBW) and the pass thresholds are judgments declared before the PhC beam results, but they are not derived from an application requirement.
* The sparse array's grating lobes make the single-aperture result of limited practical value; the real question (subarray TTD, sub-period delays) is untested.
* Whether the trim can be realized in a physical device (a constant phase of up to 2.7 rad per channel) is not addressed.
* Band edge, loss and disorder effects are outside this phase.
* The 0.95 deg residual frequency dependence of the raw lobe is due to the intercept error acting as a fixed phase; it is not a failure of the delay slope.
