# Phases 2–5 (compressed sprint): scope and pre-declared definitions

Written **before any Phase 2–5 FDTD result was inspected**. Anything changed later is recorded in the "Revisions" section at the end,
dated, with the reason. Language: American English (decision D4). Units: a = 1, c = 1, f = a/lambda, delays in a/c.

## 0. Scope (as authorized)

Minimum-complete chain: verified machinery -> W1 guided-mode delay (Phase 2) -> compact delay/bandwidth/dispersion study of three
H-family geometries (Phase 3) -> four-channel delay family (Phase 4) -> four-element array with beam-squint comparison (Phase 5).
Not done in this sprint: exhaustive geometry optimization, eight channels, Monte Carlo/tolerance study, fabrication disorder, 3-D effects.

## 1. Phase 2 — guided-mode delay of the plain r/a = 0.30 W1

* **Device.** Plain W1 (rows |j| = 1 unmodified), eps = 12, TE, supercell rows +-6 (MPB M = 6, FDTD J = 6 rows plus 2 absorbing rows).
* **Guided mode.** Hz-even class; MPB band 8 of that class (MPB `run_te_yodd`; identity cross-checked against the Phase 0.5 parity machinery in
  `2a_mpb.py`). The band decreases in frequency over k in [0, 0.5], i.e. the forward wave has Bloch wavenumber k_fwd = -kappa (v_g > 0).
* **Useful interval (rule fixed now).** Contiguous frequency interval around the minimum of n_g(f) in which n_g,MPB <= 1.10 n_g,min.
  This is chosen *by MPB alone*, before FDTD, to avoid the n_g -> infinity band edges and the n_g rise toward the upper band end.
* **Guided-mode extraction.** FDTD plane fields (Ey, Hz) are projected onto the MPB forward mode and its time-reversed partner with the
  Lorentz-reciprocity overlap (`pcttd.guided.project`). S21 = a_f(out)/a_f(in), S11 = a_b(in)/a_f(in). Mode purity is monitored by the
  relative residual of the projection. (MEEP's own `get_eigenmode_coefficients` was tried first and abandoned: for this periodic structure its
  cross-section mode solve did not reproduce the MPB supercell bands; see `pcttd/w1fdtd.py`.)
* **Lengths.** N = 6, 12, 24 periods at resolution 30 pixels/a (three lengths); N = 6, 12 at resolution 20 and 40 for convergence.
* **De-embedding.** Length differencing: (a) phase of S21(N2)/S21(N1); (b) exact wave-cascade (T-matrix) method
  T(N2) T(N1)^-1 = T_R P(dN) T_R^-1 whose eigenvalues exp(+-i k dN) are independent of both interfaces and of all multiple reflections.
  S22 follows from unitarity of the lossless reciprocal 2-port; unitarity itself (|S11|^2 + |S21|^2 = 1) is a reported diagnostic.
* **Delay estimators.** Central difference of the unwrapped phase and a local quadratic (Savitzky-Golay type, +-3 samples). Frequency
  sampling dependence is quantified by subsampling one fine-grid run (df = 5e-4) to 1e-3, 2e-3 and 4e-3.
* **Independent checks.** (1) delay-vs-length slope (3 lengths); (2) time-domain pulse arrival (centroid of |Hz|^2 envelope between two points of the lead).
* **Gate (pre-declared).** PASS if: guided mode cleanly identified (mode-projection residual small, |S21|^2 + |S11|^2 = 1 within 2%);
  S21 phase-derived delay stable against frequency spacing (change <= 1% from 5e-4 to 2e-3 over the useful interval);
  delay linear in N; and FDTD n_g within **5%** of MPB over the useful interval at the finest resolution with the error decreasing under
  refinement. The aspirational **3%** is reported separately and the continuous error is always reported. A larger but explained error does
  not falsify the pass; an unexplained one stops the sprint.

## 2. Phase 3 — compact delay / bandwidth / dispersion study

* **Geometries** (all points of the existing Phase 0.5 H_r030 sweep, radius factor f1 and outward shift s of the first rows; nothing invented):
  D0 plain (f1 = 1.0, s = 0.0); D1 (f1 = 1.2, s = -0.05); D2 (f1 = 1.2, s = -0.10). Selected from the MPB screening of eight sweep points
  (`results/phase3/3a_mpb_designs.json`): same guided band (Hz-even band 8), band edge moved up to 0.2296 (D1) and 0.2369 (D2), giving
  n_g that rises toward the edge more slowly (D1) and faster (D2). Designs with n_g > ~15 are not used: they are within ~1% of their band edge,
  where the sweep itself found the largest finite-difference sensitivity (Phase 0.5, T6), and the gap-edge/bulk-band coupling is not controlled.
  The Phase 0.5 sweep points with n_g >= 50 are excluded as numerically suspect for FDTD at the available resolution.
* **Structure.** Device section of N periods of the engineered rows inside plain-W1 leads (reference planes in the leads); interface mismatch is real.
  Lengths N = 4, 8, 12.
* **Usable band (rule fixed now).** The widest contiguous frequency interval, inside the design's guided band, over which ALL hold:
  (i) device transmission T = |S21|^2 >= 0.5 for every length used; (ii) the per-period delay tau(f) = d phi_section / d omega / dN stays within +-10% of
  its median over the interval (this defines the delay ripple); and, for the system study, the interval must also contain the array center frequency.
  Sensitivity: +-5% and +-25% are also reported. Fractional bandwidth = width / center.
* **Reported per design:** complex S21, S11, T, R, per-period group delay (and n_g = c tau), delay ripple within the usable band (max-min)/median,
  GDD per period (d tau / d omega) at the band center and its band maximum, usable fractional bandwidth, group-delay-bandwidth product
  (fractional bandwidth x n_g), and length needed for a stated delay target.
* **Pareto language.** Three points are a tradeoff illustration, not a Pareto frontier.
* **Conventional reference.** The R air-channel waveguide (square lattice, TM) is compared by MPB only tonight (delay per length, n_g flatness,
  fractional bandwidth, GDD); no FDTD for R is run. This is disclosed as an **incomplete reference comparison**.

## 3. Phase 4 — four-channel proof of concept

* Four channels of one validated structure with lengths N0 + m dN, m = 0..3, simulated **directly** (no cascade needed to build the channels);
  the cascade/T-matrix machinery is validated separately by predicting the longest Phase 2 length from the two shorter ones.
* Saved outputs are the actual complex H_m(f) (S21 between the common reference planes), never only delays.
* Reported: amplitude mismatch, incremental delay error (vs ideal m * dtau), phase error, usable common band (intersection of per-channel rules above).
* **Raw PhC** vs **phase-intercept-trimmed PhC**: trimmed applies one frequency-independent phase per channel,
  delta_m = arg H_m(f0) - arg H_0(f0) - 2 pi f0 m dtau_target, chosen at the center frequency f0 only. The two are always labeled.

## 4. Phase 5 — four-element array (pure Python)

* Uniform linear array, four isotropic elements, pitch d, steering angle theta_0 = 30 degrees, center frequency f0 (center of the common usable band).
* **Delay granularity (important, pre-declared).** One lattice period delays the guided wave by n_g/c ~ 3.9 a/c, i.e. about **one carrier cycle** (f0 n_g ~ 1).
  The ideal inter-element delay of a half-wavelength array at 30 degrees is d sin(theta_0)/c = 0.25 cycle, smaller than the smallest increment the lattice allows.
  Therefore d is **not** set to lambda_0/2: it is set so that the target step equals the simulated per-period delay step,
  d = c * dtau_step(f0) / sin(theta_0) (about 2 lambda_0 for one period per step). The array is therefore sparse (grating lobes exist);
  the design lobe is the m = 0 lobe, and the pointing angle is the maximum of |AF| inside the window |sin(theta) - sin(theta_0)| <= c/(2 f_max d)
  (half the grating-lobe spacing at the highest frequency of the band). The squint physics of the m = 0 lobe (sin theta ~ (f0/f) sin theta_0 for phase-only steering)
  is independent of d. This limitation is reported as a finding (sub-period delay increments, subarray-level TTD, or an RF-photonic mapping are the ways out).
* Cases: (1) fixed phase shifters (ideal, frequency independent) set for theta_0 at f0; (2) ideal TTD, H_m = exp(i 2 pi f m dtau_target) (stated sign convention);
  (3) raw PhC H_m(f) from Phase 4; (4) trimmed PhC.
* Metrics: theta_peak(f); R_squint = max_f |theta_PhC - theta_0| / max_f |theta_phase - theta_0| over the stated common usable band (actual value reported;
  <= 0.20 is a stretch target only); peak-gain change, half-power beamwidth, highest sidelobe, channel amplitude imbalance; patterns at lower edge, center, upper edge.

## 5. Revisions

**Revision 1 (2026-10-03, recovery after an interrupted session; before any durable Phase 2 result was committed).**
1. *Scope superseded.* The "compressed sprint" scope (Section 0: four channels, four elements, no eight-channel/robustness study, no Monte Carlo) is
   superseded by the full roadmap (Phases 2-7). The Phase 2 and Phase 3 *definitions* above (useful-interval rule, gate thresholds, usable-band rule)
   are retained unchanged; Phase 4-5 definitions are retained as the four-channel stage, to be extended to eight channels. The full-roadmap
   Phase 2 gate additionally requires a three-length slope/intercept fit, an MPB derivative-grid convergence study, FDTD-vs-MPB errors reported
   absolutely and relatively, and collected tests; these are added requirements, not relaxations.
2. *Disclosure on "written before any result was inspected".* The interrupted session had already run exploratory short-W1 runs (S21 magnitudes,
   S11 of a few percent, projection residuals 0.03-0.17, decay-tolerance phase stability) *before* this document was recovered; those numbers were
   used only to debug the projection method and set run lengths. No gate threshold above was changed after seeing them. Those exploratory values are not
   evidence; the durable Phase 2 results are re-derived from committed code.
3. *Provenance of recovered files.* After the crash, `tests/test_guided.py`, `tests/test_array.py`, `tests/test_ttdmetrics.py` and `src/pcttd/ttdmetrics.py` were found
   as zero-byte files (their content was lost), which is why pytest still reported 98 tests. Tests for `guided`, `array` and `w1modes` were rewritten
   from scratch against the recovered library code; `ttdmetrics` is to be re-implemented. `results/phase2/2a_mpb_plain.npz` was found without its
   json summary (the script had not finished) and is regenerated by re-running `experiments/phase2/2a_mpb.py`.

**Revision 2 (2026-10-04, after the Phase 2 data were inspected).** Recorded so that nothing below can be mistaken for a pre-declared rule.
1. *Decay tolerance.* The runs use `stop_when_dft_decayed`. At tol = 1e-6 (the value used in the first production runs) the stop time was not monotone
   in device length (res 40: t_end = 883, 531, 807 for N = 6, 12, 24) and the pointwise n_g noise (N = 6/12 pair) was 6.7-8.9% rms; tol = 1e-7 gave 1.6-2.0% and
   tol = 1e-8 gave 1.3% (3-length fit: 0.8-0.9% for both). The production configuration is therefore tol = 1e-7. This is a numerical-method change motivated by
   the convergence study, not a gate relaxation.
2. *Validated sub-interval.* The pre-declared useful interval (n_g <= 1.10 n_g,min) is [0.2312, 0.2742]. The unitarity criterion (|S11|^2 + |S21|^2 = 1 within 2%) is **not** met
   at N = 24 above f ~ 0.268 (defect up to 4.4% at 0.2742, growing with length, accompanied by growth of the projection residual and of the n_g bias). The
   boundary 0.268 was chosen **after** seeing these data. All quantitative Phase 3-5 results that rely on FDTD-validated guided transmission are restricted to
   [0.231, 0.268] (fractional bandwidth 14.8%); the upper 0.268-0.2742 part is reported but flagged. Cause (partly verified): a transverse-cell test (J = 8 rows instead of 6, with a
   matching M = 8 MPB profile) reduces the defect at 0.270-0.2742 by about 4x (N = 24: +0.0202 -> +0.0054) and the MPB n_g itself shifts by +0.9% between M = 6 and M = 8 there, i.e. the mode's
   lateral tail is long enough to feel the finite lateral cell near the band top. This explains the *direction* of the effect; the residual n_g bias (+2.2% against the matching M = 8 reference) is not fully explained.
3. *Time-domain pulse check* uses two carrier frequencies (0.235, 0.2564), not three.
4. *Phase 2 gate wording.* The n_g criterion (5%, aspirational 3%) is applied to the *3-length phase slope* estimator (the others are reported). Max error over the full useful interval at the
   finest resolution: 2.9% (res 40), 3.2% (res 30); inside the validated sub-interval <= 1.4%.

**Revision 3 (2026-10-09, Phase 3 start; written after the Phase 3 MPB screening 3a was regenerated but BEFORE any Phase 3 FDTD result was inspected).** This fills a gap; it
relaxes or tightens no existing criterion.
1. *Target delay for the "required length" metric (not previously fixed).* A fixed, design-independent normalization: **25 carrier cycles** at each band's own
   center frequency, tau_target = 25 / f_center (in a/c); equivalently the same physical delay for every design once each lattice is scaled to the same operating
   frequency. Reported as periods L = tau_target / tau_pp, and in free-space wavelengths at the band center, L f_center. A second target of 100 a/c is not used.
2. *Validated-domain restriction.* The Phase 3 FDTD usable band is searched inside the Phase-2-derived validated range [0.231, 0.268] (a Phase 2 pipeline
   validation limit, NOT a Phase 3 predeclared criterion; its upper edge was itself chosen post hoc, Revision 2). The same search without that restriction (analysis range
   0.22-0.275) is also reported and labeled "not validated". The primary usable band for gate purposes is the restricted one.
3. *Delay estimator.* tau_pp(f) = (1/dN) d phi_section / d omega from the three-length linear fit of the S21 phase (branch fixed by MPB of the device section),
   cross-checked by the wave-cascade phases of the pairs (4,8), (8,12), (4,12). Frequency grid 5e-4. The +-10% rule is primary; +-5% and +-25% are sensitivity only.
4. *Production numerics.* res 30, decay_tol = 1e-7, source/monitor distances as in Phase 2 (d_src = 10, d_mon = 3); res 40 confirmation for selected cases.
5. *R reference.* MPB only; usable band by the same +-10% n_g rule inside its single-mode TM air-channel band; incomplete reference comparison.

**Revision 4 (2026-10-10, Phase 4 start; written BEFORE any Phase 4 FDTD result exists or was inspected).** Fills gaps in Section 3; relaxes or tightens nothing earlier.
1. *Channel family.* One structure, D0 (plain W1, r/a = 0.30, f1 = 1.0, s = 0.0). Lengths N_m = N0 + m dN with **dN = 1 period** and **N0 = 8** (N = 8, 9, 10, 11).
   Rationale: an integer-period step is the smallest physically meaningful increment (about one carrier cycle, see Section 4); N0 = 8 lies inside the length range already
   validated in Phases 2-3 (N = 4-24), keeps all four devices at least two source/monitor-free periods long, and adds no propagation length beyond what the
   cascade check needs (it predicts N = 11 from N = 8, 9 and N = 9, 10). Chosen without reference to any Phase 4 output. The coarse one-period delay quantum is
   characterized, not hidden; no sub-period interpolation is done in Phase 4.
2. *Target delay step.* dtau_target = dN * 3.9323 a/c, where 3.9323 is the Phase 3 D0 median per-period delay over the validated +-10% usable band
   (`results/phase3/3e_D0_res30.json`, key `tau_fit|validated|tol10`, `median_tau`; res 40 gives 3.9308). Fixed now; never refit to Phase 4 data.
3. *Center frequency for trimming.* f0 = 0.2495 (center of the Phase 3 D0 validated band [0.231, 0.268]; a point of the 5e-4 frequency grid). The Phase 4 common band center is
   reported separately; if it differs from f0 by more than 0.001, the handoff additionally carries a secondary trim at that center, labeled as such.
4. *Reference planes and what H_m contains.* Input plane x = -3, output plane x = N_m + 3 (d_mon = 3), source d_src = 10, res 30 (production) and res 40 (confirmation), decay_tol = 1e-7,
   frequency grid 0.220-0.275 step 5e-4, lead profiles from the Phase 2 cache. H_m(f) = S21 = a_f(out)/a_f(in) of the full finite device, hence it contains the *entire* propagation phase over N_m + 6 periods;
   nothing is de-embedded, so the reference-plane offset is common to all channels and cancels in H_m/H_0 (the lead is plain W1, identical to the D0 section: no interface).
   Convention (unchanged): MEEP exp(-i omega t), MEEP_SIGN = +1, a delay tau gives exp(+i omega tau), tau_g = +d arg H / d omega. Ideal channel m: H_m/H_0 = exp(+i 2 pi f m dtau_target).
   Branch of the unwrapped phase fixed by the MPB Bloch phase, -2 pi kappa(f) x length (as in Phase 3).
5. *Common usable band.* Widest contiguous frequency interval inside [0.231, 0.268] on which ALL hold for all four channels: (i) T_m >= 0.5; (ii) |unitarity defect| <= 0.02 (Phase 2 criterion);
   (iii) projection residual (in and out) <= 0.10; (iv) every adjacent increment dtau_{m+1,m}(f) (central difference of the unwrapped phase of S21_{m+1}/S21_m, grid 5e-4) lies within +-10% of dtau_target
   (+-5%, +-25% and the local-quadratic estimator are sensitivities only); (v) amplitude imbalance max_m|H_m| / min_m|H_m| <= 1 dB (20 log10, amplitude, not power).
6. *Errors.* epsilon_tau_m(f) = Dtau_m(f) - m dtau_target with Dtau_m = d arg(H_m/H_0)/d omega, also normalized by m dtau_target. Raw phase error
   epsilon_phi_m(f) = arg(H_m/H_0)(f) - 2 pi f m dtau_target, unwrapped along f and placed on the 2 pi branch with epsilon_phi_m(f0) in (-pi, pi]. Trim: delta_m = epsilon_phi_m(f0) (identical to the criteria formula of
   Section 3), H_m^trim = H_m exp(-i delta_m); epsilon_phi,trim_m = epsilon_phi_m - delta_m, exactly 0 at f0. Reporting threshold for "useful broadband phase matching": max |epsilon_phi,trim| over the common band <= pi/8 (0.39 rad); it is not
   a band criterion and does not change the band.
7. *Cascade validation.* T-matrix of directly simulated shorter devices predicts the directly simulated N = 11: from (8, 9) with three steps and from (9, 10) with two steps,
   T(11) = [T(Nb) T(Na)^-1]^(11-Nb) T(Nb). Acceptance (over the common band): max |phase difference| <= 0.1 rad, max relative |S21| difference <= 3%.
8. *Resolution check (res 40 vs res 30), all four channels.* "Conclusions unchanged" means: band-mean adjacent increment changes by <= 2% (of target), common-band edges move by <= 0.002,
   max amplitude imbalance changes by <= 0.1 dB, max trimmed phase error changes by <= 0.1 rad.
9. *Phase 4 gate questions* are those listed in the Phase 4 task statement (direct simulation and numerical trust, monotonic delay, increments near target, nontrivial common band, amplitude mismatch, raw phase
   error characterized, trimming effect, cascade agreement, trustworthiness of H_m for Phase 5). PASS WITH RESTRICTIONS if a restriction (band, delay-quantum, noise) must travel with the handoff; FAIL if monotonicity, trust or a nontrivial common band is lost.

**Revision 5 (2026-10-10, written AFTER the first res-30 Phase 4 data were inspected and BEFORE any res-40 Phase 4 data were inspected).** Everything below is post hoc and is labeled as such in the report. None of it loosens the
amplitude, unitarity, transmission, residual or domain criteria, or the trim definition.
1. *Stop time (numerical method, not a criterion).* The first attempt (`stop_when_dft_decayed`, tol 1e-7, minimum run time max(150, transit)) stopped at erratic times that were not monotone in length
   (res 30: t_end = 1965, 801, 2816, 2462 for N = 8, 9, 10, 11; res 40: 636-1043 for N = 9-11), which is the Phase 2 stop-time anomaly. N = 9 was visibly under-converged (unitarity defect 1.9%, T up to 1.018, a
   spurious increments pattern, and a cascade phase mismatch > 1 rad). The attempt was discarded (kept for transparency in `results/phase4/raw_attempt1/` and `attempt1_*`) and **all** channels were re-run with a common minimum run time
   `min_after = 3500` after the source (decay tolerance 1e-7 still required), giving t_end = 3524-4586 at res 30 and a defect <= 0.34%. Raw phase errors of the first attempt agreed with the final ones to ~0.01 rad (the phase
   itself was robust; the amplitude, unitarity and the derivative were not). The same `min_after` is used at res 40.
2. *Adjacent-increment estimator.* The pre-declared pointwise central difference of the single-period phase ratio (criterion 5(iv)) is noise-limited: the rms difference between it and the local-quadratic
   (+-3 samples) estimator of the same data is 5.5-7.5% of dtau_target (final res-30 data), comparable to the +-10% tolerance, so it fragments the band into isolated samples (widest run 5.0%, 0.2425-0.2550) even though the underlying
   phase is smooth. This is the Phase 2/3 finding (pointwise derivatives of FDTD phase are noise-limited), now met at one-period spacing. The **local-quadratic estimator is therefore the primary one for criterion 5(iv)**;
   the pre-declared central-difference band is reported next to it as the conservative, noise-limited result, and the band-fit (least-squares phase slope over the band) increments and phase errors, which are not derivative-noise limited,
   are the quantitative basis for the delay statements. The tolerance (+-10% of dtau_target), the target and the other four conditions are unchanged.
3. *Erratum, cascade step counts.* Criteria Revision 4 item 7 said "three steps" for (8, 9) and "two steps" for (9, 10). The formula given there is the governing one: T(11) = [T(Nb) T(Na)^-1]^((11-Nb)/(Nb-Na)) T(Nb), i.e. two
   steps from (8, 9) and one step from (9, 10). An additional check predicts N = 10 from (8, 9) (one step). Acceptance thresholds are unchanged (<= 0.1 rad, <= 3%).

**Revision 6 (2026-10-10, Phase 5 four-element array; written BEFORE any Phase 5 array result from the Phase 4 transfer functions was computed, and before any res-40 beam result was inspected).**
1. *Array.* Four isotropic elements, x_m = m d, theta0 = 30 deg, c = 1, f0 = 0.2495, dtau_target = 3.9323 a/c (handoff values). Pitch **d = dtau_target / sin(theta0) = 7.8646 a (about 1.96 lambda0)**; d is NOT set to lambda0/2 because the one-period delay quantum (0.981 carrier cycle) would then
   correspond to 4 times the half-wave step. The array is intentionally sparse; grating lobes are expected and reported. Convention: `pcttd.array` (AF = sum_m H_m exp(-i k0 x_m sin theta); ideal TTD H_m = exp(+i 2 pi f m dtau_target); no code change to the convention).
2. *Cases.* (1) fixed phase, H_m = exp(+i 2 pi f0 m dtau_target) at all f; (2) ideal TTD; (3) raw PhC = handoff `H_raw` (amplitude and phase, as is); (4) phase-intercept-trimmed PhC = handoff `H_trim`. Case 3 and 4 are never replaced by fitted delays. Normalized-amplitude variants (|H_m| set to 1, phases kept) are diagnostics only.
3. *Domain.* All primary numbers use `mask_common` (0.231-0.268; domain-limited, post hoc upper edge). Anything else is diagnostic and marked.
4. *Design lobe.* Tracked inside |sin theta - sin theta0| <= 1/(2 f_max d), f_max = highest frequency of mask_common. Peak found on a 40001-point grid in sin theta inside this window, refined by a parabola through the three highest power samples (documented and unit-tested; angular grid error must be < 0.001 deg).
   Full pattern: |AF|^2/16 on -90..+90 deg; lobes are the local maxima of the power on a 20001-point grid in sin theta (endpoint maxima at +-90 deg included and flagged); "strongest sidelobe/grating lobe" = the strongest local maximum other than the design-lobe peak, in dB below the design-lobe peak. Grating-lobe count = local maxima within 3 dB of the design-lobe peak. HPBW is local: -3 dB crossings around the design-lobe peak, linearly interpolated, in theta.
5. *Metrics.* Pointing error theta_peak - theta0 (max, RMS, low edge, f0, high edge); R_squint = max|theta_PhC - theta0| / max|theta_fixed - theta0| (as defined; includes any absolute offset), plus a diagnostic offset-removed R_squint,spread = max|theta_PhC - theta_PhC(f0)| / max|theta_fixed - theta0|; ideal-TTD error theta_PhC - theta_ideal; peak power gain |AF|^2/16 (also in dB); amplitude diagnostic; phase-error link. R_squint <= 0.20 is a stretch target only and not a criterion.
6. *Analytic checks (must pass before PhC results are interpreted).* Fixed phase: sin theta_peak = (f0/f) sin theta0 to <= 0.01 deg over the band. Ideal TTD: |theta_peak - theta0| <= 0.001 deg over the band, peak power gain 1.
7. *Tolerance tol_pt = 0.05 x HPBW_ideal(f0)* (about 0.3 deg; fixed before looking at beam results from the PhC channels). Used as the single yardstick for: "trimmed PhC indistinguishable from ideal" (max |theta_trim - theta_ideal| <= tol_pt over the band); "raw PhC has an absolute pointing offset" (|theta_raw(f0) - theta0| > tol_pt); "amplitude mismatch matters" (change of theta_peak > tol_pt anywhere, or of peak gain > 0.1 dB).
   "PhC reduces squint" means the offset-removed spread max|theta_PhC - theta_PhC(f0)| is < 0.5 x that of fixed phase.
8. *Phase-error link.* Per-frequency least-squares slope g(f) of the residual phase error epsilon_m(f) (vs element index m) predicts a pointing shift Delta theta ~ g / (2 pi f d cos theta0); the prediction is compared with the numerical theta_PhC - theta_ideal.
9. *Resolution check (res 40 vs res 30 H arrays, same analysis).* "Beamforming conclusions unchanged" means: max over the band |theta_peak(res40) - theta_peak(res30)| <= tol_pt/3 (about 0.1 deg) for every case; change of max |pointing error| <= tol_pt/3; change of R_squint (raw and trimmed) <= 0.02 absolute; change of peak power gain <= 0.02 (normalized power); change of max |theta_trim - theta_ideal| <= tol_pt/3; change of strongest non-design lobe level <= 0.5 dB.
10. *Phase 5 gate questions.* (1) Does fixed-phase steering reproduce the analytic squint? (2) Does ideal TTD stay at theta0? (3) Does raw PhC reduce frequency-dependent squint vs fixed phase? (4) Does raw PhC have an absolute pointing offset from its Phase 4 intercept error? (5) Does the single constant trim place the lobe near
    theta0 without changing broadband delay behavior? (6) How close is trimmed PhC to ideal TTD? (7) Does amplitude mismatch matter for pointing or gain? (8) Are conclusions stable between res 30 and res 40? (9) How severe are the grating lobes from d ~ 2 lambda0? (10) Is the concept scientifically promising for squint suppression despite the coarse delay quantum?
    PASS WITH RESTRICTIONS if the channels clearly demonstrate squint suppression but use is constrained (delay quantum, grating lobes, phase trim, lossless model, validated band); FAIL only if the actual H_m give untrustworthy or physically inconsistent array behavior (analytic checks fail, or res 30/40 disagree beyond item 9).
11. *Restrictions carried from Phase 4* (all ten, see handoff and `docs/phase4_report.md`) appear in every Phase 5 conclusion. The trimmed family is always labeled "phase-intercept-trimmed PhC".

**Revision 7 (2026-10-10, Phase 6: delay-granularity mitigation and robustness; written BEFORE any Phase 6 array, subarray, Monte Carlo or res-40 result was computed. Only analytic reasoning preceded it).** Nothing earlier is changed.
1. *Common definitions.* c = 1, theta0 = 30 deg, f0 = 0.2495 (lambda0 = 4.008 a), tau_PhC = dtau_target = 3.9323 a (0.9811 cycle at f0), mask_common (0.231-0.268), conventions of `pcttd.array` (AF = sum_n W_n exp(-i k0 x_n sin theta); delay tau -> exp(+i 2 pi f tau)).
   Element pitch d_e = lambda0/2 = 2.004 a. Subarray q = 0..3 (the four Phase 4 channels, N = 8-11), Msub elements per subarray, N_el = 4 Msub, element (q, p) at x = (q Msub + p) d_e. Element weight W_(q,p)(f) = G_q(f) exp(+i 2 pi f0 p d_e sin(theta0)), G_q = the subarray channel
   transfer function. The full element-level array factor is always computed (no point-source subarrays). Power is |AF|^2 / N_el^2.
2. *Architectures.* A = Phase 5 sparse element-level PhC (d = 1.96 lambda0, 4 elements), read unchanged from `results/phase5`. B = ideal element-level TTD at d = lambda0/2 (4 elements; required step lambda0/4 = 1.002 a = 0.25 cycle; ratio tau_PhC / tau_half reported exactly), ideal reference only.
   C = hybrid element-level TTD + constant phase (Msub = 1, 4 elements at lambda0/2, one PhC channel per element, constant phase trim at f0). D = subarray-level TTD, Msub = 2, 4, 8 (8, 16, 32 elements). Required step between subarray centers tau_req(Msub) = Msub lambda0/4 = Msub x 1.002 a
   (Msub/4 cycle at f0). Reported per Msub: ratio rho = tau_PhC / tau_req, mismatch rho - 1, and f0 (tau_PhC - tau_req) in cycles.
3. *Cases for each Msub (1, 2, 4, 8).* (i) fixed phase (all N_el elements frequency-independent phase set at f0); (ii) ideal element-level TTD (every element exact TTD); (iii) ideal subarray TTD (exact tau_req between subarrays, fixed phase inside); (iv) raw PhC subarray (handoff H_raw as is, amplitude and phase);
   (v) phase-intercept-trimmed PhC, **re-trimmed to the required step**: one constant phase per channel chosen at f0 so that arg G_q(f0) - arg G_0(f0) = 2 pi f0 q tau_req(Msub) (`array.trim_phase`); the delay slope and amplitude are untouched. This is the "phase-corrected delay-mismatch case": it fixes the carrier-phase intercept only, never the delay increment;
   (vi) Phase 4 `H_trim` as is (trimmed for 3.9323 a, not for tau_req); only physically consistent for Msub = 4 and reported for the others as a diagnostic. Trimmed PhC is never called ideal TTD. Normalized-amplitude diagnostics (|H| = 1) as in Phase 5.
   *Analytic expectation stated beforehand (to be checked, not assumed):* a constant phase cannot change group delay. For a retrimmed PhC the subarray-level steering obeys sin(theta) - sin(theta0) = (f0/f - 1) sin(theta0) (1 - rho), i.e. the squint of the subarray-level factor is (1 - rho) times the fixed-phase squint; rho = 3.92 (Msub = 1), 1.96 (2), 0.981 (4), 0.49 (8). The Msub = 4 hypothesis (rho ~ 1) is accepted only if numerics confirm.
4. *Design lobe.* Tracked in |sin theta - sin theta0| <= 1/(2 f_max Msub d_e) (half the subarray-level grating-lobe spacing at f_max, clipped to [-1, 1]; for Msub = 1 this is the Phase 5 rule), 40001-point grid, parabolic refinement. Window saturation (|s_peak - s0| >= 0.98 window) is flagged.
   Full pattern -90..+90 deg on a 20001-point sin(theta) grid (`beam.lobe_report`); strongest other lobe = strongest local maximum other than the design lobe, dB below the design peak. A lobe within -10 dB of the design peak is called an *ambiguity*. For reference a uniform equal-amplitude array has its first sidelobe at about -13.3 dB.
5. *Metrics (over mask_common).* theta_peak(f); max and RMS |theta - theta0|; error at low/f0/high edge; peak-to-peak spread; R_squint = max|theta - theta0| / max|theta_fixed - theta0| (and the offset-removed spread version); difference to ideal element-level TTD; peak power gain (normalized to N_el^2); HPBW; strongest other lobe; number of lobes within 3 dB; ambiguity flag; amplitude sensitivity from the |H| = 1 diagnostic; aperture N_el d_e / lambda0.
   Yardstick tol_pt(Msub) = 0.05 x HPBW of ideal element-level TTD of that architecture at f0 (fixed from the ideal case, never from PhC data).
6. *Acceptable beam tracking.* max_f |theta_peak - theta0| <= tol_pt(Msub). *Squint suppressed* if the offset-removed spread of the PhC case is < 0.5 of the fixed-phase max error (as Revision 6); R_squint <= 0.20 is a stretch target only. *Grating-lobe problem substantially mitigated* if the strongest other lobe is <= -10 dB at every frequency of mask_common
   (Phase 5: about 0 dB) and no ambiguity exists. *Msub = 4 alignment accepted* if |rho - 1| <= 5% and the retrimmed Msub = 4 case meets the tracking and squint criteria. The pointing error caused by the 1.9% delay mismatch is read separately from (a) the retrimmed case (slope error only) and (b) Phase 4 `H_trim` as is (slope plus intercept error).
   The best architecture is the one among the evaluated PhC cases that meets tracking, squint and lobe criteria with the largest margin on max pointing error; if none qualifies, the architecture with the smallest max error and strongest-lobe level is reported as best and the gate cannot be PASS.
7. *Robustness (Part B), best architecture only, system-level, nominal = complex Phase 4 G_q(f) on mask_common.* Perturbation of channel q: G_q -> G_q a_q exp(i eps_q) exp(+i 2 pi f dtau_q). Independent uniform draws per channel, **all four channels**, half-widths: phase eps_q ~ U(-a, a), a in {0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30} deg
   (headline: 1, 2, 5, 10); delay dtau_q = e_q q tau_PhC with e_q ~ U(-e, e), e in {0.1, 0.25, 0.5, 1, 2, 3, 5, 10} % of the channel's *relative* delay q tau_PhC (headline 0.25, 0.5, 1, 2, 5; the same values also expressed in a/c, periods and cycles;
   an alternative scaling to the absolute device delay (N_q + 6) tau_pp is reported as sensitivity); two delay variants: *uncalibrated* exp(+i 2 pi f dtau) (the constant phase at f0 changes too; the headline) and *calibrated* exp(+i 2 pi (f - f0) dtau) (intercept re-trimmed, pure slope error);
   amplitude a_q = 10^(u_q/20), u_q ~ U(-A, A), A in {0.1, 0.25, 0.5, 1, 2, 3, 5, 6} dB (headline 0.1, 0.25, 0.5, 1). Combined (simultaneous, uncalibrated delay): levels L1 = (1 deg, 0.25 %, 0.1 dB), L2 = (2 deg, 0.5 %, 0.25 dB), L3 = (5 deg, 1 %, 0.5 dB), L4 = (10 deg, 2 %, 1 dB).
   Perturbations are applied to the retrimmed (Msub = 4) channels, with the trim itself not redone, i.e. errors are residual errors after calibration.
8. *Monte Carlo.* Seed 20261010 (numpy `default_rng`, one child stream per tolerance point via `SeedSequence.spawn`), **1000 samples per tolerance point** (same samples reused across tolerance values only through the spawn structure; no resampling after seeing results). Reported percentiles: median, 90th, 95th, worst observed of: max pointing error, RMS pointing error, peak gain loss
   (minimum over the band of the peak power, dB, relative to the nominal retrimmed case at the same frequency), strongest other lobe (maximum over the band, dB), and R_squint. The fast Monte Carlo engine (exact factorization of the array factor into the element-in-subarray factor times the four-term subarray sum, grid ds = 2.5e-4) is validated against `beam.analyze_case` on the nominal case before use (agreement <= 0.01 deg, <= 0.01 in power).
9. *Tolerance budget.* "System-level tolerance estimate" = the largest independent-parameter half-width a, linearly interpolated on the sweep, for which the **95th percentile of the max pointing error over mask_common (including the nominal error)** stays <= 0.5 deg, and separately <= 1.0 deg; also stated for the median. Not a fabrication specification. If the nominal error alone exceeds a limit, it is reported as unattainable.
10. *Analytic error propagation.* Phase error eps_q: shift of sin(theta) = g c / (2 pi f Msub d_e) with g the least-squares slope of eps_q over q (`beam.phase_slope_pointing` with d -> Msub d_e), Delta theta = Delta(sin theta) / cos(theta0). Delay error dtau_q: eps_q(f) = 2 pi f dtau_q, so the prediction is frequency-independent for uncalibrated delay errors
    (Delta(sin theta) = slope_q(dtau) / (Msub d_e)) and proportional to (f - f0)/f for calibrated ones. Compared with the numerical deterministic response (not just Monte Carlo statistics); agreement is reported as the max difference in degrees.
11. *Period quantization.* The set of realizable delay steps is k tau_pp (integer k periods; tau_pp(f0) = 3.9323). For Msub = 1..12 and k = 1..3 the ratio rho(Msub, k) = k tau_PhC / (Msub lambda0/4) and the implied squint fraction |1 - rho| are tabulated; fractional periods only as a labeled hypothetical. Only the 1-period channel family exists, so
    the numerical array results use k = 1; Msub = 8 with k = 2 is not simulated.
12. *Resolution check.* The best architecture is recomputed with `H_raw_res40`, `H_trim_res40` (retrim computed from the res-40 raw data). "Conclusion unchanged" means: max over the band |theta_peak(res40) - theta_peak(res30)| <= tol_pt / 3; change of max |pointing error| <= tol_pt / 3; change of R_squint <= 0.02 absolute; change of peak power <= 0.02 (normalized); change of strongest-other-lobe level <= 0.5 dB;
    best-architecture identity and gate verdict unchanged.
13. *Phase 6 gate questions* are the thirteen listed in the Phase 6 task. PASS WITH RESTRICTIONS if some coherent architecture uses the demonstrated delay behavior, squint suppression survives (criterion 6), the grating-lobe problem is substantially mitigated (criterion 6), the res 30/40 check holds, and the P95 max pointing error is <= 1.0 deg at combined level L2 (the "not too sensitive" test, declared here).
    FAIL if no architecture can use the delay quantum coherently, ambiguities remain, or the L2 test fails. PASS (no restrictions) is not available: the 2-D, lossless, uniform, no-disorder, domain-limited restrictions of Phases 4-5 travel with every Phase 6 conclusion.
14. *Scope.* No new FDTD, no fabrication disorder simulation, no eight-channel extension, no Phase 7.
