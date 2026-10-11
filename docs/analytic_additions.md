# Analytic additions: steering-angle sweep and conventional propagation references

Written for: the project author deciding how to revise the manuscript. Not for submission.
Status: analytic/post-processing only. No new FDTD, no Phase 4–6 result changed, manuscript not modified. All thresholds below are post hoc interpretive thresholds, not gates.

Scripts: `experiments/phase7/7e_angle_sweep.py`, `experiments/phase7/7f_conventional_reference.py`. Outputs: `results/final/angle_sweep.{json,npz}`, `results/final/conventional_reference.json`, `figures/final/fig10_angle_sweep.{png,pdf}`, `figures/final/fig11_match_map.{png,pdf}`. Tests: `tests/test_analytic_additions.py`.

## 1. Purpose

The closest-prior-art audit left the Msub = 4 / 30° delay-quantum match as the central result. Two questions remained: (i) is 30° an isolated coincidence or part of a usable angular window, and (ii) how much shorter is the W1 delay line than plain propagation for one carrier cycle? Both are answered here from numbers already in the repository.

## 2. Inputs and conventions

* Transfer functions: `results/phase4/phase4_handoff.npz` raw `H_raw` (res 30), common band 0.231–0.268 c/a (75 frequency points, f0 = 0.2495), delay increment 3.9323 a/c per period → q_phc = f0·τ = **0.981109** carrier cycle.
* Model: `pcttd.subarray.Hier` (exact factorized hierarchical array), `retrim`, `fixed_channels`, `ideal_channels`, `intra_fixed`, `intra_ttd` — the Phase 6 code, unchanged. Four subarrays of Msub λ0/2-spaced isotropic elements; sign convention exp(−iωt), delay τ → exp(+iωτ); power |AF|²/N_el²; design-lobe window 1/(2 f_max Msub d_e) in sin θ; 40001-point sin grid with parabolic refinement; strongest-other-lobe and subarray-replica definitions as in Phase 6 (`Hier.lobes`, `Hier.replica_levels`).
* Required step per adjacent subarray: q_req = f0·Msub·d_e·sin θ0 = Msub·sin θ0 / 2 carrier cycles.
* **Re-trimmed PhC (case A):** for each θ0 the Phase 4 raw channel responses are used with their delay slope untouched; only one constant phase per channel at f0 is recomputed so that the f0 phase step equals the required step (identical to the Phase 6 "trim"). No rescaling, no refitted delay.
* **Ideal subarray TTD (case B):** exact TTD between subarrays with the same fixed-phase weights inside the subarrays. It separates delay-quantum mismatch from intra-subarray fixed-phase squint. **Fixed phase (case C)** is the baseline for R_squint.
* The sweep grid is 5°–70° in 1° steps plus 30° and the exact match angle (67 angles).
* Check: at 30° the sweep reproduces the Phase 6 trim, fixed, ideal-subarray and R_squint values to relative 1e-9 (test).

## 3. Steering-angle analytic relation (standard array theory)

With k PhC periods per adjacent subarray increment, a re-trimmed channel family matches the required delay when

 Msub sin θ_match / 2 = k f0 τ_PhC ⇒ **sin θ_match = 2 k f0 τ_PhC / Msub = 1.96222 k / Msub.**

The match depends only on k/Msub. With ρ = k q_phc / q_req, the delay-mismatch contribution to the pointing is the Phase 6 expression sin θ_pk − sin θ0 ≈ (f0/f − 1) sin θ0 (1 − ρ). This is an application of ordinary array theory, not a new law. The match condition is a rewriting of "required delay = available delay"; the relation is not claimed as novel.

## 4. Msub/k match-angle table

| Msub | subarray pitch | k = 1 | k = 2 | k = 3 | k = 4 | k = 5 | k = 6 | k = 7 | k = 8 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.5 λ0 | none (sin = 1.96) | – | – | – | – | – | – | – |
| 2 | 1 λ0 | 78.85° | – | – | – | – | – | – | – |
| **4** | 2 λ0 | **29.38°** | 78.85° | – | – | – | – | – | – |
| 8 | 4 λ0 | 14.20° | 29.38° | 47.38° | 78.85° | – | – | – | – |
| 16 | 8 λ0 | 7.04° | 14.20° | 21.59° | 29.38° | 37.82° | 47.38° | 59.15° | 78.85° |

(– = sin θ > 1, not visible.) Exact sines: Msub = 4: 0.4906 (k = 1), 0.9811 (k = 2). Highlights:

* **Msub = 4, k = 1 matches at 29.3770°, not 30°.** At 30° the 1.89 % mismatch is a small offset from the exact match.
* Msub = 1 can never match (one period delays 0.981 cycle, an λ0/2 element needs at most 0.5).
* Matches at ≥ 78° (Msub = 2, k = 1 and Msub = 4, k = 2) are near end-fire and of no practical interest.
* Larger Msub gives matches at smaller angles in steps of 2 f0 τ / Msub in sin θ, i.e. a ladder of angles, each with a wider subarray pitch. A pitch of Msub/2 λ0 makes Msub − 1 subarray-level replicas visible at f0 (test-checked in the JSON). Fig. 11 shows the map.
* k > 1 is arithmetic only: no multi-period channels were simulated. The assumption that the Bloch-phase intercept remains trimmable by one constant for each channel is carried over from the k = 1 case and is not verified.

## 5. Msub = 4 numerical angle sweep (re-trimmed PhC)

Selected rows (full 67-angle table in `angle_sweep.json`; errors are maxima over 0.231–0.268, degrees):

| θ0 | req. cycles | ρ | mismatch % | max err | RMS err | R_squint | max Δ vs ideal sub. | 5 % HPBW | strongest other lobe dB | replica dB | HPBW(f0) | fixed-phase max err |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.174 | 5.63 | +463 | 1.763 | 0.940 | 4.39 | 1.787 | 0.32 | −12.4 | −16.6 | 6.4 | 0.40 |
| 14 | 0.484 | 2.03 | +103 | 1.078 | 0.567 | 0.94 | 1.145 | 0.33 | −12.4 | −16.6 | 6.6 | 1.15 |
| 20 | 0.684 | 1.43 | +43 | 0.625 | 0.319 | 0.372 | 0.724 | 0.34 | −12.4 | −16.6 | 6.8 | 1.68 |
| 23 | 0.781 | 1.26 | +25 | 0.396 | 0.194 | 0.202 | 0.511 | 0.35 | −12.4 | −16.6 | 6.9 | 1.96 |
| 26 | 0.877 | 1.12 | +12 | 0.163 | 0.067 | 0.072 | 0.295 | 0.35 | −12.4 | −16.6 | 7.1 | 2.26 |
| 28 | 0.939 | 1.04 | +4.5 | 0.049 | 0.025 | 0.020 | 0.149 | 0.36 | −12.4 | −16.6 | 7.2 | 2.47 |
| **29.38 (match)** | 0.981 | 1.000 | 0.0 | 0.145 | 0.083 | 0.055 | 0.046 | 0.37 | −12.4 | −16.6 | 7.3 | — |
| **30** | 1.000 | 0.981 | −1.9 | 0.189 | 0.110 | 0.070 | 0.055 | 0.37 | −12.4 | −16.6 | 7.3 | 2.69 |
| 32 | 1.060 | 0.93 | −7.4 | 0.332 | 0.201 | 0.114 | 0.186 | 0.38 | −12.4 | −16.6 | 7.5 | 2.91 |
| 36 | 1.176 | 0.83 | −16.5 | 0.669 | 0.389 | 0.196 | 0.472 | 0.39 | −12.4 | −16.6 | 7.9 | 3.41 |
| 40 | 1.286 | 0.76 | −23.7 | 1.041 | 0.590 | 0.262 | 0.813 | 0.42 | −12.4 | −16.6 | 8.3 | 3.97 |
| 46 | 1.439 | 0.68 | −31.8 | 1.670 | 0.927 | 0.335 | 1.388 | 0.46 | −12.4 | −16.6 | 9.2 | 4.98 |
| 52 | 1.576 | 0.62 | −37.7 | 2.435 | 1.328 | 0.384 | 2.085 | 0.52 | −12.4 | −16.6 | 10.4 | 6.33 |
| 55 | 1.638 | 0.60 | −40.1 | 2.896 | 1.564 | 0.401 | 2.506 | 0.56 | **−5.4** | −16.6 | 11.2 | 7.22 |
| 61 | 1.749 | 0.56 | −43.9 | 4.094 | 2.152 | 0.416 | 3.598 | 0.67 | **−0.3** | −16.6 | 13.4 | 9.85 |
| 70 | 1.879 | 0.52 | −47.8 | 7.764 | 3.689 | 0.388 | 7.000 | 1.11 | 0.0 | 0.0 | 22.2 | 20.0 |

* Minimum peak power is −0.075 dB at every angle (constant to 1e-6 dB; the validated band has flat amplitude).
* The offset-removed squint spread equals R_squint at every angle, because the re-trim sets the f0 pointing error to zero by construction.
* Range over the sweep: max error 0.049°–7.76°, RMS 0.025°–3.69°. The smallest max error on the 1° grid is 0.049° at 28° (R_squint 0.020), not at the exact match: at the match the intra-subarray fixed-phase squint (ideal-subarray error 0.145° there) dominates. The smallest deviation from ideal subarray TTD is at the exact match (0.046°, equal to the numerical-grid minimum, test-checked).
* The strongest-other-lobe level stays at the 30° value (−12.4 dB) and the strongest replica at −16.6 dB up to 52°; the strongest other lobe exceeds −10 dB from 53° (the ideal-subarray reference also reaches 0 dB somewhere in the sweep, so this is array geometry at wide angles, not a PhC effect; the mechanism was not analysed further). Replicas reach 0 dB from 64°.
* R_squint < 1 (re-trimmed PhC beats fixed phase) over 14°–70°; below ~14° the fixed phase squints less in degrees than the mismatched channel (R_squint 4.4 at 5°).
* The HPBW at f0 grows from 6.4° (5°) to 7.3° (30°) to 22° (70°); the 5 %-HPBW tolerance grows with it.

## 6. Useful steering windows (post hoc, contiguous on the 1° grid)

| criterion | window |
|---|---|
| A: R_squint ≤ 0.20 | **24°–36°** (14 grid points including the exact match) |
| B: max \|θ_PhC − θ_ideal subarray\| ≤ 0.05 HPBW_ideal(f0) | **26°–34°** (10 grid points including the exact match) |
| C: max pointing error ≤ 0.5° | **22°–34°** |
| A and B and C | 26°–34° |
| for context: max error ≤ 1° | 16°–39° |
| for context: R_squint < 1 | 14°–70° |
| strongest other lobe ≤ −10 dB | 5°–52° (not a PhC-specific limit) |

Window B is the Phase 6 5 %-HPBW logic applied to the ideal-subarray reference; HPBW_ideal(f0) is taken from the ideal element-level TTD pattern as in Phase 6. The stricter criteria give a window of about ±4° around 30° (about one HPBW wide, 7.3°). The window is asymmetric around the exact match (29.4°), because the intra-subarray squint adds to the mismatch term on one side and subtracts on the other.

## 7. Analytic-vs-numerical comparison

The analytic law is evaluated per frequency with the sweep's ρ and compared with the numerical θ_peak(f) of the re-trimmed case (every angle satisfies the law's validity condition).

* Maximum deviation over the sweep: **0.43°** (at 70°); median 0.16°.
* It is nearly angle independent, 0.127° at 5° to 0.152° at 34° and 0.19° at 50°, then rises to 0.23° at 60° and 0.43° at 70°. It is never ≤ 0.05° — the law does not reach the 0.05 HPBW level — because the law has no finite-subarray element factor (intra-subarray fixed-phase squint, 0.16° at 30°) and no pointing-curvature term.
* The law predicts the maximum error at 30° as 0.050° against 0.189° numerically (the near-match regime is where it is least useful: the residual is dominated by what the law omits). Away from the match it agrees with the numerical max error to about 0.1–0.15°: at 10° 1.47° vs 1.38°; 20° 0.72° vs 0.63°; 40° 0.92° vs 1.04°; 50° 2.01° vs 2.16°; 70° 7.33° vs 7.76°. Adding the numerical ideal-subarray error (first-order superposition) reduces the deviation at 15–30° (0.05–0.07°) but not at larger angles, so the leftover is not a simple sum. The law is therefore a good guide to the mismatch contribution (right sign and magnitude, within ~0.15° up to 40°) but not a replacement for the full array evaluation near the match.

## 8. Conventional propagation references (one carrier cycle at f0; `conventional_reference.json`)

| medium | group index n_g | L / λ0 | L / a | length relative to W1 | W1 shorter by |
|---|---|---|---|---|---|
| vacuum / air | 1 (dispersion-free) | 1.0000 | 4.008 | 3.932× | 74.6 % |
| ideal homogeneous ε = 12 | √12 = 3.4641 (assumed non-dispersive) | 0.2887 | 1.157 | 1.135× | 11.9 % |
| selected plain W1 PhC | 3.9323 (from 3.9323 a/c per period; flat ±10 % over 14.8 % band) | **0.2543** | **1.0193** (1/(f0 τ) = 1.0193 periods) | 1.000 | – |

Notes:
* For a non-dispersive medium, L = λ0/n_g gives exactly one carrier cycle of delay at f0 (arithmetic checked in the test).
* **The ε = 12 row is an ideal homogeneous-medium propagation reference. It is not a modeled dielectric waveguide.** A real dielectric waveguide has effective and group indices set by its geometry; n_g does not equal √ε. The row is an analytic lower-complexity benchmark, not a feed-network design, and no conventional line was simulated.
* Continuous length selection: air and the uniform medium yes; no integer-period quantum; no Bloch-phase intercept problem (the phase delay equals the group delay at f0 when the medium is non-dispersive, and an exact one-cycle length leaves no residual phase). The W1 simulated channels exist only in whole periods (one-period quantum 0.981 cycle) and require the constant phase trim.
* The final report already compares the W1 with an MPB-only square-lattice rod waveguide (~2.1× fewer wavelengths, i.e. an implied n_g ≈ 1.9 — inferred from that number, not recomputed here). It shows that real waveguides are not at n_g = √12; I have not extended the comparison.
* This compares length at fixed delay only. Lateral width, cladding rows, loss, couplers, bends and bandwidth are not compared.

## 9. What changes in manuscript framing

* The 30° result should be described as a **worked example at a near-match angle**, with the exact match at 29.38° and a ±4° window (several HPBW fractions) for the stated post hoc criteria when only a constant trim is reset. State the window values and that they are post hoc.
* Remove or reword wording that suggests general steering: "supports feasibility for beam-squint reduction at the matched angle" is supported; "scanning architecture" is not. Larger scan ranges would require Msub/k ladders (Section 4) or additional delay increments, which were not simulated.
* Add the match relation as standard array theory with the k/Msub table; call the angle dependence a limitation, not an advantage. Fig. 10 and Fig. 11 can go in the manuscript or the appendix.
* **Footprint:** state the lengths in Section 8. A PhC length advantage over air (3.9× shorter) holds for any slow-wave medium. Relative to the ideal ε = 12 medium, the W1 is only 11.9 % shorter. Remove statements of a large or general footprint advantage over a conventional line. The correct wording: "At the modest selected group index (3.93), a large footprint advantage over a high-index dielectric propagation reference is not demonstrated."
* The analytic squint law can be cited as standard theory with its measured agreement (deviation ≤ 0.15° up to 40°; 0.43° at 70°).

## 10. What does NOT change

* All Phase 4–6 numbers (0.189° at 30°, R_squint 0.070, −12.4 dB, −16.6 dB, 1.9 % mismatch, 0.981 cycle per period, 1.96 λ0 element-level pitch, Monte Carlo robustness).
* The "computational verification + negative-to-positive case study" framing; the element-level negative result; the Msub = 4 hierarchical positive result; the not-novel status of the relation, the architecture and the physics.
* The restrictions: validated band 0.231–0.268, 2-D lossless model, no experiment, no priority claim.

## 11. Limitations

* The sweep uses the same four res-30 channels (one mode at four lengths) at every angle; there is no independent evidence for other angles beyond the model. The res-40 check was done only at 30° in Phase 6 (not repeated).
* Case A assumes the constant phase trim can be reset per angle (a tunable phase shifter or equivalent), which adds a per-angle calibration; the 0.5°/p95 robustness tolerances from Phase 6 were derived at 30° only.
* The 1° grid and the interval edges have 1° resolution; the exact window edges lie within one grid step of the quoted values.
* The isotropic-element model has no real element pattern or mutual coupling.
* Windows depend on post hoc thresholds chosen after seeing the 30° result; they are descriptive.
* k > 1 channels, Msub = 8/16 hierarchies and the match table beyond k = 1 are arithmetic only.
* The ε = 12 reference is non-dispersive and ideal; its length advantage/disadvantage is not a statement about any real waveguide.

## Answers to the scientific questions

1. **Isolated or part of a window?** A window, but narrow: 26°–34° (criteria A, B and C all met), 22°–34° for ≤ 0.5°, 16°–39° for ≤ 1°.
2. **Degradation away from 30°:** max error 0.05° at 28° → 0.49° at 34° → 1.04° at 40° → 2.4° at 52° → 7.8° at 70°; toward lower angles 0.63° at 20°, 1.08° at 14°, 1.76° at 5° (R_squint exceeds 1 below ~13°).
3. **Exact Msub = 4, k = 1 match:** **29.377°** (sin θ = 0.4906).
4. **Other visible combinations:** see Section 4: Msub = 2 (k = 1, 78.85°), Msub = 4 (k = 2, 78.85°), Msub = 8 (14.20°, 29.38°, 47.38°, 78.85°), Msub = 16 (eight angles from 7.04° to 78.85°). Msub = 1: none. Only k = 1 channels were simulated.
5. **Scanning with re-trimming only?** Only in a narrow band about 29.4°; beyond ±5° the mismatch term dominates. It still beats fixed phase for 14°–70° in R_squint, but with pointing errors of 1–8° at the ends, i.e. not "meaningful scanning" in the HPBW-fraction sense.
6. **Description:** between a fixed-beam worked example and a narrow-angle design. Safest wording: a worked example at a near-matched angle with a narrow useful window.
7. **Footprint:** one carrier cycle needs 0.254 λ0 in the W1 vs 1.0 λ0 in air (3.93× shorter, 74.6 %) and vs 0.289 λ0 in an ideal ε = 12 homogeneous medium (11.9 % shorter, 1.135×).
8. **Remove all footprint claims?** Remove the large-advantage language; keep the factual lengths with the caveats. The air comparison is trivial (any slow-wave medium); the ideal-high-index comparison is weak.
9. **Safest arXiv framing:** unchanged in kind (computational verification + negative-to-positive case study); the angle dependence and the weak footprint comparison make the claims more modest, not different.

## Classifications

* **Angle result: C. NARROW ANGLE WINDOW** (±4° at the stricter criteria, 22°–34° at 0.5°; not broadly scannable by re-trimming alone).
* **Footprint comparison: C. WEAK PhC ADVANTAGE** (length only; 11.9 % shorter than the ideal ε = 12 reference, 74.6 % shorter than air; width and other costs not evaluated, so "no demonstrated advantage" is also defensible for area).
