# Research plan (Phase 0, as approved with modifications)

> **Note (Phase 7).** This is the approved Phase 0 text. The question uses "dispersion-engineered"; in the study the engineered first-row variants (D1, D2) were tested and rejected in Phase 3, so the final delay line is the plain W1. See `docs/phase7_audit.md`.


**Question.** Can dispersion-engineered photonic-crystal (PhC) line-defect waveguides serve as true-time-delay (TTD)
elements for broadband phased-array beamforming, and what is the delay / bandwidth / dispersion / loss trade? Mixed or
negative results are acceptable if rigorously shown.

**Added by review.**
1. *Scale invariance first.* All PhC work is in normalized units (`w = a/lambda`). No physical embodiment is chosen yet.
   After Phase 3 (delay, dispersion, bandwidth, geometry requirements known) at least two mappings are compared —
   (A) microwave/mm-wave EBG/PhC and (B) optical PhC delay element inside an RF-photonic beamformer — on physical size,
   materials, interpretation of bandwidth and modeling limits, and only then is one judged the most defensible.
2. *Conventional-waveguide baseline is core.* The question "what does the PhC buy relative to a conventional waveguide
   with the same target delay?" is answered by comparing: length for a given delay, group delay, usable fractional
   bandwidth, delay ripple, GDD, transmission/reflection, and sensitivity. The footprint-reduction vs
   bandwidth/dispersion/reflection penalty is quantified. The conventional line may win in some regimes.
3. *Three classes of thresholds.* Physical requirements (derived, derivation shown); numerical verification tolerances;
   engineering design targets. Only the first two are pass/fail for the science. Squint reduction is a *continuous*
   metric, `R_squint = max_f|theta_PhC(f)-theta_target| / max_f|theta_phase(f)-theta_target|`, reported whatever its
   value; `R_squint <= 0.2` is a stretch engineering target, not a success criterion. The project succeeds if it
   determines rigorously when, why and by how much PhC TTD does or does not improve broadband steering.
4. *Verification ladder (kept).* vacuum -> homogeneous dielectric -> analytic 1-D Bragg/TMM -> 2-D PhC vs MPB ->
   W1 FDTD-vs-MPB delay/dispersion; plus independent pulse-peak, centroid and length-scaling delay checks. The Fourier-sign
   convention is *verified numerically* on the vacuum and homogeneous cases (known positive delay), not just derived.

**Phases.** 0 plan -> 0.5 geometry selection (MPB) -> 1 baseline PhC + FDTD verification -> 2 line-defect transfer
function and group delay -> 3 dispersion / slow-light study (+ physical-mapping comparison) -> 4 delay-line family (4 then 8
channels) -> 5 phased-array model with measured transfer functions -> 6 robustness / error propagation -> 7 documentation.
Each phase ends with a gate; no phase starts before its predecessor is reviewed.

**Analytic (no MEEP):** TMM benchmark, extraction/unwrapping library, cascading, array factor, squint metrics, tolerance
propagation. **Needs MEEP/MPB:** spectra, S-parameters, pulse checks, band structures, v_g, resolution/PML studies.
