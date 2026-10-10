# Executive summary

*Written for: a technically literate reviewer (hiring manager, research scientist, professor).*

**Research question.** Radar and communication antenna arrays steer a beam by delaying the signal sent to each antenna. If the delay is made with a fixed phase shift, the beam points to a slightly different angle at each frequency (beam squint), which limits bandwidth. True-time delay fixes that, but delay lines are long. Photonic-crystal waveguides slow light down, so could they make short, wideband delay lines? I set out to find out, without assuming the answer would be yes.

**What I built.** A verified computer-simulation workflow (MEEP time-domain solver and MPB band-structure solver, plus a Python library with tests) that: checks the simulator against exact answers, extracts the delay and complex transmission of finite photonic-crystal waveguide sections, builds a four-channel delay family from separately simulated sections, and uses those simulated responses as the weights of a phased array to compute beam patterns and sensitivity to errors. Everything is traceable to result files and covered by 195+ automated tests.

**What was verified.** The solver reproduces known delays (vacuum, dielectric), an analytic multilayer stack, and the band structure of the crystal, with measured second-order convergence. The waveguide group index agrees with the band-structure calculation to about 1 % over a 14.8 % wide frequency band. That band is the range I validated, not necessarily the physical limit.

**What was found.**
* A plain waveguide (one missing row of holes) gives about one carrier cycle of delay per lattice period, accurate to 0.04 % across four channels, over that 14.8 % band. I also tried modifying the neighboring holes to slow the light more: it bought 7–14 % more delay for 43–62 % less bandwidth, so I rejected it. The final delay line is therefore the plain waveguide, not an engineered one.
* The simulated channels have a large constant phase offset (up to 2.7 rad) that has to be corrected with one fixed phase adjustment per channel. After that, the delay behavior matches ideal true-time delay closely.
* **What failed.** Using one channel per antenna did not give a practical array: because one period already delays about a full cycle, the antennas had to sit about two wavelengths apart, which produces extra beams (grating lobes) as strong as the main beam.
* **The fix.** The same delay step happens to equal what a group of four half-wavelength-spaced antennas needs when steered to 30° (1.000 cycle needed, 0.981 supplied, a 1.9 % mismatch). With one channel feeding each group of four antennas, the worst beam-pointing error over the band is 0.19° against 2.69° for a fixed phase shifter. The extra beams drop to −12.4 dB (ordinary sidelobes) with remaining replica lobes at most −16.6 dB. The match is a coincidence of this geometry, frequency and angle, not a design rule.
* Simulated error budgets (system-level, not fabrication) show the pointing stays below 0.5° for about ±7.6° of phase error or ±1.2 % of delay error without recalibration.

**Final architecture.** Four-element subarrays, one photonic-crystal delay channel plus one constant phase trim per subarray.

**What remains.** The model is two-dimensional and lossless, with no fabrication disorder, coupling between channels, bends, realistic couplers or temperature effects; there is no experiment. The constant phase trim is assumed realizable, and how is not addressed. Only one steering angle was studied. I have not compared against a conventional delay line beyond a band-structure-only estimate (the photonic crystal needs about half the length), and I have not surveyed the literature, so I make no novelty claim. Natural next steps: add loss and disorder, simulate a realistic coupler, and compare against a conventional line at the subarray level.

**Honest process notes.** Several things went wrong and were fixed in the open: a geometry sign bug (twice), an under-converged first batch of simulations that I discarded, a post hoc choice of the upper band edge, and a post hoc change of the delay estimator. They are all documented, because they are part of why the numbers can be trusted.

Further reading: `docs/final_report.md` (full report), `docs/manuscript_draft.md` (paper draft), `docs/phase7_audit.md` (claim audit), `docs/final_tables.md` (tables).
