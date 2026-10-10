# Interview / portfolio summary

Short, technically accurate answers. Numbers are from `docs/final_tables.md`.

**1. What problem were you solving?** Beam squint: with fixed phase shifters a steered beam drifts with frequency (2.69° peak error over a 14.8 % band for a 16-element array at 30°). True-time delay removes it, but delay lines are bulky. I studied, in simulation, whether photonic-crystal waveguide delay lines could do it and what the trade-offs are.

**2. Why photonic-crystal waveguides?** They slow the group velocity (here n_g about 3.9), so a given delay needs less length. The question was whether that comes with usable bandwidth and low dispersion, which I tested rather than assumed.

**3. Why W1?** A W1 is a single missing row in a hole lattice: the simplest line defect, with a guided band inside the gap and a flat group-index region. I screened a hole lattice (TE) against a rod lattice (TM) with a pre-declared rule; the hole lattice won on slow-light potential and cost, by a moderate margin. I also tried engineered first rows and found the plain W1 better for bandwidth.

**4. How did you verify MEEP?** A ladder: vacuum and homogeneous dielectric (delay equals L/c and nL/c plus the exact Yee-grid correction, second-order convergence), a Bragg stack against the transfer-matrix method, a 2-D crystal stop band against MPB, and finally the W1 group index against MPB (median bias +0.3 %, max about 1 % in the validated band). I also tested the sign convention numerically, found and fixed several extraction pitfalls (DFT frequency axis, premature termination, derivative-estimator truncation), and kept regression tests for each.

**5. What is true time delay?** A phase that is linear in frequency, φ = 2πfτ, i.e., a real delay τ. It moves the beam consistently across frequency. A constant phase is not a delay.

**6. Why does fixed phase produce beam squint?** A fixed phase per element corresponds to a delay that is correct at one frequency only. The beam angle satisfies sinθ = (f0/f) sinθ0, so it moves off target as f changes (+2.69° to −2.26° over the band).

**7. Why did the Phase 5 array fail practically?** Each lattice period delays about one carrier cycle (0.981), but a half-wave array at 30° needs 0.25 cycle per element. Using one channel per element forces a pitch of 1.96 λ0, which makes grating lobes as strong as the main beam. The delay lines themselves tracked ideal TTD to 0.05°; the limitation is the delay granularity.

**8. Why did Msub = 4 work?** A four-element λ0/2 subarray steered to 30° needs 4·(λ0/2)·sin30°/c = 1.000 cycle between adjacent subarrays; the PhC period gives 0.981 (−1.89 %). The remaining mismatch costs about ±0.05° over the band. The match is a coincidence of n_g·f0 for this device and angle; other subarray sizes mismatch by 96 % or more.

**9. Why is constant phase trim still needed?** The Bloch phase of the guided mode gives each channel a frequency-independent phase error (up to 2.7 rad). A constant phase has no group delay, so it does not hurt the delay slope but it points the beam wrongly (6.7° off) unless it is trimmed. The trim must also be re-chosen for the required step: the Phase 4 trim used unchanged gives 0.73°.

**10. What are the main numerical limitations?** 2-D lossless model; the validated band 0.231–0.268 with an upper edge chosen post hoc; noise-limited pointwise delay derivatives (so delay claims rest on band fits); only two resolutions for convergence; four channels of one mode; one steering angle; robustness is mathematical perturbation, not disorder simulation; the first channel-family run was discarded as under-converged and the estimator was changed after seeing data.

**11. What result are you most confident in?** The delay increments of the four simulated channels (0.04 % from target by band fit, amplitude imbalance 0.023 dB, cascade check 0.006 rad, resolution pair agreement) and the pointing results derived from them. Least confident: anything about real devices (loss, disorder, couplers).

**12. What would you do next with another month?** (a) Do the literature check and position the contribution. (b) Add propagation loss and a 3-D slab or a physically realizable phase trim. (c) Simulate a coupler and a small disorder ensemble. (d) Compare against a conventional delay line at the subarray level, including footprint and loss. (e) Study other steering angles and the two-period channel for Msub = 8.

**13. What would need to change for experiment/hardware?** A 3-D slab design (out-of-plane loss, light-line), realistic input/output couplers, a physical implementation of the phase trim, fabrication tolerance modeling (disorder strongly affects slow light), thermal behavior, channel isolation, and a choice of physical embodiment (microwave/mm-wave PhC or optical), noting that the one-cycle-per-period argument applies only when the PhC delays the radiated carrier, not an RF-modulated optical carrier.
