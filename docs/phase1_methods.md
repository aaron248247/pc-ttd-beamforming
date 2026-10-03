# Phase 1 methods and derivations

Units: lattice constant a = 1, c = 1; frequency f = omega a/(2 pi c) = a/lambda; delays in a/c.

## 1. Fourier-sign convention and the sign of group delay

MEEP evolves fields with time dependence exp(-i omega t) and accumulates the DFT  F(omega) = ∫ f(t) exp(+i omega t) dt, which isolates the
exp(-i omega t) component of a real signal. For a right-going wave E(x,t) = Re[E0 exp(i(kx - omega t))] this gives F(omega, x) = E0 exp(ikx).
For two reference planes x1 < x2 separated by L, the transfer function is  H(omega) = F(omega, x2)/F(omega, x1) = exp(ikL), so phi = arg H = kL
and the group delay of a wavepacket is  tau_g = L/v_g = L dk/d omega = **+ d phi / d omega**.
(With exp(+i omega t) the sign would be -d phi/d omega, the form quoted in the project brief.)
This was **not taken on faith**: tests/test_fdtd_1d.py and results/phase1/1a_vacuum.json apply both signs to the vacuum FDTD data and show that +1 gives
+L/c while -1 gives -L/c (`extract.MEEP_SIGN = +1`).

## 2. Reference planes and de-embedding

Everything is a ratio of DFT fields, so source amplitude cancels. Two situations:
* *Propagation* (Phases 1A/1B): H = F(x2)/F(x1) inside the same medium. The delay is the full L/c (vacuum) or nL/c (index n).
* *Scattering structure* (1B slab, 1C stack, 1D crystal): the structure is removed in a normalization run with the same monitor planes.
  Then t_ref = F_B(structure)/F_B(empty) is the transmission **relative to propagating through the same total length of vacuum**; its group delay
  is the *excess* delay tau_g - D/c (D = structure thickness). TMM's face-to-face t equals t_ref multiplied by exp(+i k0 D); `tmm.vacuum_referenced_t` applies the
  inverse so the two are compared on the same footing. The reflection coefficient is r_ref = (F_A(structure) - F_A(empty))/F_A(empty), because the empty run
  contains only the incident wave at monitor A.

## 3. Yee-grid (numerical) dispersion in 1-D

The leapfrog scheme for a homogeneous medium of index n gives the exact discrete dispersion  sin(k dx/2) = (n/S) sin(omega dt/2),  S = dt/dx.
Writing a = n omega dx/2 and theta_t = omega dt/2: sin(theta_x) = a(1 - theta_t^2/6 + ...) so theta_x = a[1 + (a^2 - theta_t^2)/6], i.e.

    k = n omega [1 + (omega dx)^2 (n^2 - S^2)/24]   ->   tau_g = n L [1 + (omega dx)^2 (n^2 - S^2)/8].

The delay error of a clean FDTD run is therefore positive, second order in dx, and grows as omega^2. `extract.yee_1d_delay_phase` evaluates the exact formula;
the measured vacuum, dielectric and slab results are compared with it (extraction check) and with the ideal value (discretization error).

## 4. Group-delay extraction

* `group_delay_phase`: np.unwrap of arg H, then central differences d phi/d omega. Valid only if the phase step between adjacent samples is below pi:
  tau < 1/(2 df). Aliasing cannot be detected from the wrapped phase alone (a step of 4 pi looks like zero), so `assert_sampling_adequate` imposes the
  condition *a priori* from the known maximum delay (TMM or the structure size) with a safety factor.
* `group_delay_ratio`: tau = arg(H[i+1] conj(H[i-1]))/(omega[i+1]-omega[i-1]); needs no global unwrapping but spans two samples (limit tau < 1/(4 df)).
* `group_delay_derivative`: Im(H'/H); simple but its truncation error tau (d omega tau)^2/6 limits it to d omega * tau << 1. Kept as a tested cross-check only.
* Near resonances the finite-difference estimator itself has truncation error; whenever FDTD is compared with TMM, **the same estimator is applied to both** and the estimator's own
  truncation (TMM analytic derivative vs sampled estimator) is reported separately.
* Where |t| is vanishingly small the phase is noise; delays are only compared where T > 1e-3.

## 5. Transfer-matrix method (lossless 1-D multilayer, normal incidence)

See the module docstring of `src/pcttd/tmm.py` for the matrix derivation (interface matrix D = (1/2)[[1+rho, 1-rho],[1-rho, 1+rho]], propagation diag(e^{-i delta}, e^{+i delta})).
Analytic checks used: Fresnel coefficients; vacuum layer t = e^{ikd}; half-wave layer transparent; T + R = 1 to rounding for random lossless stacks; reciprocity under layer reversal;
quarter-wave mid-gap reflectance ((n_in n_H^{2N} - n_out n_L^{2N})/(n_in n_H^{2N} + n_out n_L^{2N}))^2; infinite-lattice band edges from |cos K Lambda| = 1 with
cos K Lambda = cos d_a cos d_b - (1/2)(n_a/n_b + n_b/n_a) sin d_a sin d_b, which for the quarter-wave stack are f0[1 -+ (2/pi) asin((n_H-n_L)/(n_H+n_L))];
stop-band transmission decay exp(-2 acosh|cos K Lambda| N); the Hartman saturation of the stop-band delay.

## 6. FDTD pitfalls found and guarded against (all with regression tests)

1. `add_dft_fields` takes (fcen, df, nfreq) *or* an explicit frequency list. Passing (fmin, fmax, nfreq) silently yielded a frequency axis 1.25x too compressed
   and a delay 25% off. The driver passes explicit lists and the test asserts the DFT object's frequencies equal the requested ones.
2. The default DFT decimation is chosen from the Nyquist rate; decimation 1 is used (the automatic choice changed the delay at the 4e-4 level in vacuum).
3. `stop_when_dft_decayed` passes trivially before the wave reaches a distant monitor (a DFT that has not started has not "changed"). The driver sets the minimum run time to the
   transit time to the farthest monitor at the slowest speed present, plus a margin.
4. Point-monitor Yee-grid output is empty for a zero-size volume; the interpolated value is the average of two Yee nodes, i.e. a real factor cos(k dx/2), which does not change phase.
5. MEEP rounds cell sizes to whole pixels. The 2-D periodic height sqrt(3) a is imposed exactly by choosing the pixel count across it (resolution = ny/sqrt(3)); with
   an integer resolution the cell was silently strained by ~1%.
6. R + T = 1 is a sensitive diagnostic of ring-down termination (it was 1e-2 at decay tolerance 1e-8 for the 8-period stack and 6e-7 at 1e-14), much more so than T itself.
7. `fwidth` of a Gaussian source is the inverse temporal width: the spectral standard deviation is fwidth/(2 pi). Outside roughly 4 sigma the field is dominated by a numerical/truncation
   floor and ratios of fields are meaningless; the study records the band-edge source amplitude relative to the peak ("coverage") and the criterion derived from it.

## 7. Finite 2-D crystal vs MPB: what is and is not comparable

* The FDTD cell has height sqrt(3) a, periodic in y, with a y-uniform Hz line source: an exact normal-incidence plane wave. The crystal is then the rectangular two-hole cell that MPB solves at k = (k_x, 0).
* The source and structure are mirror symmetric about y = 0, so only y-even-Hz modes are excited. The relevant prediction is therefore **the frequency interval with no even-parity band** for any k_x in [0, 0.5]
  (the forward-mode counting of Phase 0.5 with the even-parity class), not the omnidirectional bulk gap. An odd-parity band can lie inside the interval without affecting transmission.
* The MPB gap edges are properties of an infinite crystal. A finite section of N columns has (i) a stop-band edge that approaches the band edge from inside as N grows (the transmission at the
  edge of a finite stack is intermediate, and flat bands near the edge ring for a long time), (ii) Fabry-Perot fringes in the passband from reflections at the vacuum/crystal faces
  (impedance mismatch plus the termination choice), and (iii) a deep-gap transmission that decays exponentially in N with the evanescent Bloch rate, which MPB cannot supply directly
  (it gives real k only). The comparison is therefore made on edge positions (at two attenuation levels) and on their convergence with N, resolution and termination.
