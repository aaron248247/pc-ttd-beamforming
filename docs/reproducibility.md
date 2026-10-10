# Reproducibility

## Environment

* Interpreter: `~/miniconda3/envs/MEEP-GPU/bin/python` (Python 3.10.19, Linux/WSL2). Recorded in `environment/env_versions.json`; export in `environment/environment.yml` (conda env name `MEEP-GPU`).
* Libraries: MEEP 1.33.0-beta (with MPI), MPB via `meep.mpb`, NumPy 1.24.3, SciPy 1.11.1, Matplotlib 3.10.6, pandas 2.0.3, pytest 7.4.0. (`pyproject.toml` declares `pcttd` 0.0.1, Python >= 3.10, and puts `src` on the pytest path; the package is not pip-installed, so scripts add `src` to `sys.path` themselves.)
* Nothing requires a GPU; the environment name is historical. MEEP and MPB are required only for Phases 0.5-4 and for the MEEP/MPB-marked tests.

```bash
PY=~/miniconda3/envs/MEEP-GPU/bin/python
cd ~/MEEP/pc-ttd-beamforming
$PY -m pytest -q          # expected: 206 passed, 0 failed, 0 skipped (195 at Phase 6 + 11 Phase 7); 1.5-2.5 minutes
```

## What is precomputed and tracked

| data | tracked | how to regenerate |
|---|---|---|
| processed results, handoffs, figures, tables for Phases 0.5-7 (`results/`, `figures/`, `docs/*_tables.md`) | yes | scripts below |
| raw plane-field dumps `results/phase2/raw`, `phase3/raw`, `phase4/raw`, `phase4/raw_attempt1` (about 134 MB) | **no (git-ignored)** | rerun the FDTD (hours); the processed `.npz/.json` that depend on them are tracked, so analysis, figures and Phases 5-7 do not need them |
| Phase 4 first (discarded) attempt | processed files tracked (`results/phase4/attempt1_*`), raw not | for transparency only |

## Run sequence and cost

| phase | scripts (in `experiments/phaseX/`, run in this order) | cost |
|---|---|---|
| 0.5 | `01_bulk_gap_scan.py`, `01b_bulk_bands.py`, `02_defect_sweep.py` (+ `H H_r030 extra=20`, `supp`), `02c_analyze_sweep.py`, `03_select.py`, `05_convergence.py`, `04_figures.py`, `06_report.py` | MPB, about 10-20 min per sweep on 12 cores |
| 1 | `1a_vacuum.py`, `1b_dielectric.py`, `1d_mpb_directional.py`, `1c_bragg.py`, `1d_pc2d.py`, `1d_attenuation.py`, `1d_termination.py`, `1d_edges.py`, `report.py`, `make_figures.py` | minutes to about 1.5 h per script |
| 2 | `2a_mpb.py`; `2b_fdtd.py <tag> <res> <f1> <s> <N...> [--key=value]` then `2c_analyze.py <tag> <res> [N...]`; `2d_pulse.py`, `2e_convergence.py`, `2f_sampling.py`, `report.py`, `make_figures.py` | FDTD hours. **The README did not list Phase 2 commands. The production run is inferred from the stored configuration (`results/phase2/2c_tol7_res40.json`, tag `tol7`): `2b_fdtd.py tol7 40 1.0 0.0 6 12 24 --decay_tol=1e-7` (and res 30), then `2c_analyze.py tol7 40`.** The exact original invocation was not recorded |
| 3 | `3a_mpb_designs.py`, `3c_mpb_res.py`, `3d_mpb_analyze.py`; `3b_fdtd.py` for D0/D1/D2 (res 30, 40); `3e_analyze.py`, `3g_tables.py`, `3f_figures.py` | MPB about 20 min; FDTD about 1 h (res 30, parallel), 2-3 h (res 40) |
| 4 | `4a_fdtd.py 30 --min_after=3500` and `... 40 ...`; `4b_analyze.py 30/40`; `4c_resolution_handoff.py --handoff`; `4e_tables.py`, `4d_figures.py` | about 2 h (res 30), 3.5 h (res 40), 8 parallel runs |
| 5 | `5a_array.py`, `5b_figures.py`, `5c_tables.py`, `5d_handoff.py` (needs `results/phase4/phase4_handoff.npz`) | seconds to minutes |
| 6 | `6a_subarrays.py`, `6b_robustness.py` (about 7 min on 10 cores), `6c_figures.py`, `6d_tables.py`, `6e_handoff.py` | minutes |
| 7 | `7a_final_results.py`, `7b_figures.py`, `7c_tables.py`, `7d_handoff.py` | about 1 minute total, no FDTD |

Phases 2-4 are the expensive ones; none of them is needed to regenerate the final figures, tables or handoff.

## Regenerate the final deliverables

```bash
cd experiments/phase7
$PY 7a_final_results.py     # results/final/final_results.json (149 values with sources)
$PY 7b_figures.py           # figures/final/fig01..fig09 (.png and .pdf)
$PY 7c_tables.py            # docs/final_tables.md
$PY 7d_handoff.py           # results/final/final_handoff.json
```

The final report, manuscript, executive and interview summaries are prose documents; their numbers were taken from `final_results.json`/`final_tables.md` and the Phase 4-6 reports and are checked against it by `tests/test_phase7.py` (see below).

## Tests

`tests/test_phase7.py`: final results match the Phase 6 handoff and recomputation from the Phase 4 handoff; every figure in the handoff list exists; the final handoff references existing files and records the authoritative commits; Phase 4/5/6 result files are unchanged (sha256 recorded in the final handoff); every number quoted in `docs/final_tables.md` for the headline metrics comes from `final_results.json`; README run commands point to existing scripts. Phase 1-6 tests are unchanged.

## Known omissions and caveats

* Raw FDTD dumps for Phases 2-4 are not in the repository (size). Reproducing Phase 2-4 numbers requires rerunning the FDTD; results should agree to the documented resolution-pair differences, not bit-for-bit (MEEP with MPI is not guaranteed bitwise reproducible across machines).
* Phase 2 command lines are inferred (above).
* Provenance records in earlier result files carry `git_dirty`/`code_dirty` flags and commit IDs from the time of each run; some runs were made with uncommitted code that differs from the later committed code only as documented in the phase reports.
* Monte Carlo (Phase 6) is reproducible from seed 20261010; this is tested.
