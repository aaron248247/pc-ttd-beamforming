# pc-ttd-beamforming

Computational study of **dispersion-engineered photonic-crystal (PhC) line-defect waveguides as true-time-delay
elements for broadband phased-array beamforming**, including a conventional-waveguide baseline.
The answer is not assumed to be favorable. See [docs/research_plan.md](docs/research_plan.md).

## Status

| Phase | State |
|---|---|
| 0 Research design | done ([docs/research_plan.md](docs/research_plan.md)) |
| 0.5 Geometry selection (MPB) | done and approved — [docs/phase0p5_report.md](docs/phase0p5_report.md); decisions in [docs/decisions.md](docs/decisions.md) |
| 1 FDTD verification ladder | done and accepted — [docs/phase1_report.md](docs/phase1_report.md) |
| 2 Finite W1 guided-wave verification | gate passed with a documented restriction — [docs/phase2_report.md](docs/phase2_report.md) |
| 3 Dispersion engineering / delay-bandwidth trade study | gate passed with restrictions (candidate for Phase 4: plain W1) — [docs/phase3_report.md](docs/phase3_report.md) |
| 4–7 | not started (4: TTD channel family, 5: array, 6: robustness, 7: documentation) |

Phase 3 result (3 designs, a tradeoff illustration): in the FDTD-validated range 0.231-0.268 the plain W1 gives a flat n_g of about 3.9 over 14.8% fractional bandwidth; engineered first rows (D1, D2) raise n_g by 7-14% but narrow the +-10% band to 8.4% and 5.6%. No beam-squint claim is made yet (Phases 4-5 not started).

## Environment

Uses the `MEEP-GPU` conda env (MEEP 1.33.0-beta with MPI, MPB via `meep.mpb`, NumPy 1.24.3, SciPy 1.11.1, Matplotlib 3.10.6,
pandas 2.0.3, pytest 7.4.0). Versions are recorded in `environment/env_versions.json`; `environment/environment.yml` is the export.

## Reproduce Phase 3

```bash
PY=~/miniconda3/envs/MEEP-GPU/bin/python
cd experiments/phase3
$PY 3a_mpb_designs.py && $PY 3c_mpb_res.py && $PY 3d_mpb_analyze.py       # MPB (about 20 min)
for d in "D0 1.0 0.0" "D1 1.2 -0.05" "D2 1.2 -0.10"; do set -- $d; $PY 3b_fdtd.py $1 30 $2 $3 4 8 12 & done; wait   # ~1 h; res 40 confirmation: ~2-3 h
for d in D0 D1 D2; do $PY 3e_analyze.py $d 30; $PY 3e_analyze.py $d 40; done
$PY 3g_tables.py && $PY 3f_figures.py
```

## Reproduce Phase 1 (verification ladder)

```bash
PY=~/miniconda3/envs/MEEP-GPU/bin/python
$PY -m pytest                                        # all tests (pure-Python, MPB, MEEP markers)
cd experiments/phase1
$PY 1a_vacuum.py && $PY 1b_dielectric.py             # minutes
$PY 1d_mpb_directional.py                            # a few minutes (MPB, res 32/64/128)
$PY 1c_bragg.py                                      # ~25 min (TMM vs MEEP, many studies)
$PY 1d_pc2d.py                                       # ~1 h (2-D single-length sensitivity studies)
$PY 1d_attenuation.py                                # ~1.5 h (N = 6..16 x ny = 24..48, minimum run time 20000)
$PY 1d_termination.py && $PY 1d_termination.py 12 16 # ~25 min each
$PY 1d_edges.py && $PY report.py && $PY make_figures.py
```

Methods and derivations: [docs/phase1_methods.md](docs/phase1_methods.md).

## Reproduce Phase 0.5

```bash
PY=~/miniconda3/envs/MEEP-GPU/bin/python
$PY -m pytest                                   # 30 tests (6 need MPB)
cd experiments/phase0p5
$PY 01_bulk_gap_scan.py && $PY 01b_bulk_bands.py
$PY 02_defect_sweep.py                          # H and R; add "H H_r030 extra=20" to re-solve H configs
$PY 02_defect_sweep.py H H_r030 extra=20
$PY 02_defect_sweep.py supp                     # R with wider supercell
$PY 02c_analyze_sweep.py defect_sweep.json && $PY 02c_analyze_sweep.py defect_sweep_supp.json
$PY 03_select.py && $PY 05_convergence.py H H_r030 R && $PY 04_figures.py && $PY 06_report.py
```

(`02_defect_sweep.py` uses 10 worker processes; expect ~10–20 min per sweep on 12 cores.)

## Layout

```
src/pcttd/        geometry.py (W1 supercells)  bands.py (gaps, modes, windows, GBP)  mpb_bands.py (MPB drivers)  io.py
experiments/      one folder per phase; scripts numbered in run order
tests/            analytic and regression tests; MPB closed-form checks marked `mpb`
docs/             research plan, pre-declared criteria (with dated revisions), reports
results/          machine-readable outputs with provenance; figures/ is regenerated from results/ only
```

## Data-integrity rule

Every number in the docs is either read from `results/` by a script (`06_report.py` tables) or copied from those tables.
Planned-but-not-run results are labeled planned.

## Conventions

American English spelling in all prose, comments, and documentation. Normalized units: a = 1, c = 1, frequency f = omega a / (2 pi c) = a/lambda.
