# pc-ttd-beamforming

Computational study of **dispersion-engineered photonic-crystal (PhC) line-defect waveguides as true-time-delay
elements for broadband phased-array beamforming**, including a conventional-waveguide baseline.
The answer is not assumed to be favorable. See [docs/research_plan.md](docs/research_plan.md).

## Status

| Phase | State |
|---|---|
| 0 Research design | done ([docs/research_plan.md](docs/research_plan.md)) |
| 0.5 Geometry selection (MPB) | done and approved — [docs/phase0p5_report.md](docs/phase0p5_report.md); decisions in [docs/decisions.md](docs/decisions.md) |
| 1 FDTD verification ladder | in progress |
| 2–7 | not started |

No performance claim about delay, bandwidth, or beam squint is made anywhere in this repository yet: those need the FDTD phases.

## Environment

Uses the `MEEP-GPU` conda env (MEEP 1.33.0-beta with MPI, MPB via `meep.mpb`, NumPy 1.24.3, SciPy 1.11.1, Matplotlib 3.10.6,
pandas 2.0.3, pytest 7.4.0). Versions are recorded in `environment/env_versions.json`; `environment/environment.yml` is the export.

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
