## T1. Bulk gap (declared rule, res 32)

| lattice | pol. | gap n|n+1 | r_base/a | w_lo | w_hi | rel. width |
|---|---|---|---|---|---|---|
| triangular | TE (Hz) | 1 | 0.400 | 0.246 | 0.405 | 0.490 |
| square | TM (Ez) | 1 | 0.190 | 0.289 | 0.430 | 0.391 |


Relative gap vs radius, triangular (n=1): r=0.200:0.065, r=0.240:0.145, r=0.280:0.233, r=0.320:0.326, r=0.360:0.417, r=0.385:0.467, r=0.395:0.483


## T2. Gates and sweep coverage

| config | lattice | r/a | M | extra bands | solved | unresolved (excluded) | G1 | G2 | G3 | G4 |
|---|---|---|---|---|---|---|---|---|---|---|
| H | triangular | 0.400 | 5 | 20 | 14 | 11 | PASS | PASS | PASS | PASS |
| R | square | 0.190 | 4 | 8 | 25 | 0 | PASS | PASS | FAIL | PASS |
| H_r030 | triangular | 0.300 | 5 | 20 | 25 | 0 | PASS | PASS | PASS | PASS |
| R_M6 | square | 0.190 | 6 | 8 | 25 | 0 | PASS | PASS | PASS | PASS |


## T3. Criterion metrics

| config | C1 gap | C2 SM frac | C3 strip | C4 GBP(n≥10) | GBP(n≥5) | GBP(n≥20) | C5 px/period | C6 frac | max n_g (singular at band extrema) |
|---|---|---|---|---|---|---|---|---|---|
| H | 0.490 | 0.154 | 0.9980 | 0.0483 | 0.1583 | 0.0280 | 5413 | 1.00 | 489 |
| R | 0.391 | 0.385 | 0.9646 | 0.0313 | 0.0589 | 0.0000 | 24931 | 1.00 | 393 |
| H_r030 | 0.280 | 0.138 | 0.9930 | 0.1945 | 0.1945 | 0.1373 | 9623 | 1.00 | 700 |
| R_M6 | 0.391 | 0.385 | 0.9640 | 0.0312 | 0.0583 | 0.0000 | 36011 | 1.00 | 181 |


## T4. Decision per pairing (pre-declared rule)

| pairing | eligible | per-criterion winners | selected | rule |
|---|---|---|---|---|
| declared: H vs R | {'H': True, 'R': False} | C1:H, C2:R, C3:tie, C4:H, C5:H, C6:tie | H | only eligible candidate |
| H vs R_M6 (R with converged supercell) | {'H': True, 'R_M6': True} | C1:H, C2:R_M6, C3:tie, C4:H, C5:H, C6:tie | H | most criteria won |
| H_r030 vs R (suppl. radius) | {'H_r030': True, 'R': False} | C1:R, C2:R, C3:tie, C4:H_r030, C5:H_r030, C6:tie | H_r030 | only eligible candidate |
| H_r030 vs R_M6 | {'H_r030': True, 'R_M6': True} | C1:R_M6, C2:R_M6, C3:tie, C4:H_r030, C5:H_r030, C6:tie | H_r030 | tie on count; tie-break C4 |


All pairings select the same family: **True** (H).


## T5. Plain W1 and best-GBP design per configuration

| config | design | f1/s | max modes in gap | SM-strict(even) frac | SM even-only frac | flat-window n_g | its frac. BW | its GBP | strip frac | edge frac (window max) | n_g at window center [min–max] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H | plain W1 | 1/+0.00 | 3 | 0.106 | 0.237 | — | — | — | 0.9995 | 2.1e-05 | 3.92 [3.92–21] |
| H | best GBP(≥10) | 0.8/-0.10 | 7 | 0.132 | 0.254 | 10.2 | 0.005 | 0.048 | 0.9967 | 2.1e-04 | 4.06 [3.97–489] |
| R | plain W1 | 1/+0.00 | 1 | 0.337 | 0.337 | 10.7 | 0.003 | 0.029 | 0.9727 | 1.8e-03 | 2.03 [1.79–56] |
| R | best GBP(≥10) | 0.9/+0.10 | 2 | 0.361 | 0.361 | 10.1 | 0.003 | 0.031 | 0.9695 | 2.9e-03 | 1.90 [1.58–98] |
| H_r030 | plain W1 | 1/+0.00 | 3 | 0.107 | 0.250 | 63.4 | 0.001 | 0.040 | 0.9978 | 1.1e-03 | 3.87 [3.87–4] |
| H_r030 | best GBP(≥10) | 0.9/+0.10 | 3 | 0.114 | 0.239 | 19.6 | 0.010 | 0.195 | 0.9969 | 3.1e-04 | 3.76 [3.76–4] |
| R_M6 | plain W1 | 1/+0.00 | 1 | 0.336 | 0.336 | 10.8 | 0.002 | 0.026 | 0.9720 | 4.6e-05 | 2.03 [1.79–56] |
| R_M6 | best GBP(≥10) | 0.9/+0.10 | 2 | 0.361 | 0.361 | 10.0 | 0.003 | 0.031 | 0.9686 | 8.7e-05 | 1.90 [1.58–70] |


## T6. Numerical convergence (same analysis code)

| config | f1/s | case vs base | bands ok | Δ gap_lo | Δ gap_hi | SM-even frac base → case | GBP(≥10) base → case | n_g at window center | median |v_g,MPB−FD|/v_g |
|---|---|---|---|---|---|---|---|---|---|
| H | 0.8/-0.10 | res48_tol1e-8 | True | -0.20% | +0.06% | 0.132 → 0.132 | 0.0483 → 0.0447 | 4.06 → 4.06 | 0.004 |
| H | 0.8/-0.10 | M+1 | True | -0.00% | -0.00% | 0.132 → 0.113 | 0.0483 → 0.0476 | 4.06 → 4.02 | 0.004 |
| H | 0.8/-0.10 | tol1e-9 | True | -0.00% | -0.00% | 0.132 → 0.132 | 0.0483 → 0.0483 | 4.06 → 4.06 | 0.004 |
| H | 1/+0.00 | res48_tol1e-8 | True | -0.20% | +0.06% | 0.106 → 0.107 | 0.0000 → 0.0000 | 3.92 → 3.92 | 0.003 |
| H | 1/+0.00 | M+1 | True | -0.00% | +0.00% | 0.106 → 0.106 | 0.0000 → 0.0000 | 3.92 → 3.92 | 0.003 |
| H | 1/+0.00 | tol1e-9 | True | -0.00% | -0.00% | 0.106 → 0.106 | 0.0000 → 0.0000 | 3.92 → 3.92 | 0.003 |
| H_r030 | 0.9/+0.10 | res48_tol1e-8 | True | -0.06% | -0.03% | 0.114 → 0.114 | 0.1945 → 0.1949 | 3.76 → 3.76 | 0.004 |
| H_r030 | 0.9/+0.10 | M+1 | True | +0.00% | +0.00% | 0.114 → 0.114 | 0.1945 → 0.1682 | 3.76 → 3.76 | 0.003 |
| H_r030 | 0.9/+0.10 | tol1e-9 | True | +0.00% | +0.00% | 0.114 → 0.114 | 0.1945 → 0.1945 | 3.76 → 3.76 | 0.003 |
| H_r030 | 1/+0.00 | res48_tol1e-8 | True | -0.06% | -0.03% | 0.107 → 0.108 | 0.0399 → 0.0398 | 3.87 → 3.87 | 0.003 |
| H_r030 | 1/+0.00 | M+1 | True | +0.00% | -0.00% | 0.107 → 0.109 | 0.0399 → 0.0393 | 3.87 → 3.87 | 0.003 |
| H_r030 | 1/+0.00 | tol1e-9 | True | -0.00% | -0.00% | 0.107 → 0.107 | 0.0399 → 0.0399 | 3.87 → 3.87 | 0.003 |
| R | 0.9/+0.10 | res48_tol1e-8 | True | -0.05% | +0.02% | 0.361 → 0.362 | 0.0313 → 0.0308 | 1.90 → 1.90 | 0.000 |
| R | 0.9/+0.10 | M+1 | True | +0.00% | +0.00% | 0.361 → 0.361 | 0.0313 → 0.0312 | 1.90 → 1.90 | 0.001 |
| R | 0.9/+0.10 | tol1e-9 | True | +0.00% | +0.00% | 0.361 → 0.361 | 0.0313 → 0.0313 | 1.90 → 1.90 | 0.000 |
| R | 1/+0.00 | res48_tol1e-8 | True | -0.05% | +0.02% | 0.337 → 0.337 | 0.0293 → 0.0199 | 2.03 → 2.03 | 0.000 |
| R | 1/+0.00 | M+1 | True | +0.00% | +0.00% | 0.337 → 0.336 | 0.0293 → 0.0262 | 2.03 → 2.03 | 0.000 |
| R | 1/+0.00 | tol1e-9 | True | -0.00% | +0.00% | 0.337 → 0.337 | 0.0293 → 0.0293 | 2.03 → 2.03 | 0.000 |


## Environment

```
{
 "python": "3.10.19 (main, Oct 21 2025, 16:43:05) [GCC 11.2.0]",
 "platform": "Linux-6.18.40.1-microsoft-standard-WSL2-x86_64-with-glibc2.35",
 "numpy": "1.24.3",
 "scipy": "1.11.1",
 "matplotlib": "3.10.6",
 "pandas": "2.0.3",
 "pytest": "7.4.0",
 "meep": "1.33.0-beta",
 "meep_with_mpi": true,
 "mpb_available": true,
 "git_head": "7783a200de84b17cd4c0913fc84f1a2f5033f73c"
}
```


Sweep provenance: `{"timestamp_utc": "2026-10-03T18:51:37+00:00", "git_head": "b06715e7cf8b227f79c62b50b8c50b8d38c1028b", "git_dirty": true, "meep": "1.33.0-beta", "numpy": "1.24.3", "scipy": "1.11.1"}`


## T7. Supercell edge-energy fraction over all solved designs (window max)

| config | M | min | max | designs with edge < 1e-3 | designs with window |
|---|---|---|---|---|---|
| R | 4 | 1.4e-03 | 3.2e-03 | 0 | 25 |
| R_M6 | 6 | 2.6e-05 | 1.1e-04 | 25 | 25 |
| H | 5 | 2.1e-05 | 6.2e-02 | 11 | 14 |
| H_r030 | 5 | 2.3e-04 | 1.8e-02 | 12 | 25 |
