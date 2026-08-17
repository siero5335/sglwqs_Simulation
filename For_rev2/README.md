# SGL-WQS Reviewer #2 Round 2 simulation code

This archive contains the code, manifests, and exact local source snapshots used for the Reviewer #2 Round 2 simulations. It is prepared for GitHub upload: production results, atomic `.rds` files, terminal logs, temporary files, and machine-specific absolute paths are excluded.

## Final analysis set

- Local Batteries A–H: 5,320/5,320 atomic jobs, zero failures.
- Final paired outcome-family comparison: C2 Gaussian 400 jobs plus 400 matched binary reference jobs.
- Final Battery D: 1,440 jobs with expected binomial prevalence calibrated to 0.20.
- Final adopted set after excluding the overlapping original C audit run: 5,720 jobs.

The original C and the early rare-event D code are retained only under `legacy_original_C_D/`. They are not the final manuscript evidence. Use `final_C2_D/` for the adopted C2 and D analyses.

## Directory layout

- `reviewer2_round2/R/`: common generators, wrappers, metrics, and I/O helpers for A, B, E, F, G, and H
- `reviewer2_round2/scripts/`: smoke, production, aggregation, table, and figure scripts
- `reviewer2_round2/config/`: production and smoke manifests
- `final_C2_D/reviewer2_round2/`: exact final C2 paired-family and target-prevalence-0.20 D code, manifests, and validation definitions
- `legacy_original_C_D/`: non-final original C/D scripts retained for auditability
- `source/sglwqs-envint-revision-docs/`: exact development source snapshot used by the final simulations
- `source/sglwqs_old_analysis/sglwqs-main/`: legacy source snapshot used for A1–A2
- `reference_metadata/`: final completion, QA, and source-provenance records
- `CODE_INVENTORY.csv`: analysis block-to-entry-point map
- `FILE_SHA256SUMS.csv`: SHA-256 manifest

## Software provenance

Most final simulations used SGL-WQS version `0.8.13.9001`, commit `2fdd519e520a7dad1162810643e175cd616b1154`, with source SHA-256 `fb6ed94e2a26bc71894879d4323d6e7e87f1a27b28aee31d825aa77e36adb6a0`. A1–A2 used the archived `0.8.13` analysis source. See `reference_metadata/FINAL_SOURCE_PROVENANCE.csv` for the battery-level mapping.

The scripts load the bundled snapshot with `pkgload::load_all()`. The corresponding remote branch is:

```r
remotes::install_github("siero5335/sglwqs@paper/envint-revision-docs")
```

The bundled snapshot is preferred for exact reproduction because a moving branch can change after publication.

## Dependencies

Use a recent R installation and run:

```sh
Rscript install_dependencies.R
```

The main dependencies are MASS, Matrix, qgcomp, gWQS, groupWQS, pkgload, dplyr, tidyr, purrr, readr, ggplot2, knitr, future.apply, and digest. Package/session records from the completed run are represented in the accompanying full-results archive; this code package intentionally excludes execution logs.

## Smoke tests

From the archive root:

```sh
Rscript reviewer2_round2/scripts/01_smoke_test.R
Rscript reviewer2_round2/scripts/02_smoke_gaussian_primary_benchmarks.R
Rscript reviewer2_round2/scripts/03_smoke_dense_group_signal.R
```

## Production entry points

Core A, B, and E with eight outer workers:

```sh
R2R2_WORKERS=8 bash reviewer2_round2/scripts/run_core_A_B_E.sh
```

Gaussian correlation and active-component benchmarks F–G:

```sh
R2R2_WORKERS=8 bash reviewer2_round2/scripts/run_gaussian_primary_benchmarks.sh
```

Dense group signal H:

```sh
R2R2_WORKERS=8 bash reviewer2_round2/scripts/run_dense_group_signal.sh
```

Final target-prevalence-0.20 D:

```sh
R2R2_WORKERS=8 bash final_C2_D/reviewer2_round2/scripts/run_delegated_D_target20_production.sh
```

Final exactly paired C2:

```sh
R2R2_WORKERS=1 bash final_C2_D/reviewer2_round2/scripts/run_delegated_C2_paired_production.sh
```

C2 requires the 400 completed binary reference `.rds` files under `final_C2_D/analysis/validation_calibration/output/manifests/production/split_stability/`. Those reference outputs are data, not code, and are therefore intentionally absent from this GitHub code archive. They are present in the complete results/reproducibility materials. The preflight script stops if the pairing inputs, hashes, or split membership are not correct.

## Reproducibility notes

- Job-specific data, split, method, and bootstrap seeds are stored in the manifests or generated deterministically by the helper functions.
- Parallelization is across atomic jobs. Each job constrains BLAS/OpenMP threads to one to avoid nested oversubscription.
- `R2R2_WORKERS`, `R2R2_N_BOOT`, `R2R2_FORCE`, and `RSCRIPT` can be supplied as environment variables.
- Existing valid atomic outputs are reused unless `R2R2_FORCE=true` is explicitly set.
- Validation-stage p-values are conditional on direction retention after the split-dependent screening stage. They are exploratory summaries, not unconditional confirmatory inference.
- Cross-method component attribution is not a common estimand in every method and should not be interpreted as a universal superiority test.

