# GitHub reproduction runbook

1. Extract or clone the repository and use its root as the working directory.
2. Run `Rscript install_dependencies.R` and confirm that the bundled SGL-WQS source loads.
3. Run the three smoke scripts listed in `README.md`.
4. Run A/B/E, F–G, and H independently. These blocks can be scheduled separately, but do not launch a duplicate block against the same result directory.
5. For final D, run the target-0.20 wrapper in `final_C2_D/`. It performs prevalence preflight, production, gate validation, and summary generation.
6. For C2, first place the 400 binary reference outputs at the path documented in `README.md`, then run the C2 wrapper. The preflight verifies the exact matched exposure/covariate realization and train/validation membership.
7. Run the aggregation/table/figure scripts only after all required blocks are complete.

Recommended environment variables:

```sh
export R2R2_WORKERS=8
export R2R2_N_BOOT=200
export R2R2_FORCE=false
export R2R2_SAVE_RAW_FIT=false
```

Production is restartable: valid atomic files remain intact and are skipped. Before any restart, inspect the current process list and result counts to avoid duplicate worker pools.

The full completed output is not included here. Compare a reproduced run against `reference_metadata/FINAL_BATTERY_COMPLETION.csv`, `reference_metadata/FINAL_QA_CHECKS.csv`, and the SHA/source details in `reference_metadata/FINAL_SOURCE_PROVENANCE.csv`.

