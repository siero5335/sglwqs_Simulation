# Validation-stage calibration analysis

This folder contains the reviewer-response simulation scaffold for validation-stage
conditional p-values and split stability.

The computational work is run through `targets`; the R Markdown files only read
cached target outputs and render tables/figures.
Each profile uses its own target store, e.g. `_targets_smoke` and
`_targets_production`, so smoke results cannot cause production runs to skip.

## Quick run

```sh
cd For_rev1/validation_calibration
Rscript run_pipeline.R --profile smoke --workers 8
Rscript render_reports.R --profile smoke
```

## Production run

```sh
cd For_rev1/validation_calibration
Rscript run_pipeline.R --profile production --workers 8
Rscript render_reports.R --profile production
```

The production profile is intentionally not launched by this scaffold. It uses
500 global-null replicates, 300 partial-null replicates, 200 Gaussian-null
replicates, and 100 repeated splits per fixed-data split-stability scenario.

## Package source

By default the pipeline uses the installed `sglwqs` package. To test a source
tree instead, set `SGLWQS_SOURCE`:

```sh
SGLWQS_SOURCE=/path/to/sglwqs Rscript run_pipeline.R --profile smoke --workers 8
```

The sparse group mixing parameter is explicitly passed as `asparse = 0.05`.

## Outputs

Generated artifacts are written under:

```text
output/
```

Important files include:

- `output/tables/validation_replicate_results_<profile>.csv`
- `output/tables/validation_calibration_summary_<profile>.csv`
- `output/tables/validation_split_stability_<profile>.csv`
- `output/tables/validation_runtime_projection_<profile>.csv`
- `output/html/*.html`
