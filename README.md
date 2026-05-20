# SGL-WQS Simulation Studies

This repository contains the reproducible simulation and benchmarking notebooks used to evaluate `sglwqs`, an R package for Weighted Quantile Sum (WQS) regression with Sparse Group Lasso.

The package code is maintained separately in the `sglwqs` repository. This repository is intended for public release of the simulation workflows, not for development of the package itself.

## Contents

| File | Purpose |
| --- | --- |
| `sglwqs_complete_data_100k.Rmd` | Large complete-data simulation with a binary outcome and known positive, negative, and null mixture effects. |
| `Comparison_with_Other_Methods__30-Seed_Replication.Rmd` | Benchmark using example datasets from `gWQS`, `qgcomp`, and `groupWQS`, replicated across 30 seeds. |
| `compare_methods_correlation_robustness_fixed_30seeds.Rmd` | Correlation robustness study comparing methods across exposure correlation levels. |
| `compare_mixture_methods_by_samplesize_30seeds.Rmd` | Sample-size sensitivity study across multiple values of `n`. |
| `compare_mixture_methods_43vars_30seeds.Rmd` | High-dimensional environmental-mixture simulation with 43 exposure variables. |
| `sglwqs_alpha_sensitivity_30seeds.Rmd` | Sensitivity analysis for the sparse group lasso mixing parameter `asparse`. |

## Requirements

Install the `sglwqs` package first. During local development, this can be done from a sibling checkout:

```r
install.packages("../sglwqs", repos = NULL, type = "source")
```

Install the simulation dependencies:

```r
install.packages(c(
  "MASS",
  "qgcomp",
  "gWQS",
  "groupWQS",
  "future",
  "ggplot2",
  "dplyr",
  "tidyr",
  "purrr",
  "knitr",
  "kableExtra",
  "rmarkdown",
  "tictoc",
  "tidyverse"
))
```

The `sglwqs` package also depends on packages such as `sparsegl`; those are installed with the package.

## Running the Analyses

Render any analysis from the repository root:

```r
rmarkdown::render("sglwqs_complete_data_100k.Rmd")
rmarkdown::render("compare_mixture_methods_by_samplesize_30seeds.Rmd")
```

Most notebooks run 30 independent seeds with 200 bootstrap replications per scenario. The full benchmark suite can take substantial time and memory, especially the sample-size and 43-variable studies. Several notebooks use `future::plan(multisession, workers = ...)`; reduce the worker count in the setup chunk if you are running on a constrained machine.

Generated HTML reports, knitr caches, and publication figures are intentionally ignored by Git. Figures are written to `figures/` when the notebooks are rendered.

## Reproducibility Notes

- Simulation data are generated inside the R Markdown files; no private raw data are required.
- Benchmark datasets are loaded from public R packages (`gWQS`, `qgcomp`, and `groupWQS`).
- Random seeds are declared in each notebook.
- Some numerical results may vary slightly across operating systems, R versions, package versions, and parallel backends.
- Session information is printed in the longer benchmark notebooks to help document the computational environment.

## License

This repository is released under the MIT License. See `LICENSE` for details.
