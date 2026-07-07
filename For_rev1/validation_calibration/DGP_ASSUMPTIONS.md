# DGP assumptions and provenance

The 13-exposure data-generating mechanism follows the existing paper simulation
Rmds:

- `../original_simulations/compare_mixture_methods_by_samplesize_30seeds.Rmd`
- `../original_simulations/compare_methods_correlation_robustness_fixed_30seeds.Rmd`
- the submitted minor-direction threshold description and the sensitivity Rmd in
  `../additional_sensitivity/minor_direction_threshold_sensitivity.Rmd`

Exposure groups:

- PCBs: `pcb_118`, `pcb_138`, `pcb_153`, `pcb_180`, `pcb_187`
- Metals: `cd`, `pb`, `hg`, `se`, `mg`
- PFASs: `PFOA`, `PFNA`, `PFOS`

Covariates are generated as in the existing simulation:

- `sex ~ Bernoulli(0.5)`
- `age` is a truncated rounded normal centered at 50 years
- `bmi` depends on age, sex, and normal noise, truncated to 15-45

The binary outcome uses the existing logistic scale with intercept `-3.0` and
covariate coefficients `sex * 0.30`, `(age - 50) / 10 * 0.02`, and
`(bmi - 25) / 5 * 0.05`.

Scenarios:

- Global-null binary scenarios set all exposure coefficients to zero.
- Partial-null binary scenarios use the existing 13-variable effects:
  - PCBs positive: `0.04, 0.15, 0.20, 0.08, 0.05`
  - Metals mixed: `-0.12, -0.18, 0.10, 0.14, -0.02`
  - PFASs null: `0, 0, 0`
- Gaussian sensitivity scenarios use a continuous outcome with all exposure
  coefficients set to zero.

Correlation handling:

- For global-null binary and Gaussian-null scenarios, within-group correlation is
  set to `0.20` or `0.95` as requested.
- Between-group correlation is set to `0` for global-null and Gaussian-null
  scenarios.
- For partial-null binary scenarios, the prompt varies between-group correlation
  (`0` or `0.30`) but does not specify a within-group value. This scaffold fixes
  the partial-null within-group correlation at `0.60`, matching the order of the
  existing paper DGP. If production should instead cross partial-null scenarios
  with `within_r = 0.20/0.95`, update `make_validation_scenarios()` in
  `R/02_dgp.R` before launching production.

Estimand note:

SGL-WQS validation-stage coefficients are downstream coefficients for
training-estimated WQS indices. They are not the same as the individual exposure
coefficients used in the DGP. Coverage is therefore summarized for true-null
group-directions where the null value is exactly zero. Active-direction coverage
is left as `NA` unless a separate, index-specific estimand is defined.
