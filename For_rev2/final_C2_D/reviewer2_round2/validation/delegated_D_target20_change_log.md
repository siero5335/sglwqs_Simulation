# Battery D 20% prevalence recalibration

Updated: 2026-08-06 JST

## Reason

The previous factorial binomial generator used a fixed intercept of `-3.0`. Across the 180 unique binomial Battery D datasets this produced realized event prevalence of approximately 5.6%--8.2%, although rare-event behavior was not an explicit Battery D factor. The sparse-group-lasso binomial fits consequently had only about 33--49 events in each 60% training split and were substantially slower.

## Change

- Each factorial binomial dataset now solves a dataset-specific intercept so that mean model-implied event probability is exactly 0.20.
- The Bernoulli draw remains seed-deterministic.
- The target, expected prevalence, realized prevalence, positive count, and calibrated intercept are saved in each atomic result's `data_diagnostics`.
- The production manifest records `binomial_target_prevalence = 0.20` for binomial Battery D/E jobs.
- Smoke, preflight, and final production-gate checks verify the prevalence design.
- The balanced/unbalanced group definitions and all effect profiles were left unchanged in this recalibration.

## Archived previous run

The previous rare-event run was stopped before modification and moved intact to:

`<project-root>/reviewer2_round2/results/archive/rare_event_b0_minus3_20260806_0825`

It contains 59 completed successful atomic jobs and 6 in-progress SGL-WQS checkpoints. These files must not be mixed with the 20% prevalence run.

Battery C's 400 production outputs and passed gate were not modified.

## Verification

- Smoke: all checks passed, including all four wrappers and resume behavior.
- Preflight: 1440/1440 manifest jobs; 180/180 unique binomial datasets; expected prevalence exactly 0.20; realized prevalence range 0.171--0.232; identical data across methods.
