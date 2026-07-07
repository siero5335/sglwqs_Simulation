# High-Dimensional Simulation Summary

Generated: 2026-06-29 18:08:34 JST

## Scope

- p-sensitivity: n=1000, groups=10, p=50/100/200, 30 seeds.
- group/signal sensitivity: p=100, n=1000, group count and signal strength perturbations, 70 scenario-seeds.
- correlation/signal sensitivity: p=100, n=1000, groups=10, within-group rho=0.20/0.45/0.80 and signal strength perturbations, 70 scenario-seeds.

## Metric Notes

- Active direction accuracy is direction-only among truly active components. It does not require component-level GLM significance.
- SGL-WQS TNR is based on variable-level max positive/negative bootstrap selection frequency below a threshold among null variables.
- qgcomp and gWQS produce non-sparse weights for all components; sparse TNR is therefore reported only for SGL-WQS.

## Main Interpretation

- All high-dimensional simulations completed without method failures.
- SGL-WQS remained computationally stable through p=200 and across group/signal/correlation perturbations.
- Active direction recovery was generally above chance, but SGL-WQS did not show adjusted statistically significant superiority over comparators in the p-sensitivity analysis.
- Active-variable weight concentration decreased as dimensionality increased.
- SGL-WQS bootstrap selection frequencies separated active and null components weakly in high-dimensional settings; null components often had high selection frequencies.
- Overall, these results support feasibility and robustness against outright failure rather than strong high-dimensional component-selection performance.

## p-Sensitivity: Active Direction Accuracy



|   p| SGL-WQS|  gWQS| groupWQS| qgcomp|
|---:|-------:|-----:|--------:|------:|
|  50|   0.663| 0.685|    0.715|  0.737|
| 100|   0.693| 0.722|    0.715|  0.737|
| 200|   0.678| 0.619|    0.659|  0.674|

## p-Sensitivity: Active Weight Share



|   p| SGL-WQS|  gWQS| groupWQS| qgcomp|
|---:|-------:|-----:|--------:|------:|
|  50|   0.180| 0.204|    0.194|  0.209|
| 100|   0.095| 0.101|    0.105|  0.108|
| 200|   0.049| 0.046|    0.058|  0.051|

## p-Sensitivity: Median Runtime Seconds



|   p|   SGL-WQS|      gWQS| groupWQS| qgcomp|
|---:|---------:|---------:|--------:|------:|
|  50|   935.273|  1113.880|  820.859|  0.038|
| 100|  2854.134|  4596.478| 2249.407|  0.064|
| 200| 11295.162| 24775.106| 7066.413|  0.155|

## p-Sensitivity: SGL-WQS TPR/TNR by Selection-Frequency Threshold



|   p| Threshold| Seeds| Active_Rows| Null_Rows|   TPR|   TNR|   FPR| Active_Mean_MaxSelFreq| Null_Mean_MaxSelFreq|
|---:|---------:|-----:|-----------:|---------:|-----:|-----:|-----:|----------------------:|--------------------:|
|  50|       0.5|    30|         270|      1230| 0.841| 0.208| 0.792|                  0.684|                0.653|
|  50|       0.8|    30|         270|      1230| 0.293| 0.789| 0.211|                  0.684|                0.653|
| 100|       0.5|    30|         270|      2730| 0.944| 0.078| 0.922|                  0.750|                0.718|
| 100|       0.8|    30|         270|      2730| 0.448| 0.658| 0.342|                  0.750|                0.718|
| 200|       0.5|    30|         270|      5730| 0.956| 0.075| 0.925|                  0.755|                0.727|
| 200|       0.8|    30|         270|      5730| 0.463| 0.632| 0.368|                  0.755|                0.727|

## Group/Signal Sensitivity: Active Direction Accuracy



|Scenario                  | n_groups|truth_profile | SGL-WQS|  gWQS| groupWQS| qgcomp|
|:-------------------------|--------:|:-------------|-------:|-----:|--------:|------:|
|p100_g05_baseline3_main10 |        5|baseline3     |   0.511| 0.578|    0.778|  0.667|
|p100_g05_strong3_orth5    |        5|strong3       |   0.578| 0.644|    0.778|  0.778|
|p100_g05_weak3_orth5      |        5|weak3         |   0.533| 0.467|    0.733|  0.600|
|p100_g10_baseline3_main10 |       10|baseline3     |   0.678| 0.678|    0.678|  0.767|
|p100_g10_strong3_main10   |       10|strong3       |   0.700| 0.700|    0.700|  0.811|
|p100_g10_weak3_main10     |       10|weak3         |   0.589| 0.622|    0.678|  0.689|
|p100_g20_baseline3_main10 |       20|baseline3     |   0.644| 0.600|    0.700|  0.678|
|p100_g20_strong3_orth5    |       20|strong3       |   0.667| 0.600|    0.733|  0.756|
|p100_g20_weak3_orth5      |       20|weak3         |   0.422| 0.489|    0.556|  0.511|

## Correlation/Signal Sensitivity: Active Direction Accuracy



|Scenario                   | within_rho|truth_profile | SGL-WQS|  gWQS| groupWQS| qgcomp|
|:--------------------------|----------:|:-------------|-------:|-----:|--------:|------:|
|p100_r020_baseline3_main10 |       0.20|baseline3     |   0.689| 0.733|    0.678|  0.811|
|p100_r020_strong3_orth5    |       0.20|strong3       |   0.733| 0.778|    0.689|  0.822|
|p100_r020_weak3_orth5      |       0.20|weak3         |   0.533| 0.622|    0.622|  0.733|
|p100_r045_baseline3_main10 |       0.45|baseline3     |   0.678| 0.678|    0.678|  0.767|
|p100_r045_strong3_main10   |       0.45|strong3       |   0.700| 0.700|    0.700|  0.811|
|p100_r045_weak3_main10     |       0.45|weak3         |   0.589| 0.622|    0.678|  0.689|
|p100_r080_baseline3_main10 |       0.80|baseline3     |   0.600| 0.578|    0.678|  0.667|
|p100_r080_strong3_orth5    |       0.80|strong3       |   0.578| 0.578|    0.778|  0.733|
|p100_r080_weak3_orth5      |       0.80|weak3         |   0.511| 0.489|    0.689|  0.622|

## Primary Output Files

- `p_sensitivity/`: p=50/100/200 summaries.
- `group_signal_sensitivity/`: group count and signal strength summaries.
- `correlation_signal_sensitivity/`: within-group correlation and signal strength summaries.
- `summary_table_*.csv`: compact cross-method tables used in this report.
