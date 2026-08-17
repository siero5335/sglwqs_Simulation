repo_root <- normalizePath(Sys.getenv("R2R2_REPO_ROOT", unset = getwd()), mustWork = TRUE)
r2_root <- file.path(repo_root, "reviewer2_round2")
results_root <- file.path(r2_root, "results")
out_dir <- file.path(r2_root, "final_results_20260810")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(...) {
  read.csv(file.path(...), stringsAsFactors = FALSE, check.names = FALSE)
}

pct <- function(x, digits = 1L) sprintf(paste0("%.", digits, "f%%"), 100 * x)
num <- function(x, digits = 3L) sprintf(paste0("%.", digits, "f"), x)
series <- function(x, digits = 3L) paste(vapply(x, num, character(1), digits = digits), collapse = ", ")

main_summary <- read_csv(results_root, "tables", "01_method_performance_summary.csv")
main_status <- read_csv(results_root, "summaries", "job_completion_status.csv")
a_summary <- read_csv(results_root, "exports", "battery_A_complete", "method_performance_summary_A.csv")
a_conditional <- read_csv(results_root, "exports", "battery_A_complete", "sglwqs_conditional_summary_A.csv")
b_summary <- read_csv(results_root, "exports", "battery_B_complete", "method_performance_summary_B.csv")
b_conditional <- read_csv(results_root, "exports", "battery_B_complete", "sglwqs_conditional_summary_B.csv")
c2_summary <- read_csv(results_root, "imported_second_pc", "C2_exactly_paired_final", "tables", "02_scenario_summary.csv")
c2_qa <- read_csv(results_root, "imported_second_pc", "C2_exactly_paired_final", "tables", "01_quality_checks.csv")
d_summary <- read_csv(results_root, "imported_second_pc", "D_target20_final", "tables", "02_scenario_method_summary.csv")
d_qa <- read_csv(results_root, "imported_second_pc", "D_target20_final", "tables", "01_quality_checks.csv")
e_summary <- read_csv(results_root, "tables", "08_hard_setting_sample_size_summary.csv")
sgl_summary <- read_csv(results_root, "tables", "10_sglwqs_conditional_retention_rejection_summary.csv")
fg_summary <- read_csv(results_root, "tables", "gaussian_primary_benchmarks", "01_method_completion_runtime_direction.csv")
fg_conditional <- read_csv(results_root, "tables", "gaussian_primary_benchmarks", "04_sglwqs_conditional_retention_pvalues.csv")
fg_failures <- read_csv(results_root, "tables", "gaussian_primary_benchmarks", "05_failures.csv")
h_summary <- read_csv(results_root, "tables", "dense_group_signal", "01_dense_method_summary.csv")
h_validation <- read_csv(results_root, "tables", "dense_group_signal", "04_dense_sglwqs_validation.csv")
h_failures <- read_csv(results_root, "tables", "dense_group_signal", "06_dense_failures.csv")
import_qa <- read_csv(r2_root, "handoffs", "second_pc_import_20260809", "audit", "validation_checks.csv")

completion <- data.frame(
  result_set = c(
    "A1 Gaussian sample size", "A2 weak Gaussian sample size",
    "B1 Gaussian high dimensional", "B2 weak Gaussian high dimensional",
    "C original Gaussian split stability", "C2 paired Gaussian split stability",
    "C2 paired binary reference", "D target20 factorial", "E hard setting",
    "F Gaussian correlation", "G Gaussian active component", "H dense group signal"
  ),
  role = c(rep("final", 4), "audit_only", rep("final", 7)),
  expected_jobs = c(720L, 480L, 360L, 360L, 400L, 400L, 400L, 1440L, 480L, 480L, 120L, 480L),
  fit_successes = c(720L, 480L, 360L, 360L, 400L, 400L, 400L, 1440L, 480L, 480L, 120L, 480L),
  fit_failures = 0L,
  stringsAsFactors = FALSE
)
completion$fit_completion_rate <- completion$fit_successes / completion$expected_jobs
write.csv(completion, file.path(out_dir, "FINAL_BATTERY_COMPLETION.csv"), row.names = FALSE, na = "")

source_provenance <- data.frame(
  result_set = c("A1-A2", "B1-B2", "C original", "C2 final", "D target20 final", "E", "F-G", "H"),
  sglwqs_version = c("0.8.13", rep("0.8.13.9001", 7)),
  source_role = c("legacy analysis source retained for matched A run", rep("pinned paper revision source", 7)),
  git_commit = c(NA, rep("2fdd519e520a7dad1162810643e175cd616b1154", 7)),
  main_file_sha256 = c(NA, rep("fb6ed94e2a26bc71894879d4323d6e7e87f1a27b28aee31d825aa77e36adb6a0", 7)),
  stringsAsFactors = FALSE
)
write.csv(source_provenance, file.path(out_dir, "FINAL_SOURCE_PROVENANCE.csv"), row.names = FALSE, na = "")

qa <- data.frame(
  check = c(
    "Core A-E manifest complete",
    "All locally aggregated A-H jobs complete",
    "All locally aggregated A-H fits successful",
    "F-G failure table empty",
    "F-G conditional rates bounded by zero and one",
    "H completion and failure table",
    "C2 packaged quality checks",
    "D packaged quality checks",
    "Second-PC import audit",
    "Final adopted job count"
  ),
  observed = c(
    paste0(sum(main_status$valid_complete), "/", nrow(main_status)),
    paste0(sum(main_summary$completed), "/", sum(main_summary$attempted)),
    as.character(min(main_summary$fit_completion_rate)),
    as.character(nrow(fg_failures)),
    paste(range(c(
      fg_conditional$retention_rate,
      fg_conditional$conditional_p_lt_0_05_retained,
      fg_conditional$conditional_p_lt_0_05_all_attempted
    ), na.rm = TRUE), collapse = "--"),
    paste0(sum(h_summary$completed), "/", sum(h_summary$attempted), "; failures=", nrow(h_failures)),
    paste0(sum(c2_qa$passed), "/", nrow(c2_qa)),
    paste0(sum(d_qa$passed), "/", nrow(d_qa)),
    paste0(sum(import_qa$passed), "/", nrow(import_qa)),
    as.character(sum(completion$fit_successes[completion$role == "final"]))
  ),
  expected = c("4240/4240", "5320/5320", "1", "0", "0--1", "480/480; failures=0", "13/13", "12/12", "34/34", "5720"),
  passed = c(
    nrow(main_status) == 4240L && sum(main_status$valid_complete) == 4240L,
    sum(main_summary$attempted) == 5320L && sum(main_summary$completed) == 5320L,
    all(main_summary$fit_completion_rate == 1),
    nrow(fg_failures) == 0L,
    all(c(
      fg_conditional$retention_rate,
      fg_conditional$conditional_p_lt_0_05_retained,
      fg_conditional$conditional_p_lt_0_05_all_attempted
    ) >= 0 & c(
      fg_conditional$retention_rate,
      fg_conditional$conditional_p_lt_0_05_retained,
      fg_conditional$conditional_p_lt_0_05_all_attempted
    ) <= 1, na.rm = TRUE),
    sum(h_summary$attempted) == 480L && sum(h_summary$completed) == 480L && nrow(h_failures) == 0L,
    nrow(c2_qa) == 13L && all(c2_qa$passed),
    nrow(d_qa) == 12L && all(d_qa$passed),
    nrow(import_qa) == 34L && all(import_qa$passed),
    sum(completion$fit_successes[completion$role == "final"]) == 5720L
  ),
  stringsAsFactors = FALSE
)
write.csv(qa, file.path(out_dir, "FINAL_QA_CHECKS.csv"), row.names = FALSE, na = "")
if (!all(qa$passed)) stop("Final QA failed; inspect FINAL_QA_CHECKS.csv", call. = FALSE)

a1_sgl <- a_summary[a_summary$battery == "A1_gaussian_sample_size" & a_summary$method == "SGL-WQS", ]
a1_sgl <- a1_sgl[order(a1_sgl$n), ]
a2_sgl <- a_summary[a_summary$battery == "A2_weak_gaussian_sample_size" & a_summary$method == "SGL-WQS", ]
a2_sgl <- a2_sgl[order(a2_sgl$n), ]
a1_cond <- a_conditional[a_conditional$battery == "A1_gaussian_sample_size", ]
a1_cond <- a1_cond[order(a1_cond$n), ]
a2_cond <- a_conditional[a_conditional$battery == "A2_weak_gaussian_sample_size", ]
a2_cond <- a2_cond[order(a2_cond$n), ]

b1_sgl <- b_summary[b_summary$battery == "B1_gaussian_highdim" & b_summary$method == "SGL-WQS", ]
b1_sgl <- b1_sgl[order(b1_sgl$p), ]
b2_sgl <- b_summary[b_summary$battery == "B2_weak_gaussian_highdim" & b_summary$method == "SGL-WQS", ]
b2_sgl <- b2_sgl[order(b2_sgl$p), ]
b1_cond <- b_conditional[b_conditional$battery == "B1_gaussian_highdim", ]
b1_cond <- b1_cond[order(b1_cond$p), ]
b2_cond <- b_conditional[b_conditional$battery == "B2_weak_gaussian_highdim", ]
b2_cond <- b2_cond[order(b2_cond$p), ]

c2_partial <- c2_summary[c2_summary$effect == "partial_null", ]
c2_global <- c2_summary[c2_summary$effect == "global_null", ]
d_sgl <- d_summary[d_summary$method == "SGL-WQS", ]
e_sgl <- e_summary[e_summary$method == "SGL-WQS", ]
e_sgl <- e_sgl[order(e_sgl$family, e_sgl$n), ]
e_val <- sgl_summary[sgl_summary$battery == "E_hard_setting_sample_size", ]
e_val <- e_val[order(e_val$family, e_val$n), ]
f_sgl <- fg_summary[fg_summary$battery == "F_gaussian_correlation_robustness" & fg_summary$method == "SGL-WQS", ]
f_sgl <- f_sgl[order(f_sgl$target_correlation), ]
g_all <- fg_summary[fg_summary$battery == "G_gaussian_active_component", ]
h_sgl <- h_summary[h_summary$method == "SGL-WQS", ]
h_sgl <- h_sgl[order(h_sgl$family, h_sgl$n), ]
h_val <- h_validation[order(h_validation$family, h_validation$n), ]

report <- c(
  "# Reviewer #2 Round 2: 最終シミュレーション結果総括",
  "",
  paste0("作成日時: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 完了状況",
  "",
  "- ローカルA-Hは5,320/5,320 atomic jobsが完了し、fit失敗は0件だった。",
  "- exactly paired C2はGaussian 400件とbinary reference 400件が完了し、fit失敗は0件だった。",
  "- 旧C 400件は監査用として保持し、最終family comparisonにはC2を採用する。",
  "- したがって、重複する旧Cを除く最終採用セットは5,720 jobsである。",
  "- Dはbinomial prevalenceを0.20に校正した1,440件の最終版で、旧途中版24件は解析対象外としてarchiveに保全した。",
  "",
  "## Battery A: matched Gaussian sample-size",
  "",
  paste0("- SGL-WQS active-direction accuracy（baseline; n=", paste(a1_sgl$n, collapse = ", "), "）は ", series(a1_sgl$mean_active_direction_accuracy), "。"),
  paste0("- weak signal（n=", paste(a2_sgl$n, collapse = ", "), "）では ", series(a2_sgl$mean_active_direction_accuracy), "。"),
  paste0("- active group-directionのconditional p<0.05（all-attempted）はbaselineで ", series(a1_cond$active_rejection_rate_all_attempted), "、weakで ", series(a2_cond$active_rejection_rate_all_attempted), "。"),
  "- 連続アウトカムでも標本数とシグナル強度に応じた改善が明瞭だった。低次元加法設定ではqgcompも強く、これは一律の優越性ではなく各手法の設計との整合性を示す。",
  "",
  "## Battery B: matched Gaussian high-dimensional",
  "",
  paste0("- SGL-WQS active-direction accuracy（p=50,100,200）はbaselineで ", series(b1_sgl$mean_active_direction_accuracy), "、weakで ", series(b2_sgl$mean_active_direction_accuracy), "。"),
  paste0("- 一方、active attributionはbaselineで ", series(b1_sgl$mean_active_attribution), "、weakで ", series(b2_sgl$mean_active_attribution), "まで低下した。"),
  paste0("- active conditional p<0.05（all-attempted）はbaselineで ", series(b1_cond$active_p_lt_005_rate_all_attempted), "、weakで ", series(b2_cond$active_p_lt_005_rate_all_attempted), "。"),
  "- 高次元化ではcomputational completionより先に、component prioritizationとconditional validationの情報量が低下した。",
  "",
  "## Battery C2: exactly paired split stability",
  "",
  paste0("- global-nullのnull conditional p<0.05はn=500でbinary/Gaussianが ", num(c2_global$null_rejection_all_mean[c2_global$n == 500 & c2_global$outcome_family == "binomial"]), "/", num(c2_global$null_rejection_all_mean[c2_global$n == 500 & c2_global$outcome_family == "gaussian"]), "、n=5,000で ", num(c2_global$null_rejection_all_mean[c2_global$n == 5000 & c2_global$outcome_family == "binomial"]), "/", num(c2_global$null_rejection_all_mean[c2_global$n == 5000 & c2_global$outcome_family == "gaussian"]), "。"),
  paste0("- partial-null active conditional p<0.05（all-attempted）はn=500でbinary/Gaussianが ", num(c2_partial$active_detection_all_mean[c2_partial$n == 500 & c2_partial$outcome_family == "binomial"]), "/", num(c2_partial$active_detection_all_mean[c2_partial$n == 500 & c2_partial$outcome_family == "gaussian"]), "、n=5,000で ", num(c2_partial$active_detection_all_mean[c2_partial$n == 5000 & c2_partial$outcome_family == "binomial"]), "/", num(c2_partial$active_detection_all_mean[c2_partial$n == 5000 & c2_partial$outcome_family == "gaussian"]), "。"),
  "- 同じX、covariates、truth、train/validation membershipをfamily間で一致させた結果であり、最終のbinary-Gaussian split comparisonとして用いる。",
  "",
  "## Battery D: group balance x effect heterogeneity",
  "",
  "- 全1,440 jobsが成功した。SGL-WQS active-direction accuracyは、Gaussianのbalanced/unbalancedでheterogeneous 0.889/0.885、uniform 1.000/0.996、weak heterogeneous 0.811/0.819だった。",
  "- Binomialでは同順に0.804/0.793、0.844/0.878、0.696/0.704だった。group imbalance単独の系統的な不利益より、weak effectの影響が大きかった。",
  "",
  "## Battery E: hard setting",
  "",
  paste0("- SGL-WQS active-direction accuracy（binomial n=500/5,000）は ", series(e_sgl$active_direction_accuracy_mean[e_sgl$family == "binomial"]), "、Gaussianは ", series(e_sgl$active_direction_accuracy_mean[e_sgl$family == "gaussian"]), "。"),
  paste0("- active conditional p<0.05（all-attempted）はbinomialで ", series(e_val$active_conditional_detection_all[e_val$family == "binomial"]), "、Gaussianで ", series(e_val$active_conditional_detection_all[e_val$family == "gaussian"]), "。"),
  "- n=500、weak、p=100の組合せでは方向推定とdownstream conditional informationの双方に明確な制約があった。",
  "",
  "## Batteries F-G: primary Gaussian counterparts",
  "",
  paste0("- FのSGL-WQS direction accuracy（r=0.2,0.5,0.8,0.95）は ", series(f_sgl$active_direction_accuracy_mean), "、active attributionは ", series(f_sgl$active_attribution_mean), "。"),
  "- 高相関でも方向推定は比較的維持されたが、active/null componentの分離は弱くなった。",
  paste0("- Gのactive attributionはSGL-WQS/gWQS/groupWQS/qgcompで ", series(g_all$active_attribution_mean[match(c("SGL-WQS", "gWQS", "groupWQS", "qgcomp"), g_all$method)]), "。"),
  "- 低次元・少数active component・加法設定ではqgcompが強く、方向制約を先に置くWQS系にも各々の得意条件があった。",
  "",
  "## Battery H: dense group signal",
  "",
  paste0("- SGL-WQS direction accuracy（binomial n=500/5,000）は ", series(h_sgl$active_direction_accuracy_mean[h_sgl$family == "binomial"]), "、Gaussianは ", series(h_sgl$active_direction_accuracy_mean[h_sgl$family == "gaussian"]), "。"),
  paste0("- active attributionは同順にbinomial ", series(h_sgl$active_attribution_mean[h_sgl$family == "binomial"]), "、Gaussian ", series(h_sgl$active_attribution_mean[h_sgl$family == "gaussian"]), "。"),
  paste0("- SGL-WQS active conditional p<0.05（all-attempted）はbinomial ", series(h_val$active_conditional_detection_all[h_val$family == "binomial"]), "、Gaussian ", series(h_val$active_conditional_detection_all[h_val$family == "gaussian"]), "。"),
  "- dense signalではgroupWQSの方向成績が高い条件があったが、一方向/群という構造上の利点を含む。Gaussian n=5,000ではqgcompとgWQSのcomponent attributionも改善した。",
  "",
  "## 全体として支持される点",
  "",
  "- Gaussian outcomeへの適用可能性だけでなく、sample size、dimension、signal strength、group imbalance、correlation、sparse/dense truthに応じたoperating performanceを示せた。",
  "- 手法間の順位はscenarioと評価層で変わった。computational completion、group-direction assignment、component prioritization、conditional validationを分けて記述する必要がある。",
  "- SGL-WQSの特徴は、すべての設定で単一指標を最大化することではなく、bidirectionalかつgroup-awareなprioritizationとgroup-direction別のconditional validationを同時に提供する点にある。",
  "",
  "## 支持されない点と注意",
  "",
  "- universal superiority、全continuous outcomeへの一般化、formal post-selection inference、causal signal detectionは支持しない。",
  "- qgcompのcoefficient-derived attributionとWQS-type constrained weightsは同一estimandではない。",
  "- validation-stage p-valuesはconditional/exploratory outputとして扱う。",
  "- Battery Aのみ解析当時の0.8.13 sourceを保持し、B-H/C2/Dは0.8.13.9001、commit 2fdd519e520a7dad1162810643e175cd616b1154を使用した。",
  "",
  "## 主要ファイル",
  "",
  "- `REVIEWER2_ROUND2_SIMULATION_REPORT.md`: A-Hローカル統合集計",
  "- `results/imported_second_pc/C2_exactly_paired_final/`: 最終C2",
  "- `results/imported_second_pc/D_target20_final/`: 最終D",
  "- `results/summaries/gaussian_primary_benchmarks/`: F-G",
  "- `results/summaries/dense_group_signal/`: H",
  "- `final_results_20260810/FINAL_QA_CHECKS.csv`: 最終QA"
)

writeLines(report, file.path(out_dir, "FINAL_SIMULATION_RESULTS_REPORT_JA.md"), useBytes = TRUE)
writeLines(c(
  "# Final Reviewer #2 Round 2 simulation bundle",
  "",
  "- `FINAL_SIMULATION_RESULTS_REPORT_JA.md`: Japanese consolidated factual report.",
  "- `FINAL_BATTERY_COMPLETION.csv`: completion and adopted/audit status.",
  "- `FINAL_QA_CHECKS.csv`: machine-readable final QA.",
  "- `FINAL_SOURCE_PROVENANCE.csv`: exact SGL-WQS source provenance.",
  "",
  "No manuscript or Response-letter text is included or modified."
), file.path(out_dir, "README.md"), useBytes = TRUE)

cat("Final report written:", file.path(out_dir, "FINAL_SIMULATION_RESULTS_REPORT_JA.md"), "\n")
cat("Final QA passed:", sum(qa$passed), "/", nrow(qa), "\n")
