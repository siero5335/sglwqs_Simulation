test_that("plot_validation_results() works for sglwqs_mids pooled validation results", {
  fit_mice <- structure(
    list(
      m = 5,
      n_successful = 5,
      validation = TRUE,
      pooled = list(
        validation = list(
          has_group_results = FALSE,
          group_results = NULL,
          wqs_pos = list(estimate = 0.42, se = 0.10, p_value = 0.01, fmi = 0.12, df = 20),
          wqs_neg = list(estimate = -0.25, se = 0.08, p_value = 0.03, fmi = 0.18, df = 18)
        )
      )
    ),
    class = "sglwqs_mids"
  )

  p <- plot_validation_results(fit_mice, conf_level = 0.95)

  expect_s3_class(p, "ggplot")
  expect_equal(as.character(p$data$group), c("Overall", "Overall"))
  expect_equal(as.character(p$data$direction), c("Positive", "Negative"))
  expect_match(p$labels$title, "Rubin's Rules")
})


test_that("plot_validation_results() errors for sglwqs_mids without pooled validation", {
  fit_mice <- structure(
    list(
      m = 5,
      n_successful = 5,
      validation = FALSE,
      pooled = list(validation = NULL)
    ),
    class = "sglwqs_mids"
  )

  expect_error(
    plot_validation_results(fit_mice),
    "Pooled validation information not available"
  )
})
