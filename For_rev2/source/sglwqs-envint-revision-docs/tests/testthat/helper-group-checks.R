check_canonical_groups <- function(actual, canonical) {
  if (is.list(actual) && !is.null(actual$groups)) {
    actual <- actual$groups
  }

  expect_false(is.null(actual))
  expect_identical(
    sglwqs:::.normalize_group_metadata(actual),
    sglwqs:::.normalize_group_metadata(canonical)
  )
}
