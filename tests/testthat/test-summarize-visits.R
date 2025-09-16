# tests/testthat/test-summarize-visits.R
test_that("summarize_visit_evaluations runs in both long and wide formats", {
  dt <- data.table::data.table(
    STUDY_ID = 1,
    VISIT_DATE = as.Date(c("2020-03-01", "2020-06-01")),
    ANTIGEN = c("HEPB", "HEPB"),
    DUE = c(TRUE, TRUE),
    GIVEN = c(TRUE, FALSE),
    MISSED = c(FALSE, TRUE)
  )
  
  res_long <- summarize_visit_evaluations(dt, output_format = "long", verbose = FALSE)
  res_wide <- summarize_visit_evaluations(dt, output_format = "wide", verbose = FALSE)
  
  expect_s3_class(res_long, "data.table")
  expect_s3_class(res_wide, "data.table")
  expect_true("HEPB" %in% names(res_wide))
})
