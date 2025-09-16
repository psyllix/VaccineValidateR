# tests/testthat/test-population-metrics.R
test_that("summarize_population_metrics returns expected structure", {
  patients <- data.table::data.table(
    STUDY_ID = 1,
    DOB = as.Date("2020-01-01")
  )
  visits <- data.table::data.table(
    STUDY_ID = 1,
    VISIT_DATE = as.Date("2022-01-01")
  )
  
  res <- summarize_population_metrics(
    patients = patients,
    visits = visits,
    years = 2022,
    verbose = FALSE
  )
  
  expect_s3_class(res, "data.table")
  expect_true("VISIT_YEAR" %in% names(res))
})
