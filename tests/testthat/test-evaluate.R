# tests/testthat/test-evaluate.R
test_that("evaluate_visits returns expected columns", {
  visits <- data.table::data.table(
    STUDY_ID = c(1, 1),
    VISIT_DATE = as.Date(c("2020-03-01", "2020-06-01")),
    DOB = as.Date("2020-01-01")
  )
  
  # fake antigen table to pass in
  antigens <- data.table::data.table(
    STUDY_ID = 1,
    ANTIGEN = "HEPB",
    DATE_GIVEN = as.Date("2020-01-01"),
    DELAYED = FALSE,
    DOSE_COMPLETES_SERIES = FALSE
  )
  
  res <- evaluate_visits(visits, antigens, verbose = FALSE)
  
  expect_s3_class(res, "data.table")
  expect_true(all(c("VISIT_DATE", "ANTIGEN", "DUE", "GIVEN", "MISSED") %in% names(res)))
})
