# tests/testthat/test-validate.R
test_that("validate_immunizations processes toy data", {
  toy <- data.table::data.table(
    STUDY_ID = 1,
    CVX = c("08", "08"), # HepB
    DOB = as.Date("2020-01-01"),
    PRODUCT = c("HEPATITIS B", "HEPATITIS B"), # required by function
    DATE_GIVEN = as.Date(c("2020-01-01", "2020-02-01")),
    AGE_IMM_GIVEN = c(0, 31)
  )
  
  res <- validate_immunizations(toy, verbose = FALSE)
  
  expect_type(res, "list")
  expect_true(all(c("immunizations", "antigens", "skipped_antigens") %in% names(res)))
  expect_s3_class(res$antigens, "data.table")
  expect_gt(nrow(res$antigens), 0)
})
