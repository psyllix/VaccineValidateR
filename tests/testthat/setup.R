# Ensure CVX map is loaded before tests
if (!exists("cvx") || is.null(cvx)) {
  suppressMessages({
    cvx <<- VaccineValidateR::build_cvx_map()
  })
}
assign("cvx", cvx, envir = .GlobalEnv)