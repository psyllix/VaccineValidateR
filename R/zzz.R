# R/zzz.R
# Package lifecycle hooks

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Ensure CVX mapping is available as a package constant
  package_env <- parent.env(environment())
  package_env$cvx <- build_cvx_map()
}