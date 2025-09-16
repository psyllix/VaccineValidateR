.datatable.aware = TRUE

#' Summarize Visit-Level Immunization Evaluations
#'
#' Summarize output from \code{evaluate_visits()}, producing counts and
#' percentages of visits where vaccines were due, given, or missed. Grouping
#' can be flexibly specified with flags and additional grouping columns.
#'
#' @param visit_output A \code{data.table} returned by \code{evaluate_visits()}.
#'   Must contain standardized columns defined in \code{VISIT_RETURN_COLUMNS},
#'   including \code{ANTIGEN}, \code{VISIT_DATE}, \code{DUE}, \code{GIVEN},
#'   and \code{MISSED}. Must also have attribute \code{"processed"} set.
#' @param antigens_of_interest Character vector of antigen names to include.
#'   Default is \code{NULL}, in which case all antigens in the data are used.
#' @param by_year Logical; if \code{TRUE}, group by VISIT_YEAR. Default \code{TRUE}.
#' @param by_antigen Logical; if \code{TRUE}, group by ANTIGEN. Default \code{TRUE}.
#' @param group_cols Optional character vector of additional grouping column names.
#' @param percent_var Character; which percentage column to pivot on in wide format.
#'   Default is \code{"PCT_GIVEN_IF_DUE"}. Other valid options are
#'   \code{"PCT_MISSED_IF_DUE"}, \code{"PCT_WITHIN_15_DAYS_IF_MISSED"},
#'   \code{"PCT_WITHIN_30_DAYS_IF_MISSED"}, \code{"PCT_WITHIN_90_DAYS_IF_MISSED"}.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed.
#'   Default is \code{TRUE}.
#' @param output_format Either \code{"long"} (default; stacked rows with antigen
#'   labels) or \code{"wide"} (one row per group, antigens as columns).
#'
#' @return A \code{data.table} summarizing the requested antigens:
#'   \itemize{
#'     \item In "long" format: one row per group with counts, percents, and labels.
#'     \item In "wide" format: one row per grouping set, columns for each antigen.
#'   }
#'
#' @examples
#' # Summarize all antigens, grouped by year and antigen
#' summarize_visit_evaluations(visit_output)
#'
#' # Group only by SYSTEM and SITE (no year, no antigen), wide format by missed%
#' summarize_visit_evaluations(visit_output,
#'   by_year = FALSE, by_antigen = FALSE, group_cols = c("SYSTEM","SITE"),
#'   output_format = "wide", percent_var = "PCT_MISSED_IF_DUE")
#'
#' @export
summarize_visit_evaluations <- function(visit_output,
                                        antigens_of_interest = NULL,
                                        by_year = TRUE,
                                        by_antigen = TRUE,
                                        group_cols = NULL,
                                        percent_var = "PCT_GIVEN_IF_DUE",
                                        verbose  = TRUE,
                                        output_format = c("long", "wide")) {
  output_format <- match.arg(output_format)
  
  dt <- data.table::setDT(data.table::copy(visit_output))
  
  # Check for "processed" attribute
  if (is.null(attr(dt, "processed"))) {
    stop("visit_output does not have attribute 'processed'. ",
         "Please run evaluate_visits() before summarizing.")
  }
  
  # Ensure VISIT_YEAR exists if needed
  if (by_year && !"VISIT_YEAR" %in% names(dt)) {
    dt[, VISIT_YEAR := data.table::year(VISIT_DATE)]
  }
  
  # Filter antigens if requested
  if (!is.null(antigens_of_interest)) {
    missing <- setdiff(antigens_of_interest, unique(dt$ANTIGEN))
    if (length(missing) > 0) {
      warning("Some requested antigens not found in data: ",
              paste(missing, collapse = ", "))
    }
    dt <- dt[ANTIGEN %in% antigens_of_interest]
  } else {
    antigens_of_interest <- unique(dt$ANTIGEN)
    if (verbose) message("No antigens_of_interest provided. Using all existing antigens: ",
                         paste(antigens_of_interest, collapse = ", "))
  }
  
  if (verbose) message("Starting visit-level summarization...")
  
  # Build grouping variables
  group_vars <- character()
  if (by_year)    group_vars <- c(group_vars, "VISIT_YEAR")
  if (by_antigen) group_vars <- c(group_vars, "ANTIGEN")
  if (!is.null(group_cols)) {
    missing_groups <- setdiff(group_cols, names(dt))
    if (length(missing_groups) > 0) {
      stop("Grouping columns not found in data: ",
           paste(missing_groups, collapse = ", "))
    }
    group_vars <- c(group_vars, group_cols)
  }
  if (length(group_vars) == 0) {
    stop("No grouping variables selected. Enable by_year, by_antigen, or provide group_cols.")
  }
  
  # Summarize
  summary_dt <- dt[, .(
    VISITS = .N,
    DUE    = sum(DUE, na.rm = TRUE),
    GIVEN  = sum(GIVEN, na.rm = TRUE),
    MISSED = sum(MISSED, na.rm = TRUE),
    DELAYED_PRIOR = sum(DELAYED_PRIOR_DOSE, na.rm = TRUE),
    DELAYED_VISIT = sum(DELAYED_VISIT_DOSE, na.rm = TRUE),
    DELAYED_NEXT  = sum(DELAYED_NEXT_DOSE, na.rm = TRUE),
    GIVEN_WITHIN_15_DAYS = sum(MISSED & NEXT_RCVD_15_DAYS, na.rm = TRUE),
    GIVEN_WITHIN_30_DAYS = sum(MISSED & NEXT_RCVD_30_DAYS, na.rm = TRUE),
    GIVEN_WITHIN_90_DAYS = sum(MISSED & NEXT_RCVD_90_DAYS, na.rm = TRUE)
  ), by = group_vars][order(get(group_vars[1]))]
  
  # Add proportions
  summary_dt[, PCT_GIVEN_IF_DUE := data.table::fifelse(DUE > 0, GIVEN / DUE, NA_real_)]
  summary_dt[, PCT_MISSED_IF_DUE := data.table::fifelse(DUE > 0, MISSED / DUE, NA_real_)]
  summary_dt[, PCT_WITHIN_15_DAYS_IF_MISSED := data.table::fifelse(MISSED > 0, GIVEN_WITHIN_15_DAYS / MISSED, NA_real_)]
  summary_dt[, PCT_WITHIN_30_DAYS_IF_MISSED := data.table::fifelse(MISSED > 0, GIVEN_WITHIN_30_DAYS / MISSED, NA_real_)]
  summary_dt[, PCT_WITHIN_90_DAYS_IF_MISSED := data.table::fifelse(MISSED > 0, GIVEN_WITHIN_90_DAYS / MISSED, NA_real_)]
  
  # Validate percent_var choice
  if (!percent_var %in% names(summary_dt)) {
    stop("percent_var '", percent_var, "' not found in summary. ",
         "Valid options are one of: ",
         paste(grep("^PCT_", names(summary_dt), value = TRUE), collapse = ", "))
  }
  
  if (output_format == "wide") {
    if (verbose) message("Pivoting to wide format with data.table::dcast on ", percent_var, "...")
    formula_rhs <- paste(group_vars[group_vars != "ANTIGEN"], collapse = " + ")
    formula_str <- paste0(formula_rhs, " ~ ANTIGEN")
    results_wide <- data.table::dcast(
      summary_dt,
      as.formula(formula_str),
      value.var = percent_var
    )
    return(results_wide[])
  }
  
  if (verbose) message("Completed summarization at ", lubridate::now())
  return(summary_dt[])
}
