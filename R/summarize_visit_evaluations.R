.datatable.aware = TRUE

#' Summarize Visit-Level Immunization Evaluations
#'
#' This function summarizes the output from \code{evaluate_visits()}, producing
#' counts and percentages of visits where vaccines were due, given, or missed,
#' grouped by antigen, year, and optionally a grouping column (e.g., system or subset).
#'
#' @param visit_output A \code{data.table} returned by \code{evaluate_visits()}.
#'   Must contain standardized columns defined in \code{VISIT_RETURN_COLUMNS},
#'   including \code{ANTIGEN}, \code{VISIT_DATE}, \code{DUE}, \code{GIVEN},
#'   and \code{MISSED}.
#' @param antigens_of_interest Character vector of antigen names to include.
#'   Default is \code{NULL}, in which case all antigens in the data are used.
#' @param group_col Optional column name in \code{visit_output} to group by
#'   (e.g., "SYSTEM" or "SITE"). Default is \code{NULL}, meaning no grouping.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed.
#'   Default is \code{TRUE}.
#' @param output_format Either \code{"long"} (default; stacked rows with antigen
#'   labels) or \code{"wide"} (one row per year, antigens as columns).
#'
#' @return A \code{data.table} summarizing the requested antigens:
#'   \itemize{
#'     \item In "long" format: one row per antigen-year-group with counts,
#'       percents, and antigen labels.
#'     \item In "wide" format: one row per evaluation year (and group if
#'       specified), columns for each antigen.
#'   }
#'
#' @examples
#' # Suppose visit_output is from evaluate_visits()
#' summarize_visit_evaluations(visit_output)
#'
#' # Only summarize Polio and Hib, grouped by SYSTEM
#' summarize_visit_evaluations(visit_output,
#'                             antigens_of_interest = c("POLIO","HIB"),
#'                             group_col = "SYSTEM")
#'
#' @export
summarize_visit_evaluations <- function(visit_output,
                                        antigens_of_interest = NULL,
                                        group_col = NULL,
                                        verbose  = TRUE,
                                        output_format = c("long", "wide")) {
  output_format <- match.arg(output_format)
  
  # Defensive copy
  dt <- data.table::setDT(data.table::copy(visit_output))
  
  # Ensure VISIT_YEAR exists
  if (!"VISIT_YEAR" %in% names(dt)) {
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
  
  # Build grouping vars
  group_vars <- c("VISIT_YEAR", "ANTIGEN")
  if (!is.null(group_col) && group_col %in% names(dt)) {
    group_vars <- c(group_vars, group_col)
  }
  else{stop("Group by column missing. Group by column is Case Sensitive.")}
  
  # Summarize by antigen, year, and optional group
  summary_dt <- dt[, .(
    VISITS = .N,
    DUE    = sum(DUE, na.rm = TRUE),
    GIVEN  = sum(GIVEN, na.rm = TRUE),
    MISSED = sum(MISSED, na.rm = TRUE),
    DELAYED_PRIOR = sum(DELAYED_PRIOR_DOSE, na.rm = TRUE),
    DELAYED_VISIT = sum(DELAYED_VISIT_DOSE, na.rm = TRUE),
    DELAYED_NEXT  = sum(DELAYED_NEXT_DOSE, na.rm = TRUE),
    GIVEN_WITHIN_15_DAYS = sum(MISSED&NEXT_RCVD_15_DAYS, na.rm = TRUE),
    GIVEN_WITHIN_30_DAYS = sum(MISSED&NEXT_RCVD_30_DAYS, na.rm = TRUE),
    GIVEN_WITHIN_90_DAYS = sum(MISSED&NEXT_RCVD_90_DAYS, na.rm = TRUE)
  ), by = group_vars][order(VISIT_YEAR, ANTIGEN)]
  
  # Add proportions
  summary_dt[, PCT_GIVEN_IF_DUE := data.table::fifelse(DUE > 0, GIVEN / DUE, NA_real_)]#can be given when not due!
  summary_dt[, PCT_MISSED_IF_DUE := data.table::fifelse(DUE > 0, MISSED / DUE, NA_real_)]
  summary_dt[, PCT_WITHIN_15_DAYS_IF_MISSED := data.table::fifelse(MISSED > 0, (GIVEN_WITHIN_15_DAYS) / MISSED, NA_real_)]
  summary_dt[, PCT_WITHIN_30_DAYS_IF_MISSED := data.table::fifelse(MISSED > 0, (GIVEN_WITHIN_30_DAYS) / MISSED, NA_real_)]
  summary_dt[, PCT_WITHIN_90_DAYS_IF_MISSED := data.table::fifelse(MISSED > 0, (GIVEN_WITHIN_90_DAYS) / MISSED, NA_real_)]
  
  if (output_format == "wide") {
    if (verbose) message("Pivoting to wide format with data.table::dcast...")
    
    if (!is.null(group_col) && group_col %in% names(dt)) {
      results_wide <- data.table::dcast(
        summary_dt,
        paste0("VISIT_YEAR + ", group_col, " ~ ANTIGEN"),
        value.var = "PCT_GIVEN_IF_DUE"
      )
    } else {
      results_wide <- data.table::dcast(
        summary_dt,
        VISIT_YEAR ~ ANTIGEN,
        value.var = "PCT_GIVEN_IF_DUE"
      )
    }
    return(results_wide[])
  }
  
  if (verbose) message("Completed summarization at ", lubridate::now())
  return(summary_dt[])
}
