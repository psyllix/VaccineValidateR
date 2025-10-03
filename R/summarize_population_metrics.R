#' Summarize Population-Level Immunization Metrics
#'
#' This function summarizes patient-level population immunization metrics
#' (created by \code{population_metrics()}) across evaluation years, optionally
#' returning results in long or wide format. Metrics are grouped by evaluation
#' year (e.g., age 2 for CIS, age 13 for ADOL).
#'
#' @param patients A \code{data.table} of patient-level records that has already
#'   been processed by \code{population_metrics()} (contains metric columns and
#'   attached attribute \code{population_metric_definitions}).
#' @param visits A \code{data.table} of patient visits including at least
#'   \code{STUDY_ID} and \code{VISIT_DATE}.
#' @param years A numeric vector of evaluation years to include in the summary.
#' @param metrics_of_interest Character vector of metric names (e.g.,
#'   \code{"HEDIS_ADOL2"}, \code{"UTD_CIS10"}). These must exist as columns in
#'   \code{patients}. If not specified, will return all possible metrics on the patient table.
#' @param id_col Column name in \code{patients} identifying unique patients.
#'   Default is \code{"STUDY_ID"}.
#' @param dob_col Column name in \code{patients} containing date of birth.
#'   Default is \code{"DOB"}.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed.
#'   Default is \code{TRUE}.
#' @param output_format Either \code{"long"} (default; stacked rows with metric
#'   labels) or \code{"wide"} (one row per year, metrics as columns).
#'
#' @return A \code{data.table} summarizing the requested metrics:
#'   \itemize{
#'     \item In "long" format: one row per metric-year with counts, percents,
#'       and metric labels.
#'     \item In "wide" format: one row per evaluation year, columns for each
#'       metric.
#'   }
#'
#' @details
#' Evaluation age is inferred from the metric name:
#' \itemize{
#'   \item CIS metrics (\code{"CIS"}) → evaluated at age 2
#'   \item ADOL metrics (\code{"ADOL"}) → evaluated at age 13
#'   \item All other metrics → default to age 18
#' }
#'
#' @export
summarize_population_metrics <- function(patients, 
                                         visits, 
                                         years, 
                                         metrics_of_interest = NULL,
                                         id_col   = "STUDY_ID", 
                                         dob_col  = "DOB",
                                         verbose  = TRUE,
                                         output_format = c("long", "wide")) {
  output_format <- match.arg(output_format)
  
  # Pull metric definitions from attribute
  definitions <- attr(patients, "population_metric_definitions")
  if (is.null(definitions)) {
    stop("Patient table is missing metric definitions. ",
         "Please run population_metrics() first.")
  }
  # If metrics_of_interest not supplied, default to all available
  if (is.null(metrics_of_interest)) {
    metrics_of_interest <- names(definitions)
    if (verbose) message("metrics_of_interest not provided. Using all available metrics: ",
                         paste(metrics_of_interest, collapse = ", "))
  }
  # Ensure data.table
  patients <- data.table::setDT(data.table::copy(patients))
  visits   <- data.table::setDT(data.table::copy(visits))
  
  if (verbose) message("Starting population metric summarization for years: ",
                       paste(years, collapse = ", "))
  
  
  if (verbose) message("Mapping visit dates... this may take some time.")
  if (verbose) message("Mapping last visit year... Please wait ~1 minute...")
  
  # collapse to one row per id_col
  visit_map <- visits[, .(LAST_VISIT_YEAR = max(VISIT_DATE, na.rm = TRUE),
                          FIRST_VISIT_YEAR = min(VISIT_DATE, na.rm = TRUE)
                          ), by = id_col]
  
  # convert to year just once
  visit_map[, LAST_VISIT_YEAR := data.table::year(LAST_VISIT_YEAR)]
  
  # join back onto patients
  patients <- visit_map[patients, on = id_col]
  
  if (verbose) message("Mapped last visit year onto patients")
  
  # add birth year (vectorized)
  patients[, BIRTH_YEAR := data.table::year(get(dob_col))]
  if (verbose) message("Established birth year")
  results_list <- list()
  for (m in metrics_of_interest) {
    if (!m %in% names(patients)) {
      warning(sprintf("Metric %s not found in patients table, skipping.", m))
      next
    }
    
    # Make a temporary copy to avoid modifying main patients table
    pat_tmp <- data.table::copy(patients)
    
    # Choose evaluation age based on metric
    if (grepl("CIS", m)) {
      eval_age <- 2
    } else if (grepl("ADOL", m)) {
      eval_age <- 13
    } else {
      eval_age <- 18
    }
    
    # Derive evaluation year
    pat_tmp[, EVAL_YEAR := BIRTH_YEAR + eval_age]
    
    # Only keep patients still in cohort at evaluation year
    pat_tmp[, IN_COHORT := (!is.na(LAST_VISIT_YEAR) & LAST_VISIT_YEAR >= EVAL_YEAR&!is.na(LAST_VISIT_YEAR)& FIRST_VISIT_YEAR<EVAL_YEAR)]
    
    # Summarize by evaluation year
    dt <- pat_tmp[IN_COHORT == TRUE & EVAL_YEAR %in% years,
                  .(IS_UTD = sum(get(m), na.rm = TRUE),
                    N      = .N),
                  by = EVAL_YEAR][order(EVAL_YEAR)]
    
    dt[, `:=`(PERCENT = IS_UTD / N, METRIC = m)]
    
    # Attach definition if available
    if (m %in% names(definitions)) {
      data.table::setattr(dt, "definition", definitions[[m]])
    }
    
    results_list[[m]] <- dt
    
    if (verbose) {
      def_msg <- if (m %in% names(definitions)) paste0(" (", definitions[[m]], ")") else ""
      message("Evaluating metric: ", m, def_msg)
      for (yr in dt$EVAL_YEAR) {
        message("  Year ", yr, ": ",
                dt[EVAL_YEAR == yr, IS_UTD], "/", 
                dt[EVAL_YEAR == yr, N], " (",
                round(100 * dt[EVAL_YEAR == yr, PERCENT], 1), "%)")
      }
    }
  }
  
  # Combine into one tall table
  results_long <- data.table::rbindlist(results_list, fill = TRUE, use.names = TRUE)
  
  if (output_format == "wide") {
    if (verbose) message("Pivoting to wide format using data.table::dcast...")
    results_wide <- data.table::dcast(
      results_long,
      EVAL_YEAR ~ METRIC,
      value.var = "PERCENT"
    )
    return(results_wide)
  }
  
  if (verbose) message("Completed summarization at ", Sys.time())
  return(results_long)
}
