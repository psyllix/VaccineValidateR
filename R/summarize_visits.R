.datatable.aware = TRUE

#' Summarize Visit-Level Immunization Evaluations
#'
#' Summarize output from \code{evaluate_visits()}, producing counts and
#' percentages of visits where vaccines were due, given, or missed.
#' Grouping can be flexibly specified by time period and user-defined columns.
#'
#' @param visit_output A \code{data.table} returned by \code{evaluate_visits()}.
#'   Must contain standardized columns defined in \code{VISIT_RETURN_COLUMNS},
#'   including \code{ANTIGEN}, \code{VISIT_DATE}, \code{DUE}, \code{GIVEN},
#'   and \code{MISSED}. Must also have attribute \code{"processed"} set.
#' @param antigens_of_interest Character vector of antigen names to include.
#'   Default is \code{NULL}, in which case all antigens in the data are used.
#' @param time_slice Character; defines how to group visits by time period.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"NONE"} — no time-based grouping (all visits grouped together).
#'     \item \code{"YEAR"} — calendar year of \code{VISIT_DATE}.
#'     \item \code{"MONTH"} — calendar month of \code{VISIT_DATE}, labeled as YYYY-MM.
#'     \item \code{"BIMONTH"} — two-month intervals (Jan–Feb, Mar–Apr, etc.), labeled as YYYY-MM (ending month).
#'     \item \code{"QUARTER"} — three-month quarters (Q1–Q4), labeled as YYYY-Q\emph{n}.
#'     \item \code{"FOURMONTH"} — four-month periods (Jan–Apr, May–Aug, Sep–Dec), labeled as YYYY-P\emph{n}.
#'     \item \code{"VIRAL"} — viral season grouping (August–July), where visits from
#'       September–December are assigned to the current year and visits from January–July
#'       are assigned to the previous year (e.g., September 2024–April 2025 → 2024 season).
#'   }
#'   Default is \code{"NONE"}.
#' @param group_cols Optional character vector of additional grouping column names.
#'   These columns must exist in \code{visit_output}.
#' @param first_per_group Logical; if \code{TRUE}, restricts to the first visit
#'   per \code{STUDY_ID} within each grouping set. Default \code{TRUE}.
#' @param percent_var Character; which percentage column to pivot on in wide format.
#'   Default is \code{"PCT_GIVEN_IF_DUE"}. Other valid options are
#'   \code{"PCT_MISSED_IF_DUE"}, \code{"PCT_WITHIN_15_DAYS_IF_MISSED"},
#'   \code{"PCT_WITHIN_30_DAYS_IF_MISSED"}, and \code{"PCT_WITHIN_90_DAYS_IF_MISSED"}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages. Default \code{TRUE}.
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
#' # Summarize all antigens by year and antigen
#' summarize_visits(visit_output)
#'
#' # Summarize by system, site, and viral season
#' summarize_visits(visit_output,
#'   time_slice = "VIRAL",
#'   group_cols = c("SYSTEM","SITE"),
#'   output_format = "wide",
#'   percent_var = "PCT_MISSED_IF_DUE")
#'
#' @export
summarize_visits <- function(visit_output,
                                        antigens_of_interest = NULL,
                                        time_slice = c("NONE", "YEAR", "MONTH","BIMONTH","QUARTER","FOURMONTH","VIRAL"),
                                        group_cols = NULL,
                                        first_per_group = TRUE,
                                        percent_var = "PCT_GIVEN_IF_DUE",
                                        verbose  = TRUE,
                                        output_format = c("long", "wide")) {
  #match arguments to first option
  output_format <- match.arg(output_format)
  time_slice <- match.arg(time_slice)
  #initialize data tables
  dt <- data.table::setDT(data.table::copy(visit_output))
  
  # Check for "processed" attribute
  if (is.null(attr(dt, "processed"))) {
    stop("visit_output does not have attribute 'processed'. ",
         "Please run evaluate_visits() before summarizing.")
  }
  
  # Ensure VISIT_YEAR exists if needed
  # Ensure VISIT_YEAR exists if needed
  if (!"VISIT_YEAR" %in% names(dt)) {
    dt[, VISIT_YEAR := year_from_date(VISIT_DATE)]
  }
  if (!"VISIT_MONTH" %in% names(dt)) {
    dt[, VISIT_MONTH := month_from_date(VISIT_DATE)]
  }
  
  if (time_slice == "YEAR") {
    time_var <- "VISIT_YEAR"
    
  } else if (time_slice == "MONTH") {
    dt[, VISIT_MONTH_LABEL := sprintf("%04d-%02d", VISIT_YEAR, VISIT_MONTH)]
    time_var <- "VISIT_MONTH_LABEL"
    
  } else if (time_slice == "BIMONTH") {
    dt[, VISIT_BIMONTH := sprintf("%04d-%02d", VISIT_YEAR, ceiling(VISIT_MONTH / 2) * 2)]
    time_var <- "VISIT_BIMONTH"
    
  } else if (time_slice == "QUARTER") {
    dt[, VISIT_QUARTER := sprintf("%04d-Q%d", VISIT_YEAR, ceiling(VISIT_MONTH / 3))]
    time_var <- "VISIT_QUARTER"
    
  } else if (time_slice == "FOURMONTH") {
    dt[, VISIT_FOURMONTH := sprintf("%04d-P%d", VISIT_YEAR, ceiling(VISIT_MONTH / 4))] # P = period
    time_var <- "VISIT_FOURMONTH"
    
  } else if (time_slice == "VIRAL") {
    dt[, VISIT_SEASON := fifelse(VISIT_MONTH > 7L, VISIT_YEAR, VISIT_YEAR - 1L)]
    time_var <- "VISIT_SEASON"
    
  } else {
    time_var <- NULL
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
  group_vars <- c(group_vars, "ANTIGEN")#ensure group by antigen
  # Add time slice variable if applicable
  if (!is.null(time_var)) {
    group_vars <- c(group_vars, time_var)
  }
  
  if (!is.null(group_cols)) {
    missing_groups <- setdiff(group_cols, names(dt))
    if (length(missing_groups) > 0) {
      stop("Grouping columns not found in data: ",
           paste(missing_groups, collapse = ", "))
    }
    group_vars <- c(group_vars, group_cols)
  }
  if (length(group_vars) == 0) {
    stop("No grouping variables selected. Use time_slice or provide group_cols.")
  }
  # Collapse to first event in each group if requested
  if (first_per_group) {
    if (verbose) message("Restricting to first event per STUDY_ID within grouping vars.")
    dt <- dt[order(VISIT_DATE)][, .SD[1], by = c("STUDY_ID", group_vars)]
  }
  #key for performance
  data.table::setkeyv(dt, group_vars)
  ###### SUMMARIZATION STEPS ############
  summary_dt <- dt[, .(
    VISITS = .N,
    DUE    = sum(DUE, na.rm = TRUE),
    GIVEN  = sum(GIVEN, na.rm = TRUE),
    MISSED = sum(MISSED, na.rm = TRUE),
    DELAYED_PRIOR = sum(DELAYED_PRIOR_DOSE, na.rm = TRUE),
    DELAYED_VISIT = sum(DELAYED_VISIT_DOSE, na.rm = TRUE),
    DELAYED_NEXT  = sum(DELAYED_NEXT_DOSE, na.rm = TRUE),
    NO_FUTHER_DOSES  = sum(is.na(NEXT_GIVEN_DATE), na.rm = TRUE),
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
    formula_rhs <- paste(setdiff(group_vars, "ANTIGEN"), collapse = " + ")
    formula_str <- paste0(formula_rhs, " ~ ANTIGEN")
    results_wide <- data.table::dcast(summary_dt, as.formula(formula_str), value.var = percent_var)
    return(results_wide[])
  }
  
  if (verbose) message("Completed summarization at ", Sys.time())
  return(summary_dt[])
}
