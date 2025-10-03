# R/constants.R

#' Vaccine schedule reference dates
#'
#' Important implementation or guidance change dates used in dose validation.
#'
#' @format Constants of class \code{Date}.
#' \describe{
#'   \item{RSV_DATE}{2023-09-01, Beyfortus release (prior = palivizumab).}
#'   \item{MENB_DATE}{2024-10-24, MenB guidance change.}
#'   \item{OPV_DATE}{2016-04-01, OPV rule update.}
#'   \item{HPV_DATE}{2016-04-01, HPV schedule update.}
#'   \item{PCV_DATE}{2010-02-24, PCV13 release.}
#' }
#' @keywords internal
RSV_DATE   <- as.Date("2023-09-01")  # Beyfortus release, prior doses are palivizumab
MENB_DATE  <- as.Date("2024-10-24")  # MenB guidance change: 6mo apart or 3-dose 0,1,5mo
OPV_DATE   <- as.Date("2016-04-01")  # OPV rules updated; older IPV doses not validated
HPV_DATE   <- as.Date("2016-04-01")  # HPV schedule update for 3 dose to 2 dose, not relevant for calculations
PCV_DATE   <- as.Date("2010-02-24")  # PCV13 release; sites transitioned from PCV10
#' @format Numeric scalar.
#' @keywords internal
YEAR_LENGTH <- 365.25
MIN_INTERVAL_DEFAULT <- 24#applies grace period
MIN_INTERVAL_LIVE <- 28
MIN_INTERVAL_COVID_PFIZER <-17

#' Supported antigens
#'
#' List of antigens the package is capable of interpreting. Used internally.
#'
#' @return Character vector of antigen names.
#' @keywords internal
SYSTEM_ANTIGENS <-c('POLIO','HIB','PCV','HEPA','HEPB','HPV','MCV','MENB','TETANUS','ROTA','COVID','RSV','MMR','VZV','INFLUENZA')

#' Seasonal Antigens
#'
#' Internal constant listing antigens that follow seasonal rules
#' (e.g., influenza, COVID, RSV).
#'
#' @format A character vector of antigen names.
#' @examples
#' VaccineValidateR:::SEASONAL_ANTIGENS
#'
#' @keywords internal
SEASONAL_ANTIGENS <- c("INFLUENZA", "COVID", "RSV")

#' CVX Mapping
#'
#' A named list of CVX codes grouped by antigen.
#'
#' This object is built automatically at package load by \code{.onLoad()} 
#' using \code{\link{build_cvx_map}} and provides convenient access to 
#' antigen-specific CVX code sets (e.g., \code{CVX$HIB}, \code{CVX$POLIO}).
#'
#' @details
#' The mapping combines reference CVX–antigen definitions bundled in 
#' \code{extdata/} with optional local LIM→CVX mapping if provided by 
#' the user at build time.
#'
#' @format A named list of integer vectors. Each element corresponds 
#'   to an antigen or grouped set of antigens (e.g., \code{HIB}, 
#'   \code{PCV}, \code{PPV23}, \code{LIVE_NON_ENTERAL}).
#'
#' @examples
#' # Access Hib CVX codes
#' CVX$HIB
#'
#' # Access all group names
#' names(CVX)
#'
#' @seealso \code{\link{build_cvx_map}}, \code{\link{supported_antigens}}
#'
#' @usage CVX
#' @export
cvx <- NULL

#' Antigen return schema
#'
#' Constants defining the expected structure of antigen-level validation
#' outputs, including column names and human-readable definitions.
#'
#' @format
#' \describe{
#'   \item{ANTIGEN_RETURN_COLUMNS}{Character vector of all required column names.}
#'   \item{ANTIGEN_RETURN_DEFS}{Named character vector mapping each column to a
#'         short description of its meaning.}
#' }
#'
#' These constants are used internally for validation, consistency checks, and
#' population-level metric calculations.
#'
#' \tabular{ll}{
#'   \strong{Column} \tab \strong{Definition} \cr
#'   STUDY_ID           \tab Unique patient identifier \cr
#'   PRODUCT            \tab Immunization product name \cr
#'   CVX                \tab CVX vaccine code \cr
#'   ANTIGEN            \tab Standardized antigen grouping \cr
#'   TABLE_INDEX        \tab Row index from original processing table \cr
#'   DATE_GIVEN         \tab Date the immunization was given \cr
#'   AGE_IMM_GIVEN      \tab Age in days at immunization \cr
#'   ABS_ADMIN_COUNTER  \tab Absolute order of doses for the patient \cr
#'   ADMIN_COUNTER      \tab Relative order within product grouping \cr
#'   DOSE_COUNTER       \tab Relative order within antigen grouping \cr
#'   INTERVAL           \tab Interval from prior dose (days) \cr
#'   LIVE_INTERVAL      \tab Interval from prior live dose (days) \cr
#'   VALID              \tab TRUE if dose met minimum interval/age rules \cr
#'   NEXT_DOSE_MIN      \tab Earliest allowable date for next dose \cr
#'   NEXT_DOSE_RECOMMENDED \tab Recommended date for next dose \cr
#'   DOSE_COMPLETES_SERIES \tab TRUE if this dose completed the series \cr
#'   SERIES_COMPLETE    \tab TRUE if antigen is complete at this point \cr
#'   DELAYED            \tab TRUE if this dose was delayed \cr
#'   AGE_FIRST_DOSE     \tab Age in days when first dose of series was given \cr
#' }
#'
#' These constants are used internally for validation, consistency checks, and
#' population-level metric calculations.
#'
#' @keywords internal
# Human-readable definitions for each column (optional, can help with metadata/validation)
ANTIGEN_RETURN_DEFS <- c(
  STUDY_ID = "Unique patient identifier",
  PRODUCT = "Immunization product name",
  CVX = "CVX vaccine code",
  ANTIGEN = "Standardized antigen grouping",
  TABLE_INDEX = "Row index from original processing table",
  DATE_GIVEN = "Date the immunization was given",
  AGE_IMM_GIVEN = "Age in days at immunization",
  ABS_ADMIN_COUNTER = "Absolute order of doses for the patient",
  ADMIN_COUNTER = "Relative order within product grouping",
  DOSE_COUNTER = "Relative order within antigen grouping",
  INTERVAL = "Interval from prior dose (days)",
  LIVE_INTERVAL = "Interval from prior live dose (days)",
  VALID = "TRUE if dose met minimum interval/age rules",
  NEXT_DOSE_MIN = "Earliest allowable date for next dose",
  NEXT_DOSE_RECOMMENDED = "Recommended date for next dose",
  DOSE_COMPLETES_SERIES = "TRUE if this dose completed the series",
  SERIES_COMPLETE = "TRUE if antigen is complete at this point",
  DELAYED = "TRUE if this dose was delayed",
  AGE_FIRST_DOSE = "Age in days when first dose of series was given"
)
# Universal antigen return columns
ANTIGEN_RETURN_COLUMNS <- c(
  'STUDY_ID','PRODUCT','CVX','ANTIGEN','TABLE_INDEX',
  'DATE_GIVEN','AGE_IMM_GIVEN',
  'ABS_ADMIN_COUNTER','ADMIN_COUNTER','DOSE_COUNTER','INTERVAL',
  'LIVE_INTERVAL','VALID','NEXT_DOSE_MIN','NEXT_DOSE_RECOMMENDED',
  'DOSE_COMPLETES_SERIES','SERIES_COMPLETE','DELAYED','AGE_FIRST_DOSE'
)

#' Invalid Dose Columns
#'
#' Standardized column set for invalid dose records. Used consistently
#' across functions that identify, store, and report invalid or
#' excluded immunization doses.
#'
#' @format A character vector of column names:
#' \describe{
#'   \item{STUDY_ID}{Unique patient identifier}
#'   \item{PRODUCT}{Immunization product name}
#'   \item{CVX}{CVX vaccine code}
#'   \item{ANTIGEN}{Standardized antigen grouping}
#'   \item{TABLE_INDEX}{Row index from original processing table}
#'   \item{DATE_GIVEN}{Date the immunization was given}
#'   \item{AGE_IMM_GIVEN}{Age in days at immunization}
#'   \item{ABS_ADMIN_COUNTER}{Absolute order of doses for the patient}
#'   \item{INTERVAL}{Interval from prior dose (days)}
#'   \item{LIVE_INTERVAL}{Interval from prior live dose (days)}
#'   \item{COUNTER}{Internal counter for validation cycles}
#'   \item{CYCLE}{Cycle ID if dose was processed in multiple iterations}
#'   \item{SERIES_COMPLETE}{TRUE if antigen series marked complete}
#'   \item{NOTES}{Free-text notes on why the dose was flagged invalid}
#' }
#'
#' @keywords internal
INVALID_DOSE_COLUMNS <- c(
  'STUDY_ID','PRODUCT','CVX','ANTIGEN','TABLE_INDEX',
  'DATE_GIVEN','AGE_IMM_GIVEN',
  'ABS_ADMIN_COUNTER','INTERVAL','LIVE_INTERVAL',
  'COUNTER','CYCLE','SERIES_COMPLETE','NOTES'
)
# Human-readable definitions for each invalid dose column
INVALID_DOSE_DEFS <- c(
  STUDY_ID        = "Unique patient identifier",
  PRODUCT         = "Immunization product name",
  CVX             = "CVX vaccine code",
  ANTIGEN         = "Standardized antigen grouping",
  TABLE_INDEX     = "Row index from original processing table",
  DATE_GIVEN      = "Date the immunization was given",
  AGE_IMM_GIVEN   = "Age in days at immunization",
  ABS_ADMIN_COUNTER = "Absolute order of doses for the patient",
  INTERVAL        = "Interval from prior dose (days)",
  LIVE_INTERVAL   = "Interval from prior live dose (days)",
  COUNTER         = "Internal counter for validation cycles",
  CYCLE           = "Cycle ID if dose was processed in multiple iterations",
  SERIES_COMPLETE = "TRUE if antigen series marked complete",
  NOTES           = "Free-text notes on why the dose was flagged invalid"
)

#' Visit-level return columns
#'
#' Names of the standard columns returned by \code{evaluate_visits()}.
#'
#' @format A character vector of column names.
#' @keywords internal
VISIT_RETURN_COLUMNS <- c(
  'STUDY_ID','VISIT_ID','VISIT_DATE','ANTIGEN',
  'IS_VIRAL_SEASON','DOB','AGED_OUT','COMPLETED_PREVIOUSLY',
  'NEXT_DOSE_MIN','NEXT_DOSE_RECOMMENDED','DUE','GIVEN','MISSED',
  'DELAYED_PRIOR_DOSE','DELAYED_VISIT_DOSE','DELAYED_NEXT_DOSE',
  'NEXT_GIVEN_DATE','NEXT_RCVD_15_DAYS','NEXT_RCVD_30_DAYS','NEXT_RCVD_90_DAYS'
)

#' Visit-level column definitions
#'
#' Human-readable definitions of each column in \code{VISIT_RETURN_COLUMNS}.
#'
#' @format A named character vector where names correspond to column names
#'   and values are descriptions.
#' @keywords internal
VISIT_RETURN_DEFS <- c(
  STUDY_ID = "Unique identifier for the patient/subject",
  VISIT_ID        = "Unique visit identifier",
  VISIT_DATE = "Date of the clinical visit",
  ANTIGEN = "Vaccine antigen being evaluated",
  IS_VIRAL_SEASON = "Flag for viral season (Sept to Apr = TRUE)",
  DOB = "Date of birth of patient",
  AGED_OUT = "Flag if patient has exceeded eligible age for vaccine",
  COMPLETED_PREVIOUSLY = "TRUE if the series was completed before this visit",
  NEXT_DOSE_MIN = "Earliest allowable date for the next dose",
  NEXT_DOSE_RECOMMENDED = "Recommended date for the next dose",
  DUE = "TRUE if the vaccine dose is due at this visit",
  GIVEN = "TRUE if the vaccine dose was given at this visit",
  MISSED = "TRUE if dose was due but not given at this visit",
  DELAYED_PRIOR_DOSE = "TRUE if the previous dose was given late",
  DELAYED_VISIT_DOSE = "TRUE if the dose given at this visit was delayed",
  DELAYED_NEXT_DOSE = "TRUE if the next scheduled dose was given late",
  NEXT_GIVEN_DATE = "Date of the next dose given (if available)",
  NEXT_RCVD_15_DAYS="TRUE if an immuniztion was given within 15 days of a visit, not inclusive of visit date",
  NEXT_RCVD_30_DAYS="TRUE if an immuniztion was given within 30 days of a visit, not inclusive of visit date",
  NEXT_RCVD_90_DAYS="TRUE if an immuniztion was given within 90 days of a visit, not inclusive of visit date"
)