# R/utils.R
# Internal helper functions for VaccineValidateR
# Do not export

.datatable.aware <- TRUE

#' Convert weeks to days with optional grace period
#' @param x Number of weeks
#' @param grace Days to subtract (default 4)
#' @return Days
#' @keywords internal
wk_with_grace <- function(x, grace = 4) {
  x * 7 - grace
}

#' Convert months to days with optional grace period
#' @param x Number of months
#' @param grace Days to subtract (default 4)
#' @return Days
#' @keywords internal
mon_with_grace <- function(x, grace = 4) {
  floor(YEAR_LENGTH/ 12 * x) - grace
}

#' Convert years to days with optional grace period
#' @param x Number of years
#' @param grace Days to subtract (default 4)
#' @return Days
#' @keywords internal
yr_with_grace <- function(x, grace = 4) {
  floor(YEAR_LENGTH* x) - grace
}

#' Identity wrappers for no-grace calculations
#' @keywords internal
wk_no_grace <- function(x) wk_with_grace(x, grace = 0)

#' @keywords internal
mon_no_grace <- function(x) mon_with_grace(x, grace = 0)

#' @keywords internal

#' @keywords internal
yr_no_grace <- function(x) yr_with_grace(x, grace = 0)

#' @keywords internal
days_no_grace <- function(x) x

#' @keywords internal
days_with_grace <- function(x, grace = 4) x - grace

#' Safe min with NA handling
#' @param x Numeric vector
#' @param na.rm Logical, remove NAs (default TRUE)
#' @return Minimum or NA
#' @keywords internal
min_na <- function(x, na.rm = TRUE) {
  if (length(x) > 0) min(x, na.rm = na.rm) else NA
}

#' Safe max with NA handling
#' @param x Numeric vector
#' @param na.rm Logical, remove NAs (default TRUE)
#' @return Maximum or NA
#' @keywords internal
max_na <- function(x, na.rm = TRUE) {
  if (length(x) > 0) max(x, na.rm = na.rm) else NA
}
#' Safe max with NA handling
#' @param x Numeric vector
#' @param na.rm Logical, remove NAs (default TRUE)
#' @return Maximum or 0
#' @keywords internal
max0<-function(x,na.rm = TRUE){as.numeric(if(!is.infinite(suppressWarnings(max(x,na.rm=TRUE)))){max(x,na.rm=TRUE)}else {0})}


#' Returns the list of Antigens the system is capable of interpreting
#' @export
#' @return Character vector of the system antigens
supported_antigens<-function(){
  SYSTEM_ANTIGENS
}


#' Build CVX Mapping List
#'
#' Utility function to create a structured list of CVX code groups
#' for each antigen of interest. This function is called automatically
#' at package load to populate the global constant \code{\link{CVX}}.
#' 
#' Most users will not need to call this directly. To refresh the
#' global mapping (e.g., with a custom LIM to CVX file), use
#' \code{\link{update_cvx_map}} instead.
#'
#' @param lim_cvx_map Optional \code{data.table} containing local immunization
#'   identifiers (LIM) mapped to CVX codes. Must include at least \code{PRODUCT}
#'   and \code{CVX}.
#'
#' @return A named list of CVX code vectors (cvx), each corresponding to an antigen
#'   or grouped antigen set (e.g., \code{HIB}, \code{POLIO}, \code{COVID},
#'   \code{LIVE_NON_ENTERAL}).
#'
#' @details
#' The mapping combines reference CVX-antigen definitions bundled in 
#' \code{extdata/} with an optional local LIM to CVX mapping provided by
#' the user. The result is assigned to \code{\link{CVX}} when the package
#' loads, making the mapping available globally.
#'
#' @examples
#' # Normally not needed, CVX is available after library load: CVX$HIB
#'
#' # Rebuild manually (not usually recommended)
#' custom_map <- build_cvx_map()
#' custom_map$POLIO
#'
#' # Supported way to refresh the global mapping
#' \dontrun{
#' update_cvx_map()
#' }
#'
#' @seealso \code{\link{cvx}}, \code{\link{update_cvx_map}}, \code{\link{supported_antigens}}
#' @export
build_cvx_map <- function(lim_cvx_map = NULL) {
  # load package reference files
  cvx_files <- list.files(system.file("extdata", package = "VaccineValidateR"),
                          full.names = TRUE)
  cvx_map <- data.table::rbindlist(lapply(cvx_files, data.table::fread))
  
  # merge with user-provided LIM to CVX map if present
  if (!is.null(lim_cvx_map)) {
    cvx_map <- data.table::rbindlist(list(cvx_map, lim_cvx_map),
                                     fill = TRUE, use.names = TRUE)
  }
  
  # construct list of groups
  cvx <- list(
    # core antigens
    TETANUS   = cvx_map[grepl("TD|DT", ANTIGEN, ignore.case = TRUE), CVX],
    POLIO     = cvx_map[grepl("IPV|OPV|POLIO", ANTIGEN, ignore.case = TRUE), CVX],
    HIB       = cvx_map[grepl("HIB", ANTIGEN, ignore.case = TRUE), CVX],
    PCV       = cvx_map[grepl("PREVNAR|PCV", ANTIGEN, ignore.case = TRUE), CVX],
    PPV23     = cvx_map[grepl("PPV23", ANTIGEN, ignore.case = TRUE), CVX],
    HEPB      = cvx_map[grepl("HEPB", ANTIGEN, ignore.case = TRUE), CVX],
    HEPA      = cvx_map[grepl("HEPA", ANTIGEN, ignore.case = TRUE), CVX],
    MMR       = cvx_map[grepl("MMR", ANTIGEN, ignore.case = TRUE), CVX],
    VZV       = cvx_map[grepl("VARICELLA", ANTIGEN, ignore.case = TRUE), CVX],
    HPV       = cvx_map[grepl("HPV", ANTIGEN, ignore.case = TRUE), CVX],
    MCV       = cvx_map[grepl("MCV|MENACPOLY", ANTIGEN, ignore.case = TRUE), CVX],
    MENB      = cvx_map[grepl("MENB", ANTIGEN, ignore.case = TRUE), CVX],
    INFLUENZA = cvx_map[grepl("INFLUENZA", ANTIGEN, ignore.case = TRUE), CVX],
    LIVE_FLU  = cvx_map[grepl("NASALINFLUENZA", ANTIGEN, ignore.case = TRUE), CVX],
    ROTA      = cvx_map[grepl("ROTA", ANTIGEN, ignore.case = TRUE), CVX],
    RSV       = cvx_map[grepl("RSV", ANTIGEN, ignore.case = TRUE), CVX],
    
    # COVID families
    COVID         = cvx_map[grepl("COVID19", ANTIGEN, ignore.case = TRUE), CVX],
    COVID_BIV     = cvx_map[grepl("COVID19BIVALENT", ANTIGEN, ignore.case = TRUE), CVX],
    COVID_2023    = cvx_map[grepl("COVID2023", ANTIGEN, ignore.case = TRUE), CVX],
    MODERNA       = cvx_map[grepl("MODERNA", ANTIGEN, ignore.case = TRUE), CVX],
    PFIZER        = cvx_map[grepl("PFIZER", ANTIGEN, ignore.case = TRUE), CVX],
    JANSSEN       = cvx_map[grepl("JANSSEN", ANTIGEN, ignore.case = TRUE), CVX],
    AZ_COVID      = cvx_map[grepl("COVID19RECOMBINANTMULTI", ANTIGEN, ignore.case = TRUE), CVX],
    
    # specialized sets
    ORAL_POLIO    = cvx_map[grepl("OPV", ANTIGEN, ignore.case = TRUE), CVX],
    HIB3          = cvx_map[grepl("HIB3", ANTIGEN, ignore.case = TRUE), CVX],
    HIB4          = cvx_map[grepl("HIB", ANTIGEN, ignore.case = TRUE) & 
                              !grepl("HIB3", ANTIGEN, ignore.case = TRUE), CVX],
    ROTA2         = cvx_map[grepl("ROTA2DOSE", ANTIGEN, ignore.case = TRUE), CVX],
    ADULT_HEPA    = cvx_map[grepl("ADULTHEPA", ANTIGEN, ignore.case = TRUE), CVX],
    ADOL_HEPB2    = cvx_map[grepl("ADOLHEPB", ANTIGEN, ignore.case = TRUE), CVX],
    ADULT_HEPB    = cvx_map[grepl("ADULTHEPB", ANTIGEN, ignore.case = TRUE), CVX],
    ADULT_HEPB2   = cvx_map[CVX %in% c(189), CVX], # 2-dose adult HepB
    BEXSERO       = cvx_map[grepl("MENBBEXSERO", ANTIGEN, ignore.case = TRUE), CVX],
    TRUMENBA      = cvx_map[grepl("MENBTRUMENBA", ANTIGEN, ignore.case = TRUE), CVX],
    YELLOWFEVER   = cvx_map[grepl("YELLOWFEVER", ANTIGEN, ignore.case = TRUE), CVX],
    SMALLPOX      = cvx_map[grepl("SMALLPOX", ANTIGEN, ignore.case = TRUE), CVX],
    
    # grouped sets
    LIVE_NON_ENTERAL = cvx_map[grepl("MMR|VARICELLA|NASALINFLUENZA|YELLOWFEVER|SMALLPOX", 
                                     ANTIGEN, ignore.case = TRUE), CVX]
  )
  
  return(cvx)
}
#' Update the CVX Mapping
#'
#' Rebuilds the global \code{\link{CVX}} constant with an optional
#' local LIM to CVX map. This updates the copy inside the package
#' namespace so that all downstream functions will see the new
#' mapping.
#'
#' @param lim_cvx_map Optional \code{data.table} containing local
#'   immunization identifiers (LIM) mapped to CVX codes. Must include
#'   at least \code{PRODUCT} and \code{CVX}.
#'
#' @return Invisibly returns the updated \code{\link{CVX}} object.
#' @examples
#' # Refresh CVX mapping with default references
#' update_cvx_map()
#'
#' # Refresh with custom local LIM to CVX mapping
#' my_map <- data.table::data.table(LIM = "X123", PRODUCT = "Test", CVX = 999)
#' update_cvx_map(my_map)
#'
#' # Verify it updated
#' cvx$CUSTOM <- 999
#'
#' @seealso \code{\link{cvx}}, \code{\link{build_cvx_map}}
#' @export
update_cvx_map <- function(lim_cvx_map = NULL) {
  new_map <- build_cvx_map(lim_cvx_map)
  assign("cvx", new_map, envir = parent.env(environment()))
  invisible(new_map)
}

#' Validate Antigen Data Table
#'
#' Internal utility to confirm that a provided antigen-level \code{data.table}
#' matches the expected schema from \code{validation()}.
#'
#' Checks that the table is a \code{data.table}, contains the \code{VALID}
#' column, and includes all required fields from
#' \code{\link{ANTIGEN_RETURN_COLUMNS}}.
#'
#' @param antigens A \code{data.table} of antigen-level immunization records,
#'   typically from \code{validation()}.
#'
#' @return Invisibly returns \code{TRUE} if validation passes.
#'   Otherwise, throws an error with a list of missing columns.
#'
#' @keywords internal
#' @examples
#' check_antigen_table(antigens)
check_antigen_table <- function(antigens) {
  if (!"data.table" %in% class(antigens) || !"VALID" %in% names(antigens)) {
    stop("Antigen data has not been validated (missing VALID column).")
  }
  missing_cols <- setdiff(ANTIGEN_RETURN_COLUMNS, names(antigens))
  if (length(missing_cols) > 0) {
    stop("Antigen data must come from validation() output. Missing columns: ",
         paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}
#' Build Antigen List
#'
#' Utility function to derive the set of antigens to evaluate, given a target list
#' and the system defaults. Handles "ALL" and live-virus logic.
#'
#' @param antigen_to_eval Character vector of requested antigens (or "ALL").
#'   Defaults to "ALL".
#' @param verbose Logical; whether to message progress. Default = TRUE.
#'
#' @return Character vector of antigens to evaluate. If none, returns character(0).
#' @keywords internal
build_antigen_list <- function(antigen_to_eval = c("ALL"), verbose = TRUE) {
  # start with system list
  antigens_list <- SYSTEM_ANTIGENS
  antigen_to_eval <- toupper(antigen_to_eval)#upperfy the data just in case
  # restrict to requested antigens unless "ALL"
  if (!"ALL" %in% antigen_to_eval) {
    antigens_list <- intersect(antigens_list, antigen_to_eval)
  }
  
  # ensure all live-virus antigens are included if any live vaccine is requested
  if (any(c("MMR", "VZV", "INFLUENZA") %in% antigens_list)) {
    antigens_list <- c(antigens_list, "MMR", "VZV", "INFLUENZA")
  }
  
  # force unique
  antigens_list <- unique(antigens_list)
  
  # message
  if (verbose) {
    if (length(antigens_list) == 0L) {
      message("No matching antigens requested - empty set will be returned.")
    } else {
      message("Antigen list developed - will evaluate for: ",
              paste(antigens_list, collapse = ", "))
    }
  }
  
  return(antigens_list)
}
#' Compute Next September 1 After a Given Date
#'
#' Given a vector of \code{Date} objects, returns the next September 1
#' on or after each input date. Used internally for immunization schedule
#' calculations (e.g., seasonal influenza logic).
#'
#' @param d A vector of class \code{Date}.
#'
#' @return A vector of \code{Date} values, each corresponding to
#' the September 1 following (or equal to) the year of \code{d}.
#'
#' @examples
#' # Internal use only
#' test_dates <- as.Date(c("2020-01-15", "2020-09-01", "2020-12-15"))
#' VaccineValidateR:::next_sept1(test_dates)
#'
#' @keywords internal
next_sept1 <- function(d) {
  stopifnot(inherits(d, "Date"))
  y <- year_from_date(d)#extract year (actual)
  m <- month_from_date(d)#extract month (actual)
  y_out <- ifelse(m > 8L, y+1, y)#(sept-->next sept)
  as.Date(ISOdate(y_out, 9L, 1L))
}
#' Extract Calendar Year from Date
#'
#' Internal helper to return the calendar year as an integer from a
#' \code{Date} vector. 
#'
#' @param d A vector of class \code{Date}.
#'
#' @return An integer vector of calendar years.
#'
#' @examples
#' VaccineValidateR:::year_from_date(as.Date(c("2000-01-01", "2025-09-18")))
#'
#' @keywords internal
year_from_date <- function(d) {
  stopifnot(inherits(d, "Date"))
  as.POSIXlt(d)$year + 1900L#starts at 1900 as offset because computers are silly
}

#' Extract Calendar Month from Date
#'
#' Internal helper to return the calendar month (1-12) as an integer from a
#' \code{Date} vector.
#' 
#' @param d A vector of class \code{Date}.
#'
#' @return An integer vector of calendar months (1-12).
#'
#' @examples
#' VaccineValidateR:::month_from_date(as.Date(c("2000-01-01", "2025-09-18")))
#'
#' @keywords internal
month_from_date <- function(d) {
  stopifnot(inherits(d, "Date"))
  as.POSIXlt(d)$mon + 1L #this at least makes sense!
}
#' @keywords internal
#' @noRd
make_a_date <- function(DT, col) {
  # Expect col to exist and contain ISO8601 YYYY-MM-DD strings
  DT[, (col) := fasttime::fastDate(as.character(get(col)))]
}

#' Ensure DOB is present in a data.table
#'
#' Adds a properly formatted DOB column to the main table (`primary_dt`),
#' using either an existing DOB column or a supplemental `patients` table.
#'
#' @param primary_dt A \code{data.table} (e.g., immunization_data or visit_data).
#' @param patients Optional \code{data.table} of patient info with \code{STUDY_ID} and DOB.
#' @param study_id_column_name Name of the study ID column in the data (default \code{"STUDY_ID"}).
#' @param date_of_birth_column_name  Name of the DOB column in the input data/patients (default \code{"DOB"}).
#' @param verbose Logical, print progress messages (default \code{TRUE}).
#'
#' @return \code{data.table} = \code{primary_dt} with a standardized \code{DOB} column
#'         of type \code{Date}, merged from patients if necessary.
#'
#' @details
#' Cases handled:
#' \enumerate{
#'   \item If \code{primary_dt} already has DOB as Date:  trust and renamed to "DOB" if needed.
#'   \item Else, if \code{patients} supplied: DOB merged in by STUDY_ID, coerced to Date if needed.
#'   \item Else:  attempt to derive DOB directly from \code{primary_dt}, coercing if needed.
#' }
#'
#' @examples
#' library(data.table)
#' visits <- data.table(STUDY_ID = 1, VISIT_DATE = Sys.Date())
#' pats   <- data.table(STUDY_ID = 1, DOB = as.Date("2020-01-01"))
#' apply_dob(visits, pats)
#'
#' @export
apply_dob <- function(primary_dt,
                            patients = NULL,
                            study_id_column_name = "STUDY_ID",
                            date_of_birth_column_name = "DOB",
                            verbose = TRUE) {
  #primary checks
  stopifnot(data.table::is.data.table(primary_dt))
  # capture name of the input variable for clearer messages
  primary_name <- deparse(substitute(primary_dt))
  # Defensive copy
  dt <- data.table::copy(primary_dt)
  # Standardize STUDY_ID in primary
  if (!study_id_column_name %in% names(dt)) {    stop("Primary table must include study_id_column_name: ", study_id_column_name)}
  if (study_id_column_name != "STUDY_ID") {data.table::setnames(dt, study_id_column_name, "STUDY_ID")}
  
  # --- Case 1: DOB present in primary_dt and already a Date ---
  if (date_of_birth_column_name %in% names(dt) && inherits(dt[[date_of_birth_column_name]], "Date")) {
    if (date_of_birth_column_name != "DOB") data.table::setnames(dt, date_of_birth_column_name, "DOB")
    if (verbose) message("DOB present in ",primary_name," as Date; patients not accessed if supplied. No further processing needed")
    return(dt)
  }
  
  # --- Case 2: Use patients file if provided ---
  if (!is.null(patients)) {
    if (verbose) message("Using patients table to map DOB.")
    
    pats <- data.table::copy(patients)
    
    # Standardize STUDY_ID
    if (!study_id_column_name %in% names(pats)) stop("Patients must include study_id_column_name.")
    if (study_id_column_name != "STUDY_ID") data.table::setnames(pats, study_id_column_name, "STUDY_ID")
    
    # Standardize DOB
    if (!date_of_birth_column_name %in% names(pats)) stop("Patients must include a DOB column.")
    if (date_of_birth_column_name != "DOB") data.table::setnames(pats, date_of_birth_column_name, "DOB")
    
    # Normalize DOB type
    if (!inherits(pats$DOB, "Date")) {make_a_date(pats, "DOB")}
    
    # Merge
    data.table::setkey(pats, STUDY_ID)
    data.table::setkey(dt, STUDY_ID)
    dt <- pats[dt]
    if (verbose) message("DOB mapping from patients completed @", Sys.time())
    return(dt)
  }
  
  # --- Case 3: Derive DOB from primary_dt ---
  if (!date_of_birth_column_name %in% names(dt)) {stop("Need DOB in ",primary_name," if no patients dataset provided.")}
  if (verbose) message("Deriving DOB from ",primary_name," (no patients dataset). Least optimal solution.")
  if (date_of_birth_column_name != "DOB") data.table::setnames(dt, date_of_birth_column_name, "DOB")
  if (!inherits(dt$DOB, "Date")) make_a_date(dt, "DOB")
  #T/C If this is fast enough swap case 2 and case 3 if DOB is not a date on patients table.
  #return snake
  return(dt)
}

