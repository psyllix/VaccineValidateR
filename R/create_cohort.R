.datatable.aware = TRUE

#' Generate Cohort List with Age Calculations
#'
#' This function generates a cohort list from a dataset containing patient identifiers and 
#' dates of birth. It validates required columns, optionally renames them, ensures date consistency, 
#' and calculates age in both days and years.
#'
#' @param patients A \code{data.table} containing at least patient identifiers and date of birth.
#' @param date_of_birth_column_name String. Name of the DOB column in \code{patients}. Default is \code{"DOB"}.
#' @param study_id_column_name String. Name of the patient identifier column. Default is \code{"STUDY_ID"}.
#' @param reference_date Optional. A \code{Date} used for age calculation. Default is \code{Sys.Date()}.
#' @param verbose Logical. Whether to print progress messages. Default is \code{TRUE}.
#' 
#' @examples
#' library(data.table)
#'
#' # Example patient data
#' dt <- data.table(
#'   STUDY_ID = c(1, 2, 3),
#'   DOB = c("2020-01-15", "2019/07/30", "15-Mar-2018")
#' )
#'
#' # Create cohort with default reference date (today)
#' cohort <- create_cohort(dt)
#' print(cohort)
#'
#' # Specify a fixed reference date for reproducibility
#' cohort_fixed <- create_cohort(dt, reference_date = as.Date("2025-01-01"))
#' print(cohort_fixed)
#' @return A \code{data.table} with one row per patient, including:
#'   \itemize{
#'     \item \code{DOB}: earliest date of birth (as Date)
#'     \item \code{AGE_DAYS}: age in days at \code{reference_date}
#'     \item \code{AGE_YEARS}: age in years (decimal, using YEAR_LENGTH)
#'   }
#' @export

create_cohort<-function(patients,date_of_birth_column_name='DOB',study_id_column_name='STUDY_ID',reference_date=Sys.Date(),verbose=TRUE){
  patients<-data.table::setDT(data.table::copy(patients))
  if(verbose) message("Starting patient cohort list development. Data validation step...")
  #validate the study id column name
  if(!study_id_column_name %in% colnames(patients)){
    stop("Source data table requires a patient identifier column name to create a patient list. Default is STUDY_ID, otherwise please specify.")
  }else if(study_id_column_name!='STUDY_ID'){
    patients<-patients[,STUDY_ID:=get(study_id_column_name)]
  }
  #validate the DOB column name
  if(!date_of_birth_column_name %in% colnames(patients)){
      stop("Source data table requires a date of birth column name to create a patient list. Default is DOB, otherwise please specify")
  }else if(date_of_birth_column_name!='DOB'){
    patients<-patients[,DOB:=get(date_of_birth_column_name)]
  }
  # Ensure reference date is valid
  if (is.na(reference_date)) {
    reference_date <- Sys.Date()
  } else if (!inherits(reference_date, "Date")) {
    stop("reference_date must be a Date.")
  }
  
  
  # Coerce DOB to Date BEFORE aggregation
  if (!inherits(patients$DOB, "Date")) {
    # Try common formats; keep it dependency-free
    patients[, DOB := as.Date(
      DOB,
      tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d-%b-%Y", "%Y/%m/%d", "%m-%d-%Y")
    )]
  }
  
  # Warn & drop rows with unparseable DOB
  if (anyNA(patients$DOB)) {
    n_bad <- sum(is.na(patients$DOB))
    warning("Dropped ", n_bad, " row(s) with missing/unparseable DOB.")
    patients <- patients[!is.na(DOB)]
  }
  #Move onto creating the cohort
  if (verbose) message("Data validated. Generating cohort.")
  patients<-data.table::as.data.table(patients[,.(DOB=min(DOB)),by=STUDY_ID])
  patients[,AGE_DAYS:=as.integer(reference_date-DOB)]
  patients[,AGE_YEARS:=AGE_DAYS/YEAR_LENGTH]
  if (verbose) message("Cohort list development complete. N = ", nrow(patients))
  return(patients)
}