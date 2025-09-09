.datatable.aware = TRUE
#' This function generates a cohort list from a dataset containing patient identifiers and 
#' dates of birth. It validates the presence of required columns, optionally renames them, 
#' ensures date consistency, and calculates age in both days and years. 
#' 
#' @param source_data data.table with at least a patient identifier and a date of birth column.
#' @param date_of_birth_column_name Name of the DOB column in 'source_data'. Default = "DOB".
#' @param study_id_column_name Name of the patient identifier column. Default = "STUDY_ID".
#' @param reference_date Optional Date to use for age calculation. Default = today's date.
#' @param verbose Logical. Whether to print progress messages. Default = TRUE.
#'
#' @return data.table with unique patients, DOB (as Date), AGE_DAYS, and AGE_YEARS
#' @export

create_cohort<-function(source_data,date_of_birth_column_name='DOB',study_id_column_name='STUDY_ID',reference_date=NA,verbose=TRUE){
  source_data<-data.table::copy(source_data)
  if(verbose) message("Starting patient cohort list development. Data validation step...")
  #validate the study id column name
  if(!study_id_column_name %in% colnames(source_data)){
    stop("Source data table requires a patient identifier column name to create a patient list. Default is STUDY_ID, otherwise please specify.")
  }else if(study_id_column_name!='STUDY_ID'){
    source_data<-source_data[,STUDY_ID:=get(study_id_column_name)]
  }
  #validate the DOB column name
  if(!date_of_birth_column_name %in% colnames(source_data)){
      stop("Source data table requires a data of birth column name to create a patient list. Default is DOB, otherwise please specify")
  }else if(date_of_birth_column_name!='DOB'){
    source_data<-source_data[,DOB:=get(date_of_birth_column_name)]
  }
  # Ensure reference date is valid
  if (is.na(reference_date)) {
    reference_date <- Sys.Date()
  } else if (!inherits(reference_date, "Date")) {
    stop("reference_date must be a Date.")
  }
  
  #create the cohort with DOB returned as a date and age in days/years calculated from the time point
  if (verbose) message("Data validated. Generating cohort.")
  p<-data.table::as.data.table(source_data[,.(DOB=min(DOB)),by=STUDY_ID])
  if(!inherits(p$DOB,"Date")){
    p$DOB <- as.Date(p$DOB,format="%Y-%m-%d")
  }
  p[,AGE_DAYS:=as.integer(reference_date-DOB)]
  p[,AGE_YEARS:=AGE_DAYS/365.25]
  if (verbose) message("Cohort list development complete. N = ", nrow(p))
  return(p)
}