.datatable.aware = TRUE

#' Evaluate Visits for Immunization Status
#'
#' This function evaluates a list of visits to determine whether each immunization that could have been given
#' was given, and if not, whether it was delayed, at the antigen level.
#'
#' @param visit_data A \code{data.table} of visit records. Must contain at least a study ID, DOB (or be linked
#'   via \code{patients}), and a visit date column.
#' @param antigens A \code{data.table} of antigen-level results, typically the \code{antigens} element from
#'   \code{validate_immunizations()}. Required. Must include \code{ANTIGEN}, \code{DATE_GIVEN}, \code{VALID},
#'   and series/delay fields (see \code{ANTIGEN_RETURN_COLUMNS}).
#' @param patients A \code{data.table} of patient info, with at least \code{STUDY_ID} and \code{DOB}.
#'   Required if DOB is not present in \code{visit_data}.
#' @param full_return Logical. Default is \code{FALSE}. If \code{TRUE}, returns all visit–antigen combinations.
#'   If \code{FALSE}, excludes rows where the patient has aged out, the series was completed previously,
#'   or the visit occurred before the earliest allowable dose date.
#' @param antigen_to_eval Character vector of antigens to evaluate. Default is \code{"ALL"}.
#' @param study_id_column_name Name of the column in \code{visit_data} containing study IDs.
#'   Default is \code{"STUDY_ID"}.
#' @param date_of_birth_column_name Name of the column in \code{visit_data} for dates of birth.
#'   Default is \code{"DOB"}.
#' @param visit_date_column_name Name of the column in \code{visit_data} for visit dates.
#'   Default is \code{"VISIT_DATE"}.
#' @param verbose Logical. If \code{TRUE}, prints messages during processing. Default is \code{TRUE}.
#'
#' @details
#' When \code{full_return = FALSE}, the output excludes visits where:
#' \itemize{
#'   \item The patient has aged out for that antigen
#'   \item The series was already completed before the visit
#'   \item The visit occurred before the minimum allowable dose date
#' }
#' This keeps only visit–antigen combinations that could meaningfully be evaluated.
#'
#' The returned \code{data.table} always contains the standardized set of columns defined in
#' \code{VISIT_RETURN_COLUMNS}, plus all columns originally present in \code{visit_data}.
#' Temporary internal columns (e.g., \code{JOIN_DATE}, \code{TOMORROW}, \code{YESTERDAY}) are dropped.
#'
#' @note Using \code{full_return = TRUE} can produce very large outputs and run slowly, since it
#' evaluates every possible visit–antigen pair regardless of clinical relevance. Use with caution
#' for large cohorts.
#'
#' @examples
#' library(data.table)
#'
#' # Minimal toy data
#' visits <- data.table(
#'   STUDY_ID = 1,
#'   VISIT_DATE = as.Date("2025-01-01"),
#'   DOB = as.Date("2024-01-01")
#' )
#' # antigens is an output from \code{\link{validate_immunizations}}
#'
#' # Evaluate visits (full return = FALSE trims unnecessary rows)
#' evaluate_visits(visits, antigens, patients = NULL, full_return = FALSE)
#'
#' # Keep all possible visit–antigen pairs
#' evaluate_visits(visits, antigens, patients = NULL, full_return = TRUE)
#'
#' @return A \code{data.table} containing:
#'   \itemize{
#'     \item All standardized columns in \code{VISIT_RETURN_COLUMNS}
#'     \item All original columns from \code{visit_data}
#'   }
#' Temporary evaluation fields are excluded.
#'
#' @seealso \code{\link{validate_immunizations}}, \code{VISIT_RETURN_COLUMNS}
#' @export

evaluate_visits<-function(visit_data
                          ,antigens#should only ever be the output from validate_immunizations
                          ,patients=NULL
                          ,full_return=FALSE
                          ,antigen_to_eval=c('ALL')
                          ,study_id_column_name='STUDY_ID'
                          ,date_of_birth_column_name='DOB'
                          ,visit_date_column_name='VISIT_DATE'
                          ,verbose=TRUE
                          ){
  # confirm the antigen data is validated
  check_antigen_table(antigens)
  #confirm uppercase
  antigen_to_eval <- toupper(antigen_to_eval)
  #create internal copy of immunization and patient data (if passed into function)
  visit_data<-data.table::setDT(data.table::copy(visit_data))
  #antigens and patients are only references, not changed - do not need to create a copy
  if (!is.null(patients)) {
    patients <- data.table::setDT(data.table::copy(patients))
  }
  
  #EVALUATE INCOMING TABLES FOR COMPLETENESS ######
  if(verbose) message(paste0("Starting evaluation process for visits. Validating inputs."))
 
  #VISIT STUDY ID validation
  if(!study_id_column_name %in% colnames(visit_data)){stop("Column for study ID is required in immunization_data. Either use STUDY_ID or declare the column name. If a Patients file is being used study id column name must match accross tables.")}
  if(study_id_column_name!='STUDY_ID'){
    visit_data<-visit_data[,STUDY_ID:=get(study_id_column_name)]
    if(!missing(patients) && !is.null(patients)){
      patients[, STUDY_ID := get(study_id_column_name)]
    }
  }
  # Ensure patient table has unique STUDY_ID
  if (!is.null(patients)) {
    dup_ids <- patients[, .N, by = STUDY_ID][N > 1, STUDY_ID]
    if (length(dup_ids) > 0) {
      stop("Patients table must have unique STUDY_ID values. Found duplicates for: ",
           paste(dup_ids, collapse = ", "))
    }
  }
  #VISIT DATE
  if(!visit_date_column_name %in% colnames(visit_data)){stop("Column for visit date is required in visit_data. Either use VISIT_DATE or declare the column name.")}
  if(visit_date_column_name!='VISIT_DATE'){visit_data<-visit_data[,VISIT_DATE:=get(visit_date_column_name)]}
  #DOB
  if(!date_of_birth_column_name %in% colnames(visit_data)){
    if(is.null(patients) || all(is.na(patients))){
      stop("Column for DOB is required in visit_data or in patients. Either use DOB or declare the column name. Patients table is only requierd if no DOB in visit_data")
    }
    else{
      if(!date_of_birth_column_name %in% colnames(patients)){
        stop("Column for DOB is required in visit_data or in patients. Either use DOB or declare the column name. Patients is only requierd if no DOB in visit_data")
      }
      if(date_of_birth_column_name!='DOB'){patients<-patients[,DOB:=get(date_of_birth_column_name)]}
      #now that DOB is confirmed, apply the DOB to visit table
      visit_data<-patients[,.(DOB=max(DOB,na.rm = TRUE))][visit_data, on=.(STUDY_ID)]
    }
  }
  if(date_of_birth_column_name!='DOB'&!'DOB' %in% colnames(visit_data)){visit_data<-visit_data[,DOB:=get(date_of_birth_column_name)]}
  if(verbose) message(paste0("validate_immunizations check complete - formatting data. This may take a while..."))
  
  #select the antigens we will evaluate against
  antigens_list<-SYSTEM_ANTIGENS
  if(!identical(antigen_to_eval, "ALL")){
    antigens_list<-intersect(antigens_list, antigen_to_eval)
    if (length(antigens_list) == 0L) {
      stop("No valid antigens in investigation list. Use all or see supported_antigens()")
    }
  }
  antigens<-antigens[ANTIGEN %in% antigens_list]
  #PREPARE RETURN storage
  visit_return<-vector("list",length(antigens_list))
  
  #confirm dates are in the correct format
  if(!inherits(visit_data$VISIT_DATE,'Date')){
    visit_data$VISIT_DATE <- as.Date(visit_data$VISIT_DATE,format="%Y-%m-%d")
  }
  if(!inherits(visit_data$DOB,'Date')){
    visit_data$DOB <- as.Date(visit_data$DOB,format="%Y-%m-%d")
  }
  if(!inherits(antigens$DATE_GIVEN,'Date')){
    antigens$DATE_GIVEN <- as.Date(antigens$DATE_GIVEN,format="%Y-%m-%d")
  }
  visit_data[,VISIT_MONTH:=(lubridate::month(VISIT_DATE))]
  #create past and future dates for merging
  visit_data[,TOMORROW:=VISIT_DATE+1]
  visit_data[,YESTERDAY:=VISIT_DATE-1]
  # Viral season = September (9) through April (4)
  visit_data[,IS_VIRAL_SEASON:=(VISIT_MONTH>=9|VISIT_MONTH<=4)]
  
  #prep the antigens for joining to visit_data
  antigens[,JOIN_DATE:=DATE_GIVEN]
  antigen_split <- split(antigens, by="ANTIGEN", keep.by=TRUE)
  if(verbose) message(paste0("Data formatting completed, moving onto visit-antigen validate_immunizations."))
  #for each antigen in the antigen list create a set of variables and apply that to the visit_data
  cols_to_keep <- union(VISIT_RETURN_COLUMNS, setdiff(names(visit_data), VISIT_RETURN_COLUMNS))
  
  for(i in 1:length(antigens_list)){
    if(verbose) message(paste0("Evaluating status of visits for ",antigens_list[i]," vaccination status."))
    # Pull this antigen’s subset
    visit_antigen_eval<-visit_data#need all visits each time - even aged out visits
    this_antigen <- antigen_split[[antigens_list[i]]]
    # Skip if missing or empty
    if (is.null(this_antigen) || nrow(this_antigen) == 0) {
      if (verbose) message("  No records found for ", antigens_list[i], ". Skipping evaluation. Check CVX mapping.")
      next
    }
    if(verbose) message(paste0("--Start with yesterday - calculate next date an immunization is due and previous dose information."))
    #create new instance of visit_data for the antigen in question
    visit_antigen_eval[,JOIN_DATE:=YESTERDAY]
    #evaluate the merge from previous day visits
    visit_antigen_eval<-this_antigen[ANTIGEN==antigens_list[i],.(
      STUDY_ID
     ,JOIN_DATE
     ,DATE_GIVEN
     ,DELAYED_LAST_DOSE=DELAYED
     ,NEXT_DOSE_MIN
     ,NEXT_DOSE_RECOMMENDED
     ,COMPLETED_PREVIOUSLY=DOSE_COMPLETES_SERIES#DID THE LAST DOSE COMPLETE THE SERIES
   )][visit_antigen_eval, on = .(STUDY_ID,JOIN_DATE), roll = Inf]
    visit_antigen_eval[is.na(COMPLETED_PREVIOUSLY),COMPLETED_PREVIOUSLY:=FALSE]
    if(verbose) message(paste0("--Determine if a dose was given at the day of the visit and if so, was it delayed."))
    #now only on the visit date grab if an antigen dose was given - roll is not need
    visit_antigen_eval[,JOIN_DATE:=VISIT_DATE]
    visit_antigen_eval<-this_antigen[ANTIGEN==antigens_list[i],.(
      STUDY_ID
      ,JOIN_DATE
      ,DELAYED_THIS_DOSE=DELAYED
      ,GIVEN_AT_VISIT=(DATE_GIVEN==JOIN_DATE)
      ,DELAYED_CURRENT_DOSE=DELAYED
    )][visit_antigen_eval,on = .(STUDY_ID,JOIN_DATE)]
  visit_antigen_eval[is.na(GIVEN_AT_VISIT),GIVEN_AT_VISIT:=FALSE]
   #APPLY AGE MINIMUNMS TO DATE MIN/REC times
   if(antigens_list[i] %in% c('HIB','PCV','POLIO','TETANUS','ROTA')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+wk_no_grace(6),NEXT_DOSE_RECOMMENDED=DOB+wk_no_grace(8))]
   }
    else if(antigens_list[i] %in% c('COVID')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(1),NEXT_DOSE_RECOMMENDED=DOB+mon_no_grace(6))]
    }
   else if(antigens_list[i] %in% c('VZV','MMR','HEPA')){
            visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(1),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(1))]
   }
   else if(antigens_list[i] %in% c('HEPB','RSV')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(0),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(0))]
   }
   else if(antigens_list[i] %in% c('HPV')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(9),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(9))]
   }
   else if(antigens_list[i] %in% c('MCV')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(10),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(11))]
   }
    else if(antigens_list[i] %in% c('MENB')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(16),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(16))]
    }
   else if(antigens_list[i] %in% c('INFLUENZA')){
     visit_antigen_eval[is.na(DATE_GIVEN)==TRUE,`:=`(NEXT_DOSE_MIN=pmax(DOB+mon_no_grace(6),as.Date(paste0(lubridate::year(DOB)+1,"-09-01")))
                                                                 ,NEXT_DOSE_RECOMMENDED=pmax(DOB+mon_no_grace(6),as.Date(paste0(lubridate::year(DOB)+1,"-09-01"))))]
   }
   #apply maximum ages
   if(antigens_list[i] %in% c('HIB','PCV')){
     visit_antigen_eval[,AGE_OUT:=DOB+yr_no_grace(5)]
   }
  else if(antigens_list[i] %in% c('RSV')){
    visit_antigen_eval[,AGE_OUT:=DOB+yr_no_grace(2)]
  }
   else if(antigens_list[i] %in% c('ROTA')){
      visit_antigen_eval[,AGE_OUT:=DOB+mon_no_grace(8)]
   }
   else
     visit_antigen_eval[,AGE_OUT:=DOB+yr_no_grace(20)]
   
  #IF NOT FULL APPLY THE RETURN LIMIT HERE (EARLIEST POSSIBLE POINT TO REMOVE VISITS)
  if(full_return==FALSE){
    visit_antigen_eval<-visit_antigen_eval[VISIT_DATE<AGE_OUT&VISIT_DATE>=NEXT_DOSE_MIN&!COMPLETED_PREVIOUSLY]
  }
  
  #JOIN TOMORROW
   if(verbose) message(paste0("--Add in content from the next possible antigen given."))
    
   visit_antigen_eval[,JOIN_DATE:=TOMORROW]
   #data.table::setkey(visit_antigen_eval,STUDY_ID,JOIN_DATE)
   visit_antigen_eval<-this_antigen[ANTIGEN==antigens_list[i],.(
     STUDY_ID
     ,JOIN_DATE
     ,DELAYED_NEXT_DOSE=DELAYED
     ,NEXT_DATE_GIVEN=DATE_GIVEN
   )][visit_antigen_eval, on = .(STUDY_ID,JOIN_DATE), roll = -Inf]
   
   
   #DETERMINE IF IS DUE AT THE VISIT
   visit_antigen_eval[,IS_DUE:=(VISIT_DATE<AGE_OUT#CANNOT BE DUE IF TOO OLD
                            &!COMPLETED_PREVIOUSLY#CANNOT BE DUE IF COMPLETED PREVIOUSLY
                            &!is.na(NEXT_DOSE_RECOMMENDED)#MUST HAVE A NEXT RECOMMENDED DATE TO BE DUE - ALL VISITS SHOULD HAVE A RECOMMENDED NEXT DOSE IF NOT COMPLETED AND NOT AGED OUT
                            &(NEXT_DOSE_RECOMMENDED<=VISIT_DATE)#AND THE DATE RECOMMENDED IS BEFORE OR ON THE VISIT DATE
   )]
   #DETERMINE CONCEPTS FOR EACH ANTIGEN
   visit_antigen_eval[,ANTIGEN:=antigens_list[i]]
   visit_antigen_eval[,AGED_OUT:=!(VISIT_DATE<AGE_OUT)]#SHOULD ALWAYS BE FALSE IF NOT FULL RETURN
   #visit_antigen_eval[,LAST_GIVEN_DATE:=ifelse(GIVEN_AT_VISIT==TRUE,VISIT_DATE,DATE_GIVEN)]#DATE THE LAST IMMUNIZATION WAS GIVEN, NA IF no previous dose, INCLUDES THE DAY OF VISIT
   visit_antigen_eval[,DELAYED_PRIOR_DOSE:=DELAYED_LAST_DOSE]#FLAG THAT THE LAST DOSE WAS GIVEN DELAYED, NA if there was no doses before/during visit
   visit_antigen_eval[,DELAYED_VISIT_DOSE:=DELAYED_THIS_DOSE]#FLAG THAT THE LAST DOSE WAS GIVEN DELAYED, NA if there was no doses before/during visit
   visit_antigen_eval[,DELAYED_NEXT_DOSE:=DELAYED_NEXT_DOSE]#FLAG THAT THE NEXT DOSE WAS GIVEN DELAYED, NA if there was no future dose given
   visit_antigen_eval[,DUE:=IS_DUE]#RE CASTING TO DECREASE CONFUSION AS IS_DUE EVOKES NOT GIVEN VIBES
   visit_antigen_eval[,GIVEN:=(GIVEN_AT_VISIT)]
   visit_antigen_eval[,MISSED:=(IS_DUE&!GIVEN_AT_VISIT)]#DUE and !GIVEN
   visit_antigen_eval[,NEXT_GIVEN_DATE:=NEXT_DATE_GIVEN]#NEXT DATE A VACCINATION IS GIVEN, IF NA NO FUTURE DOSES
   visit_antigen_eval[,NEXT_RCVD_15_DAYS:=(!is.na(NEXT_GIVEN_DATE)&NEXT_GIVEN_DATE<=VISIT_DATE+15)]
   visit_antigen_eval[,NEXT_RCVD_30_DAYS:=(!is.na(NEXT_GIVEN_DATE)&NEXT_GIVEN_DATE<=VISIT_DATE+30)]
   visit_antigen_eval[,NEXT_RCVD_90_DAYS:=(!is.na(NEXT_GIVEN_DATE)&NEXT_GIVEN_DATE<=VISIT_DATE+90)]
   #prep for return - keep any extra columns in original data set
   visit_return[[i]] <- visit_antigen_eval[, ..cols_to_keep]
   
  }
  if(verbose) message(paste0("Evaluation completed for all visits. Preparing return object as visit-antigen relation table."))
  #create the visit-antigen relation table for return
  visit_return<-data.table::rbindlist(visit_return,fill=TRUE)
  
  #Ensure all expected return columns exist
  for (col in VISIT_RETURN_COLUMNS) {
    if (!col %in% names(visit_return)) {
      visit_return[, (col) := NA]
      warning("Added missing column: ", col, " (filled with NA).")
    }
  }
  data.table::setattr(visit_return, "processed", TRUE)
  if(verbose) message("See ?VISIT_RETURN_DEFS for explanations of return columns.")
  #return the relation table
  return(visit_return)
}

