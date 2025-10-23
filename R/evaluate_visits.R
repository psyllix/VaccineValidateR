.datatable.aware <- TRUE

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
#' @param return_mode Character; controls which visit–antigen rows are returned.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"FULL"} — return *all* evaluated visit–antigen combinations,
#'           including those where the patient was too young, aged out, or the series was completed.
#'     \item \code{"POSSIBLE"} — return only visit–antigen combinations that were
#'           at least theoretically possible (after minimum age and interval),
#'           but exclude those where the antigen was no longer recommended or series already complete.
#'     \item \code{"DUE"} — return only visit–antigen combinations where the immunization
#'           would be recommended at that visit.
#'   }
#'   Default is \code{"POSSIBLE"}.
#' @param antigen_to_eval Character vector of antigens to evaluate. Default is \code{"ALL"}.
#' @param study_id_column_name Name of the column in \code{visit_data} containing study IDs.
#'   Default is \code{"STUDY_ID"}.
#' @param date_of_birth_column_name Name of the column in \code{visit_data} for dates of birth.
#'   Default is \code{"DOB"}.
#' @param visit_date_column_name Name of the column in \code{visit_data} for visit dates.
#'   Default is \code{"VISIT_DATE"}.
#' @param years A numeric vector of evaluation years to include test.
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
                          ,return_mode=c("POSSIBLE","FULL","DUE")#DELAYED
                          ,antigen_to_eval=c('ALL')
                          ,study_id_column_name='STUDY_ID'
                          ,date_of_birth_column_name='DOB'
                          ,visit_date_column_name='VISIT_DATE'
                          ,years = NULL
                          ,verbose=TRUE
                          ){
  # confirm the antigen data is validated
  check_antigen_table(antigens)
  #confirm uppercase
  antigen_to_eval <- toupper(antigen_to_eval)
  #create internal copy of immunization and patient data (if passed into function)
  visit_data<-data.table::setDT(data.table::copy(visit_data))
  #EVALUATE INCOMING TABLES FOR COMPLETENESS ######
  if(verbose) message(paste0("Starting evaluation process for visits. Validating inputs."))
 
  #VISIT STUDY ID validation
  if(!study_id_column_name %in% colnames(visit_data)){stop("Column for study ID is required in immunization_data. Either use STUDY_ID or declare the column name. If a Patients file is being used study id column name must match accross tables.")}
  if(study_id_column_name!='STUDY_ID'){visit_data<-visit_data[,STUDY_ID:=get(study_id_column_name)]}
  #APPLY DOB to VISIT table
  visit_data<-apply_dob(visit_data
                        ,patients=patients
                        ,study_id_column_name=study_id_column_name
                        ,date_of_birth_column_name=date_of_birth_column_name)
  
  #select the antigens we will evaluate against - confirm that only allowed antigens exist
  antigens_list<-SYSTEM_ANTIGENS
  if(!identical(antigen_to_eval, "ALL")){
    antigens_list<-intersect(antigens_list, antigen_to_eval)
    if (length(antigens_list) == 0L) {
      stop("No valid antigens in investigation list. Use all or see supported_antigens()")
    }
  }
  antigens<-antigens[ANTIGEN %in% antigens_list]
  #PREPARE RETURN storage
  visits_evaluated<-vector("list",length(antigens_list))
  
  #confirm dates are in the correct format
  if(!inherits(visit_data$VISIT_DATE,'Date')){
    make_a_date(visit_data,"VISIT_DATE")
  }
  # Optional year filter
  if (!is.null(years)) {
    if (verbose) message("Restricting visits to years: ", paste(years, collapse = ", "))
    visit_data <- visit_data[year_from_date(VISIT_DATE) %in% years]
    if (nrow(visit_data) == 0L) {
      warning("No visits remain after year filter. Returning empty table.")
      return(data.table())
    }
  }  
  
  if(!inherits(visit_data$DOB,'Date')){
    if(verbose) message("DATE mapping needed:DOB")
    make_a_date(visit_data,"DOB")
  }
  #date given should already exist on antigens, but good to check
  if(!inherits(antigens$DATE_GIVEN,'Date')){
    if(verbose) message("DATE mapping needed:DATE_GIVEN")
    make_a_date(antigens,"DATE_GIVEN")
  }
  if(verbose) message("Extract and determine dates of interest.")
  visit_data[,VISIT_MONTH:=(month_from_date(VISIT_DATE))]
  #create past and future dates for merging
  visit_data[,TOMORROW:=VISIT_DATE+1]
  visit_data[,YESTERDAY:=VISIT_DATE-1]
  # Viral season = September (9) through April (4)
  visit_data[,IS_VIRAL_SEASON:=(VISIT_MONTH>=9|VISIT_MONTH<=4)]
  
  #prep the antigens for joining to visit_data
  if(verbose) message("Prepare the antigens for joining. Antigen groups")
  antigens[,JOIN_DATE:=DATE_GIVEN]
  data.table::setkey(antigens, ANTIGEN, STUDY_ID, DATE_GIVEN)
  #antigen_split <- split(antigens, by="ANTIGEN", keep.by=TRUE)
  if(verbose) message(paste0("Data formatting completed, moving onto visit-antigen validate_immunizations."))
  #for each antigen in the antigen list create a set of variables and apply that to the visit_data
  cols_to_keep <- union(VISIT_RETURN_COLUMNS, setdiff(names(visit_data), VISIT_RETURN_COLUMNS))
  
  for(ant in antigens_list){
    if(verbose) message(paste0("Evaluating status of visits for ",ant," vaccination status. Starting ",Sys.time()))
    # Pull this antigen’s subset
    visit_antigen_eval <- data.table::copy(visit_data) #need all visits each time - even aged out visits
    #this_antigen <- antigen_split[[ant]]
    this_antigen <- antigens[.(ant)]
    # Skip if missing or empty
    if (is.null(this_antigen) || nrow(this_antigen) == 0) {
      if (verbose) message("  No records found for ", ant, ". Skipping evaluation. Check CVX mapping.")
      next
    }
    if(verbose) message(paste0("--Start with yesterday - calculate next date an immunization is due and previous dose information."))
    #create new instance of visit_data for the antigen in question
    visit_antigen_eval[,JOIN_DATE:=YESTERDAY]
    #evaluate the merge from previous day visits
    visit_antigen_eval<-this_antigen[ANTIGEN==ant,.(
      STUDY_ID
     ,JOIN_DATE
     ,DATE_GIVEN
     ,LAST_GIVEN=DATE_GIVEN
     ,DELAYED_LAST_DOSE=DELAYED
     ,NEXT_DOSE_MIN
     ,NEXT_DOSE_RECOMMENDED
     ,COMPLETED_PREVIOUSLY=DOSE_COMPLETES_SERIES#DID THE LAST DOSE COMPLETE THE SERIES
   )][visit_antigen_eval, on = .(STUDY_ID,JOIN_DATE), roll = Inf]
    visit_antigen_eval[is.na(COMPLETED_PREVIOUSLY),COMPLETED_PREVIOUSLY:=FALSE]
    if(verbose) message(paste0("--Determine if a dose was given at the day of the visit and if so, was it delayed."))
    #now only on the visit date grab if an antigen dose was given - roll is not need
    visit_antigen_eval[,JOIN_DATE:=VISIT_DATE]
    visit_antigen_eval<-this_antigen[ANTIGEN==ant,.(
      STUDY_ID
      ,JOIN_DATE
      ,DELAYED_THIS_DOSE=DELAYED
      ,GIVEN_AT_VISIT=(DATE_GIVEN==JOIN_DATE)
      ,DELAYED_CURRENT_DOSE=DELAYED
    )][visit_antigen_eval,on = .(STUDY_ID,JOIN_DATE)]
  visit_antigen_eval[is.na(GIVEN_AT_VISIT),GIVEN_AT_VISIT:=FALSE]
   #APPLY AGE MINIMUNMS TO DATE MIN/REC times
   if(ant %in% c('HIB','PCV','POLIO','TETANUS','ROTA')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+wk_no_grace(6),NEXT_DOSE_RECOMMENDED=DOB+wk_no_grace(8))]
   }
    else if(ant %in% c('COVID')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(1),NEXT_DOSE_RECOMMENDED=DOB+mon_no_grace(6))]
    }
   else if(ant %in% c('VZV','MMR','HEPA')){
            visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(1),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(1))]
   }
   else if(ant %in% c('HEPB','RSV')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(0),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(0))]
   }
   else if(ant %in% c('HPV')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(9),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(9))]
   }
   else if(ant %in% c('MCV')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(10),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(11))]
   }
    else if(ant %in% c('MENB')){
      visit_antigen_eval[is.na(DATE_GIVEN),`:=`(NEXT_DOSE_MIN=DOB+yr_no_grace(16),NEXT_DOSE_RECOMMENDED=DOB+yr_no_grace(16))]
    }
   else if(ant %in% c('INFLUENZA')){
     visit_antigen_eval<-visit_antigen_eval[IS_VIRAL_SEASON==TRUE]#remove non viral season
     visit_antigen_eval[is.na(DATE_GIVEN)==TRUE,NEXT_DOSE_MIN:=pmax(DOB+mon_no_grace(6),next_sept1(DOB))]
     visit_antigen_eval[is.na(DATE_GIVEN)==TRUE,NEXT_DOSE_RECOMMENDED:=NEXT_DOSE_MIN]
   }
   #apply maximum ages
   if(ant %in% c('HIB','PCV')){
     visit_antigen_eval[,AGE_OUT:=DOB+yr_no_grace(5)]
   }
  else if(ant %in% c('RSV')){
    visit_antigen_eval[,AGE_OUT:=DOB+yr_no_grace(2)]
  }
   else if(ant %in% c('ROTA')){
      visit_antigen_eval[,AGE_OUT:=DOB+mon_no_grace(8)]
   }
   else
     visit_antigen_eval[,AGE_OUT:=DOB+yr_no_grace(20)]
  
  #DETERMINE IF IS DUE AT THE VISIT
  
  visit_antigen_eval[,IS_DUE:=(VISIT_DATE<AGE_OUT#CANNOT BE DUE IF TOO OLD
                               &!COMPLETED_PREVIOUSLY#CANNOT BE DUE IF COMPLETED PREVIOUSLY
                               &!is.na(NEXT_DOSE_RECOMMENDED)#MUST HAVE A NEXT RECOMMENDED DATE TO BE DUE - ALL VISITS SHOULD HAVE A RECOMMENDED NEXT DOSE IF NOT COMPLETED AND NOT AGED OUT
                               &(NEXT_DOSE_RECOMMENDED<=VISIT_DATE)#AND THE DATE RECOMMENDED IS BEFORE OR ON THE VISIT DATE
  )]
  visit_antigen_eval[is.na(IS_DUE),IS_DUE:=FALSE]
  # Apply return_mode filters ------------------------------------------
  if (return_mode == "POSSIBLE") {
    if (verbose) message("Returning visits where immunization could have been given (after min date, before age-out, not completed).")
    visit_antigen_eval <- visit_antigen_eval[
      VISIT_DATE < AGE_OUT &
        VISIT_DATE >= NEXT_DOSE_MIN &
        !COMPLETED_PREVIOUSLY
    ]
    
  } else if (return_mode=="DUE") {
    if (verbose) message("Returning visits where vaccine was due at time of visit.")
    visit_antigen_eval <- visit_antigen_eval[IS_DUE == TRUE]
  }
  #JOIN TOMORROW
   if(verbose) message(paste0("--Add in content from the next possible antigen given."))
    
   visit_antigen_eval[,JOIN_DATE:=TOMORROW]
   #data.table::setkey(visit_antigen_eval,STUDY_ID,JOIN_DATE)
   visit_antigen_eval<-this_antigen[ANTIGEN==ant,.(
     STUDY_ID
     ,JOIN_DATE
     ,DELAYED_NEXT_DOSE=DELAYED
     ,NEXT_DATE_GIVEN=DATE_GIVEN
   )][visit_antigen_eval, on = .(STUDY_ID,JOIN_DATE), roll = -Inf]
   #JOIN TOMORROW
   #DETERMINE CONCEPTS FOR EACH ANTIGEN
   visit_antigen_eval[,ANTIGEN:=ant]
   visit_antigen_eval[,AGED_OUT:=!(VISIT_DATE<AGE_OUT)]#SHOULD ALWAYS BE FALSE IF NOT FULL RETURN
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
   visits_evaluated[[ant]] <- visit_antigen_eval[, ..cols_to_keep]
   
  }
  if(verbose) message(paste0("Evaluation completed for all visits. Preparing return object as visit-antigen relation table."))
  #create the visit-antigen relation table for return
  visits_evaluated<-data.table::rbindlist(visits_evaluated,fill=TRUE)
  
  #Ensure all expected return columns exist
  for (col in VISIT_RETURN_COLUMNS) {
    if (!col %in% names(visits_evaluated)) {
      visits_evaluated[, (col) := NA]
      warning("Added missing column: ", col, " (filled with NA).")
    }
  }
  #add the processed attribute to the table to allow next step to proceed
  data.table::setattr(visits_evaluated, "processed", TRUE)
  if(verbose) message("See ?visits_evaluated_DEFS for explanations of return columns.")
  #return the relation table
  return(visits_evaluated)
}

