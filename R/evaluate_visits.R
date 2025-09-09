.datatable.aware = TRUE

#' This function takes a list of visits and evaluates whether each immunization that could have been given was given, and if not given, whether it was delayed on the antigen level.
#' 
#' @param visit_data data table with of visits.
#' @param antigens The third output from the validation function. Access via output$antigens. 
#' @param patients (Optional) Defaults to NA. Specifies a subset of patients for evaluation. If NA, all patients in visit_data are considered.
#' @param full_return Default is FALSE. If true, will return all visits-antigen combinations. If false, will not return combinations where it would be too soon, the patient has aged out, or the series is completed previously. 
#' @param antigen_to_eval Default is ALL. Replace with a list containing antigens of interest. See the Antigens help topic for a list of antigens that can be included.
#' @param study_id_column_name A string ('STUDY_ID' by default) that denotes the column name in visit_data containing study IDs.
#' @param date_of_birth_column_name A string ('DOB' by default) specifying the column name in visit_data for patients' dates of birth.
#' @param visit_date_column_name A string ('VISIT_DATE' by default) indicating the column name in visit_data where visit dates are recorded.
#' @param verbose turn logging on/off.
#' @return A table containing all visit-level antigen evaluations for DUE, MISSED, GIVEN and the status of prior and future doses (DELAYED/NOT DELAYED) as well as the date of the NEXT_DOSE.
#'
evaluate_visits<-function(visit_data
                          ,antigens
                          ,patients=NA
                          ,full_return=FALSE
                          ,antigen_to_eval=c('ALL')
                          ,study_id_column_name='STUDY_ID'
                          ,date_of_birth_column_name='DOB'
                          ,visit_date_column_name='VISIT_DATE'
                          ,verbose=TRUE
                          ){
  visit_data<-data.table::copy(visit_data)
  #EVALUATE INCOMING TABLES FOR COMPLETENESS ######
  if(verbose) message(paste0("Starting evaluation process for visits. Validating inputs."))
  return_columns<-c('STUDY_ID','VISIT_DATE','ANTIGEN','SYSTEM','VISIT_MODALITY','IS_VIRAL_SEASON','DOB','AGED_OUT','COMPLETED_PREVIOUSLY','NEXT_DOSE_MIN','NEXT_DOSE_RECOMMENDED','DUE','GIVEN','MISSED','DELAYED_PRIOR_DOSE','DELAYED_VISIT_DOSE','DELAYED_NEXT_DOSE','NEXT_GIVEN_DATE')
  return_columns_definitions <- c(
    "Unique identifier for the patient/subject",
    "Date of the clinical visit",
    "Vaccine antigen being evaluated",
    "System source of the record (if applicable)",
    "Type of visit (e.g., clinic, telehealth)",
    "Flag for viral season (Sept to Apr = TRUE)",
    "Date of birth of patient",
    "Flag if patient has exceeded eligible age for vaccine",
    "TRUE if the series was completed before this visit",
    "Earliest allowable date for the next dose",
    "Recommended date for the next dose",
    "TRUE if the vaccine dose is due at this visit",
    "TRUE if the vaccine dose was given at this visit",
    "TRUE if dose was due but not given at this visit",
    "TRUE if the previous dose was given late",
    "TRUE if the dose given at this visit was delayed",
    "TRUE if the next scheduled dose was given late",
    "Date of the next dose given (if available)")
  #VISIT STUDY ID
  if(!study_id_column_name %in% colnames(visit_data)){stop("Column for study ID is required in immunization_data. Either use STUDY_ID or declare the column name. If a Patients file is being used study id column name must match accross tables.")}
  if(study_id_column_name!='STUDY_ID'){
    visit_data<-visit_data[,STUDY_ID:=get(study_id_column_name)]
    if(!missing(patients) && !is.null(patients)){
      patients[, STUDY_ID := get(study_id_column_name)]
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
  #confrim the antigen data is validated
  if(!'data.table' %in% class(antigens)&& !'VALID' %in% colnames(antigens)){
    stop("Antigen data has not been validated.")
  }
  if(verbose) message(paste0("Validation check complete - formatting data. This may take a while..."))
  ####### INTERNAL HELPER FUNCITONS #####
  wk_no_grace <-function(x) {wk_with_grace(x,grace=0)}
  yr_no_grace <-function(x) {yr_with_grace(x,grace=0)}
  mon_no_grace <-function(x) {mon_with_grace(x,grace=0)}
  days_no_grace <-function(x) {days_with_grace(x,grace=0)}
  days_with_grace<-function(x,grace=4) {x-grace}
  wk_with_grace<-function(x,grace=4) {x*7-grace}
  mon_with_grace<-function(x,grace=4) {floor(365.25/12*x)-grace}
  yr_with_grace<-function(x,grace=4) {floor(365.25*x)-grace}
  min_na <- function(x,na.rm = TRUE) {if (length(x)>0) min(x,na.rm=na.rm) else NA}
  max_na <- function(x,na.rm = TRUE) {if (length(x)>0) max(x,na.rm=na.rm) else NA}
  #remove antigens we don't care about for determining due dates
  
  antigens_list<-c('POLIO','HIB','PCV','HEPA','HEPB','HPV','MCV','MENB','TETANUS','ROTA','COVID','RSV','MMR','VZV','INFLUENZA')
  if(!identical(antigen_to_eval, "ALL")){
      antigens_list<-antigen_to_eval
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
  if(!inherits(visit_data$DATE_GIVEN,'Date')){
    antigens$DATE_GIVEN <- as.Date(antigens$DATE_GIVEN,format="%Y-%m-%d")
  }
  visit_data[,VISIT_MONTH:=(lubridate::month(VISIT_DATE))]
  #create past and future dates for merging
  visit_data[,TOMORROW:=VISIT_DATE+1]
  visit_data[,YESTERDAY:=VISIT_DATE-1]
  visit_data[,IS_VIRAL_SEASON:=(VISIT_MONTH>=9|VISIT_MONTH<=4)]
  
  #prep the antigens for joining to visit_data
  antigens[,JOIN_DATE:=DATE_GIVEN]
  antigen_split <- split(antigens, by="ANTIGEN", keep.by=TRUE)
  #data.table::setkey(antigens,STUDY_ID,JOIN_DATE)
  if(verbose) message(paste0("Data formating completed, moving onto visit-antigen validation."))
  #for each antigen in the antigen list create a set of variables and apply that to the visit_data
  
  for(i in 1:length(antigens_list)){
    if(verbose) message(paste0("Evaluating status of visits for ",antigens_list[i]," vaccination status."))
    if(verbose) message(paste0("Start with yesterday - calculate next date an immunization is due and previous dose information."))
    
    #create new instance of visit_data for the antigen in question
    visit_antigen_eval<-visit_data
    this_antigen <- antigen_split[[antigens_list[i]]]
    visit_antigen_eval[,JOIN_DATE:=YESTERDAY]
    #data.table::setkey(visit_antigen_eval,STUDY_ID,JOIN_DATE)
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
    if(verbose) message(paste0("Determine if a dose was given at the day of the visit and if so, was it delayed."))
    #now only on the visit date grab if an antigen dose was given - roll is not need
    visit_antigen_eval[,JOIN_DATE:=VISIT_DATE]
    #data.table::setkey(visit_antigen_eval,STUDY_ID,JOIN_DATE)
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
   if(verbose) message(paste0("Add in content from the next possible antigen given."))
    
   visit_antigen_eval[,JOIN_DATE:=TOMORROW]
   #data.table::setkey(visit_antigen_eval,STUDY_ID,JOIN_DATE)
   visit_antigen_eval<-this_antigen[ANTIGEN==antigens_list[i],.(
     STUDY_ID
     ,JOIN_DATE
     ,DELAYED_NEXT_DOSE=DELAYED
     ,NEXT_DATE_GIVEN=DATE_GIVEN
   )][visit_antigen_eval, on = .(STUDY_ID,JOIN_DATE), roll = -Inf]
   
   #LOCAL VARIABLES 
   
   #DETERMINE IF IS DUE AT THE VISIT
   visit_antigen_eval[,IS_DUE:=(VISIT_DATE<AGE_OUT#CANNOT BE DUE IF TOO OLD
                            &!COMPLETED_PREVIOUSLY#CANNOT BE DUE IF COMPLETED PREVIOUSLY
                            &!is.na(NEXT_DOSE_RECOMMENDED)#MUST HAVE A NEXT RECOMMENDED DATE TO BE DUE - ALL VISITS SHOULD HAVE A RECOMMENDED NEXT DOSE IF NOT COMPLETED AND NOT AGED OUT
                            &(NEXT_DOSE_RECOMMENDED<=VISIT_DATE)#AND THE DATE RECOMMENDED IS BEFORE OR ON THE VISIT DATE
   )]
   #DETERMINE CONCEPTS FOR EACH ANTIGEN
   #BASED ON DAY BEFORE VISIT
   visit_antigen_eval[,ANTIGEN:=antigens_list[i]]
   visit_antigen_eval[,AGED_OUT:=!(VISIT_DATE<AGE_OUT)]
   #visit_antigen_eval[,LAST_GIVEN_DATE:=ifelse(GIVEN_AT_VISIT==TRUE,VISIT_DATE,DATE_GIVEN)]#DATE THE LAST IMMUNIZATION WAS GIVEN, NA IF no previous dose, INCLUDES THE DAY OF VISIT
   visit_antigen_eval[,COMPLETED_PREVIOUSLY:=COMPLETED_PREVIOUSLY]#FLAG THAT THE LAST DOSE GIVEN COMPLETED THE SERIES
   visit_antigen_eval[,DELAYED_PRIOR_DOSE:=DELAYED_LAST_DOSE]#FLAG THAT THE LAST DOSE WAS GIVEN DELAYED, NA if there was no doses before/during visit
   visit_antigen_eval[,DELAYED_VISIT_DOSE:=DELAYED_THIS_DOSE]#FLAG THAT THE LAST DOSE WAS GIVEN DELAYED, NA if there was no doses before/during visit
   visit_antigen_eval[,DELAYED_NEXT_DOSE:=DELAYED_NEXT_DOSE]#FLAG THAT THE NEXT DOSE WAS GIVEN DELAYED, NA if there was no future dose given
   visit_antigen_eval[,DUE:=IS_DUE]
   visit_antigen_eval[,GIVEN:=(GIVEN_AT_VISIT)]
   visit_antigen_eval[,MISSED:=(IS_DUE&!GIVEN_AT_VISIT)]
   visit_antigen_eval[,NEXT_GIVEN_DATE:=NEXT_DATE_GIVEN]#NEXT DATE A VACCINATION IS GIVEN, IF NA NO FUTURE DOSES
   #prep for return
   visit_return[[i]]<-visit_antigen_eval[,..return_columns]
  }
  #columns_to_remove<-c('TOMORROW','YESTERDAY')
  if(verbose) message(paste0("Evaluation completed for all visits. Preparing return object as visit-antigen relation table."))
  visit_return<-data.table::rbindlist(visit_return,fill=TRUE)
  if(verbose) message(paste0("Column in table: ",return_columns, ": ",return_columns_definitions))
  #return the relation table as EOF
  return(visit_return)
}

