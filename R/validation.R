.datatable.aware = TRUE

#' Evaluate Immunizations Against ACIP Schedule
#'
#' This function takes a list of immunizations and evaluates them against minimum intervals based on the ACIP schedule, determining the immunity status across a list of antigens.
#' If any live virus products are included, all live-virus-containing products are evaluated.
#' Either a patient list or dates of birth within a patient list are required. A patient list will not be returned with the immunization data.
#' Either a date the immunization was given or age at immunization in days is required. Using age in days is notably faster for large datasets.
#' Either CVX codes or a map between local immunization antigen data and CVX codes is required.
#' Default column names can be overwritten using the individual parameters. Please use the same column names across files!
#'
#' @param immunization_data A \code{data.table} of immunizations containing at minimum \code{STUDY_ID}, \code{PRODUCT}, and either \code{DATE_GIVEN} or \code{AGE_IMM_GIVEN}. 
#'   If there is no patient table, a \code{DOB} column is required. \code{CVX}.
#' @param patients Optional \code{data.table} of patient info, with DOB. 
#'   Default is \code{NULL}, in which case DOB must be present in \code{visit_data}.
#' @param lim_cvx_map Optional \code{data.table} mapping local immunization identifiers (\code{LIM}) to CVX codes. Required if \code{CVX} is not present in the immunization data. Must include \code{PRODUCT} and \code{CVX}.
#'   Default is \code{NULL}, in which case the system will assume that immunization data has a valid CVX' 
#' @param antigens_to_eval Character vector of antigens to evaluate. Default is \code{"ALL"}.
#' @param verbose Logical. Default is \code{TRUE}. If \code{TRUE}, console messages are printed after each step.
#' @param immunization_code_column_name Column name containing CVX codes. Default is \code{"CVX"}.
#' @param local_immunization_identifier_column_name Column name containing local immunization identifiers. Default is \code{"LIM"}.
#' @param immunization_product_name_column_name Column name containing product names. Default is \code{"PRODUCT"}.
#' @param immunization_date_given_column_name Column name containing the date the immunization was given. Default is \code{"DATE_GIVEN"}.
#' @param age_at_immunization_column_name Column name containing age at immunization in days. Default is \code{"AGE_IMM_GIVEN"}.
#' @param reference_date Date cutoff for immunizations. Default is \code{Sys.Date()}.
#'
#' @return A named list with three \code{data.table} elements:
#' \describe{
#'   \item{immunizations}{Product-level immunization data after processing. 
#'     Typical columns include:
#'     \itemize{
#'       \item \code{STUDY_ID} – patient identifier
#'       \item \code{PRODUCT} – vaccine product name
#'       \item \code{CVX} – vaccine CVX code (integer)
#'       \item \code{DATE_GIVEN} – immunization administration date
#'       \item \code{AGE_IMM_GIVEN} – age in days at administration
#'       \item Antigen flags (e.g., \code{POLIO}, \code{MMR}, \code{HIB}) indicating antigen mapping
#'     }}
#'   \item{antigens}{Antigen-level dose data including validation checks.
#'     Typical columns include:
#'     \itemize{
#'       \item \code{STUDY_ID}, \code{ANTIGEN}, \code{DOB}, \code{PRODUCT}, \code{CVX}
#'       \item \code{DATE_GIVEN}, \code{AGE_IMM_GIVEN}
#'       \item \code{DOSE_COUNTER}, \code{ABS_ADMIN_COUNTER}, \code{TABLE_INDEX}
#'       \item Validation fields: \code{VALID}, \code{SERIES_COMPLETE}, \code{DELAYED}, \code{NOTES}
#'       \item Interval tracking: \code{INTERVAL}, \code{LIVE_INTERVAL}, \code{CYCLE}, \code{NEXT_DOSE_MIN}, \code{NEXT_DOSE_RECOMMENDED}
#'     }}
#'   \item{skipped_antigens}{Immunizations that were excluded from evaluation due to failing rules
#'     (e.g., given before birth, given too early, live vaccine spacing violation). 
#'     Contains the same structure as \code{antigens} plus \code{NOTES} indicating the reason.}
#' }
#' @export
  validate_immunizations <- function(immunization_data
                         ,patients=NULL
                         ,lim_cvx_map
                         ,antigen_to_eval=c('ALL')
                         ,verbose = TRUE
                         ,date_of_birth_column_name='DOB'
                         ,study_id_column_name='STUDY_ID'
                         ,immunization_code_column_name='CVX'
                         ,local_immunization_identifier_column_name='LIM'
                         ,immunization_product_name_column_name='PRODUCT'
                         ,immunization_date_given_column_name='DATE_GIVEN'
                         ,age_at_immunization_column_name='AGE_IMM_GIVEN'
                         ,reference_date=Sys.Date()
                         ){
      #Initialization
      
      if(verbose){message(paste0("Starting validatation of: ",immunization_data[,.N]," immunization @",Sys.time()))}
      immunization_data<-data.table::setDT(data.table::copy(immunization_data))
      #pull attributes immediately to get submitted values
      processed_attr <- attr(immunization_data, "processed")
      antigens_used  <- attr(immunization_data, "antigens_used")
      antigen_to_eval <- toupper(antigen_to_eval)
      
      #Validation
      if(verbose) message("validation of immunizations may take some time depending on your data and system.")
      if(verbose){message("Please wait.... Status updates will occur routinely if verbose is True.")}
      
      #Check for study id 
      if(!study_id_column_name %in% colnames(immunization_data)){stop("Column for study ID is required in immunization_data. Either use STUDY_ID or declare the column name. Must be the same accross files.")}
      if(study_id_column_name!='STUDY_ID'){data.table::setnames(immunization_data,study_id_column_name,"STUDY_ID")}
      
      #Check for Product
      if(!immunization_product_name_column_name %in% colnames(immunization_data)){stop("Column for product name is required. Product names are used for disambinuation. Either use PRODUCT or declare the product column name.")}
      if(immunization_product_name_column_name!='PRODUCT'){
        data.table::setnames(immunization_data,immunization_product_name_column_name,"PRODUCT")
      }
  
      #Product Codes/LIM validate_immunizations
      if (!immunization_code_column_name %in% names(immunization_data)) {
        if (is.null(lim_cvx_map)) {
          stop("Need CVX column in immunization_data or lim_cvx_map with ",
               local_immunization_identifier_column_name, " and ", immunization_code_column_name)
        }
        if (!(local_immunization_identifier_column_name %in% names(lim_cvx_map) &&
              immunization_code_column_name %in% names(lim_cvx_map))) {
          stop("lim_cvx_map must include columns: ",
               local_immunization_identifier_column_name, " and ", immunization_code_column_name)
        }
      }
      # ----------------------------
      # Handle patient DOB in wither imm_data or patient file
      # ----------------------------
      
      if (date_of_birth_column_name %in% names(immunization_data) && inherits(immunization_data[[date_of_birth_column_name]], "Date")) {
        # Case 1: immunization_data has DOB column already as Date -> trust it
        if (date_of_birth_column_name != "DOB") {data.table::setnames(immunization_data, date_of_birth_column_name, "DOB")}
        if (verbose) message("DOB present in immunization_data as Date; patients file not accessed if supplied.")
        
      } else if (!is.null(patients)) {
        # Case 2: use patients data file to create the date on immunization - it may or may not have a Date format
        if (verbose) message("Using patient dataset for DOB mapping.")
        # Ensure STUDY_ID
        if (!study_id_column_name %in% names(patients)) {stop("Patient data must include a study ID column.")}
        if (study_id_column_name != "STUDY_ID") {data.table::setnames(patients, study_id_column_name, "STUDY_ID")}
        # Ensure DOB column
        if (!date_of_birth_column_name %in% names(patients)) {stop("Patient data must include a DOB column.")}
        if (date_of_birth_column_name != "DOB") {data.table::setnames(patients, date_of_birth_column_name, "DOB")}
        # Normalize DOB
        if (!inherits(patients$DOB, "Date")) {make_a_date(patients,"DOB")} 
        # Merge
        data.table::setkey(patients, STUDY_ID)
        data.table::setkey(immunization_data, STUDY_ID)
        immunization_data <- patients[immunization_data]
        if(verbose){message(paste0("DOB mapping completed @",Sys.time()))}
      } else {
        # Case 3: derive patients from immunization_data
        if (!date_of_birth_column_name %in% names(immunization_data)) {stop("Need DOB in immunization_data if no separate patient dataset.")}#no reason to progress
        
        if (verbose) message("Deriving DOB from immunization_data (no patient dataset provided) - sub optimal path.")
        #make patients - never use again
        patients <- unique(immunization_data[, .(
          STUDY_ID = get(study_id_column_name),
          DOB  = get(date_of_birth_column_name)
        )])
        #no need to check since if DOB in immunizations was date would not get here
        make_a_date(patients,"DOB")
        data.table::setkey(patients, STUDY_ID)
        data.table::setkey(immunization_data, STUDY_ID)
        immunization_data <- patients[immunization_data]
        if(verbose){message(paste0("DOB mapping completed @",Sys.time()))}
      }
      
      if(verbose){message(paste0("Incoming data check passed and patient list exists: ",Sys.time()))}
      if(verbose){message(paste0("Starting immunization pre-processing... This may take a few minutes..."))}
      
      ####LIMIT DATA IF DATA FILE ALREADY PROCESSED#####
      if (is.null(processed_attr) ||!isTRUE(processed_attr)) {
        immunization_data<-select_immunizations(
          immunization_data,
          ,antigen_to_eval=antigen_to_eval
          ,verbose=FALSE
          ,immunization_date_given_column_name=immunization_date_given_column_name
          ,age_at_immunization_column_name=age_at_immunization_column_name
          ,immunization_code_column_name=immunization_code_column_name
          ,local_immunization_identifier_column_name=local_immunization_identifier_column_name)
        if(verbose){message(paste0("Preprocessing completed."))}
      }
      else{if (verbose) message("Skipping preprocessing — immunization_data already marked as processed.")}
      # Immunization date processing - preference for Age modifications
      if(verbose){message(paste0("Immunization Date processing starting at @ ",Sys.time()))}
      
      #now create or format the DATE_GIVEN field - FIRST TIME INTENSIVE PROCESS - PREFER THAT WE GET AGE_IMM_GIVEN, 
      if(age_at_immunization_column_name %in% colnames(immunization_data)){
        if(verbose){message("Optimal path - AGE + patient DOB")}
        immunization_data[,DATE_GIVEN:=DOB+get(age_at_immunization_column_name)]#automatically gets us dates - no reformatting
        if (age_at_immunization_column_name != "AGE_IMM_GIVEN"){ 
          data.table::setnames(immunization_data,age_at_immunization_column_name,"AGE_IMM_GIVEN")
        }
      } else{
          if(!immunization_date_given_column_name %in% colnames(immunization_data)){stop("Column for age immunization given or date immunization given required in immunization data. Age immunization given is preferred.")}
          if(immunization_date_given_column_name!='DATE_GIVEN'){data.table::setnames(immunization_data,immunization_date_given_column_name,"DATE_GIVEN")}
          #check if date and make date if not
          if(!inherits(immunization_data$DATE_GIVEN,'Date')){make_a_date(immunization_data,"DATE_GIVEN")}
          if(verbose){message("suboptimal path - Date given - DOB")}
          immunization_data[,AGE_IMM_GIVEN:=as.integer(DATE_GIVEN-DOB)]
        }
    gc()
    
    if(verbose){message("Starting immunization data processing at ",Sys.time(),". This may take a few minutes...")}
    
    #### PART 2 #####
    # CONVERT to ANTIGENS from immunization_data and drop tables
    antigens_list <- antigens_used#from the select_immuniations table - at this point is completed
    
    antigens_v <- lapply(antigens_list, function(ant) {
      dt <- immunization_data[get(ant) == TRUE,
                              .(STUDY_ID, ANTIGEN = ant, DOB, PRODUCT, CVX,
                                DATE_GIVEN, AGE_IMM_GIVEN)
      ]
      if (nrow(dt) > 0) dt else NULL
    })
    
    antigens <- data.table::rbindlist(antigens_v, use.names = TRUE, fill = TRUE)
    rm(antigens_v)
    # add columns to Antigens table
    antigens[, `:=`(
      IS_LIVE        = CVX %in% cvx$LIVE_NON_ENTERAL,
      NOTES          = as.character(NA),
      VALID          = NA,
      DOSE_COUNTER   = as.numeric(0),
      IS_HIB4        = CVX %in% cvx$HIB4,
      IS_ROTA2       = CVX %in% cvx$ROTA2,
      IS_BEXSERO     = CVX %in% cvx$BEXSERO,
      IS_TRUMENBA    = CVX %in% cvx$TRUMENBA,
      IS_ADULT_HEPA  = CVX %in% cvx$ADULT_HEPA,
      SERIES_COMPLETE      = FALSE,
      DELAYED             = FALSE,
      COUNTER             = as.numeric(0),
      CYCLE               = as.numeric(0),
      NEXT_DOSE_MIN       = as.Date(NA),
      NEXT_DOSE_RECOMMENDED = as.Date(NA),
      DOSE_COMPLETES_SERIES = FALSE,
      LIVE_INTERVAL       = 0,
      LAST_LIVE           = as.numeric(0)
    )]
    ##### DEVELOP ANTIGEN LIST FROM IMMUNIZAITIONS  ####
    if(verbose){message(paste0("Antigen list developed: ",Sys.time()))}
    if(verbose){message(paste0("Total antigens (before duplicate removal): ",antigens[,.N]))}
    ##### Clean up From Immunization Processing #####
    rm(immunization_data)
    gc()
    ##### PREPARE ANTIGEN LIST FOR EVALUATION #####
    #Remove all duplication age-antigen combinations
    if(verbose){message("Extracting data for validation. Please wait...")}
    antigens_count<-data.table::as.data.table(unique(antigens,by=c('STUDY_ID','ANTIGEN','AGE_IMM_GIVEN')))
    if(verbose){message(paste0("Unique dose administrations extracetd for each antigen: ",Sys.time()))}
    if(verbose){message(paste0("Unique dose administrations identifed: ",antigens_count[,.N]))}
    # PREPARE list to collect skipped doses
    skipped_list <- list()
    #clean up space
    rm(antigens)
    gc()
    if(verbose){message(paste0("Next step is resource intensive. Please Wait...."))}
    #adjust the counter to evaluate dose #1
    #order the table to support assignment of intervals between doses
    data.table::setkey(antigens_count, STUDY_ID, ANTIGEN, AGE_IMM_GIVEN)
    #create the antigen variables
    antigens_count[,ADMIN_COUNTER:=seq_len(.N),by=.(STUDY_ID, ANTIGEN)]
    antigens_count[,ABS_ADMIN_COUNTER:=ADMIN_COUNTER]
    antigens_count[,TABLE_INDEX:=.I]
    data.table::setindex(antigens_count,TABLE_INDEX)
    #calculate the intervals - optimize based on once for whole table than backing out the cross steps
    antigens_count[, INTERVAL := AGE_IMM_GIVEN - data.table::shift(AGE_IMM_GIVEN)]
    antigens_count[ADMIN_COUNTER==1, INTERVAL := NA_integer_]# reset first rows in each group
    if(verbose){message(paste0("There are ",antigens_count[,.N]," remaining @",Sys.time(),". Moving onto preliminary immunization removals."))}
    
    #create an data set for removed from processing antigens
    skipped_group<-antigens_count[AGE_IMM_GIVEN<0]
    if(nrow(skipped_group)> 0L){
      skipped_group[,NOTES:="Before birth."]#information relevant to handling
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if(verbose){message("Administrations discarded as before birth doses as of today: ",skipped_group[,.N])}
    }
    #future immunization dates skipped
    skipped_group<-antigens_count[AGE_IMM_GIVEN>reference_date]
    if(nrow(skipped_group)> 0L){
      skipped_group[, NOTES := paste0("After ", reference_date)]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if(verbose){message("Administrations discarded as future doses as of ",reference_date,": ",skipped_group[,.N])}
    }
    #immunizaitons skipped removed from count
    
    ### PRE-emptive removal of nonlive given before the earliest allowed date. Minimum intervals are retained and removed after live vaccine validate_immunizations
    skipped_group<-antigens_count[(ANTIGEN=="MCV"&AGE_IMM_GIVEN<yr_with_grace(10))|
                     (ANTIGEN=="MENB"&AGE_IMM_GIVEN<yr_with_grace(16))|
                     (ANTIGEN=="HEPA"&AGE_IMM_GIVEN<yr_with_grace(1))|
                     (ANTIGEN=="HPV"&AGE_IMM_GIVEN<yr_with_grace(9))|
                     (ANTIGEN=="INFLUENZA"&AGE_IMM_GIVEN<mon_with_grace(6))|
                     (ANTIGEN=="COVID"&AGE_IMM_GIVEN<mon_with_grace(6))|
                     (ANTIGEN=="TETANUS"&AGE_IMM_GIVEN<wk_with_grace(6))|
                     (ANTIGEN=="ROTA"&AGE_IMM_GIVEN<wk_with_grace(4))|
                     (ANTIGEN=="HIB"&AGE_IMM_GIVEN<wk_with_grace(6))|
                     (ANTIGEN=="POLIO"&AGE_IMM_GIVEN<wk_with_grace(6))|
                     (ANTIGEN=="PCV"&AGE_IMM_GIVEN<wk_with_grace(6))]
    if(nrow(skipped_group)> 0L){
      skipped_group[,NOTES:="Given before minimum age on US Schedule and not non-enteral live vaccine."]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if(verbose){message(paste0("Administrations discarded as given before first dose allowed: ",skipped_group[,.N]))}
    }
    #ROTA after 8 months (immunize.org as reference)
    skipped_group<-antigens_count[(ANTIGEN=="ROTA"&AGE_IMM_GIVEN>mon_no_grace(8))]
    if(nrow(skipped_group)> 0L){
      skipped_group[,NOTES:="Given after max age on US schedule."]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if(verbose){message(paste0("Administrations discarded as after max age on US schedule: ",skipped_group[,.N]))}
    }
    
    
    #RSV removed from counter if given before release of beyfortus. We do not validate on Palivizumab
    skipped_group<-antigens_count[(ANTIGEN =='RSV'&DATE_GIVEN<RSV_DATE)]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=("Palivizumab doses removed based on RSV date for Beyfortus starting."))]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if(verbose){message(paste0("Administrations discarded as historic Palivizumab doses - not Beyfortus: ",skipped_group[,.N]))}
    }
    
    #process step complete message
    if(verbose){message(paste0("There are ",antigens_count[,.N]," remaining @",Sys.time(),". Moving onto Live Virus processing. Please wait..."))}
    
    ##### --> LIVE VIRUS EVALUATION #####
    # remove Live Vaccines given too close together, nothing is validated as true here - just interval based removals
    # Keep vaccines given too soon as they could potentially interact
    # ACIP guidance is ANY Live vaccine (not enteral) immunization given within 28 days of another Live Vaccine (not enteral) does not generate a sufficient immune response and needs to be repeated
    # the first vaccine given does not need to be repeated
    # return table to ordered table after completion
    if (antigens_count[IS_LIVE == TRUE, .N] == 0) {
      if (verbose) message("No live vaccines detected in dataset. Skipping live vaccine interval checks.")
    } 
    else {
      # Proceed only if there are live doses
      data.table::setkey(antigens_count, STUDY_ID, IS_LIVE, AGE_IMM_GIVEN)
      antigen_live_map <- antigens_count[IS_LIVE == TRUE,.N,by = .(STUDY_ID, AGE_IMM_GIVEN)]
      #set interval
      antigen_live_map[, LIVE_INTERVAL := AGE_IMM_GIVEN - data.table::shift(AGE_IMM_GIVEN, 1L, type = "lag"),by = STUDY_ID]
      # replace NA with 0 for first row per STUDY_ID
      antigen_live_map[is.na(LIVE_INTERVAL), LIVE_INTERVAL := 0]
      antigen_live_map <- antigen_live_map[, .(LIVE_INTERVAL = min(LIVE_INTERVAL, na.rm = TRUE)),by = .(STUDY_ID, AGE_IMM_GIVEN)]
      antigens_count <- antigen_live_map[antigens_count, on = .(STUDY_ID, AGE_IMM_GIVEN)]
      #handle removed
      skipped_group<-antigens_count[IS_LIVE==TRUE&LIVE_INTERVAL>0&LIVE_INTERVAL<MIN_INTERVAL_LIVE]
      if(nrow(skipped_group)> 0L){
        skipped_group[,`:=`(NOTES=("Live vaccine given too soon after other live vaccine."))]
        skipped_list[[length(skipped_list) + 1L]]<-skipped_group
        antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
        if (verbose) {
          live_skipped<-skipped_group[,.N]
          message("Live removal: Total of ", live_skipped," removed. ", nrow(antigens_count), " remaining (total) @", Sys.time())
          message("Live Virus administrations discarded: ", live_skipped)
        }
        
      }
      if (verbose) {message(nrow(antigens_count), " administrations remain to be validated at antigen level starting at @", Sys.time())}
      rm(antigen_live_map)
      data.table::setkey(antigens_count, STUDY_ID, ANTIGEN, AGE_IMM_GIVEN)
    }
    #APRIORI REMOVE imm with INTERVAL < shortest minterval for Antigen
    # HEP A dose inadvertently given less than 6 months after the skipped dose, it does not need to be repeated again as long as the interval between the initial HepA vaccine and the most recent dose is at least 6 calendar months.
    # HEP B TIME BETWEEN DOSE 2 (VALID) and DOSE 3 (VALID) is 8 weeks. An skipped interval admin does not reset this clock as per ACIP.
    skipped_group<-antigens_count[(!ANTIGEN %in% c('COVID','HEPA','HEPB')&INTERVAL<MIN_INTERVAL_DEFAULT)#HEP A and HEP B ignore extra doses for minimum intervals
                                            |(ANTIGEN =='COVID'&INTERVAL<MIN_INTERVAL_COVID_PFIZER)#COVID has a potential minimum interval of 21 days
                                            ]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=("Interval between doses to soon for second dose to be valid."))]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if(verbose){message(paste0("Administrations removed -- given too soon, no additional checking needed: ",skipped_group[,.N]))}
    }
    
    ##### COMPLETED ANTIGEN PRE-PREOCESSING ####
    if(verbose){message(paste0("Beggining the evaluation of antigens by dose counter with first dose in each series: ", Sys.time()))}
    if(verbose){message(paste0("Applying final admin counter. This may take some time. Please wait..."))}
    #### FINAL REBASE BEFORE DOSE CYCLE CALCULATIONS #####
    antigens_count[,`:=`(ADMIN_COUNTER=seq_len(.N)),by=.(STUDY_ID,ANTIGEN)]
    if(verbose){message(paste0("Starting Dose 1"))}
    #PREPARE RETURN storage
    antigens_valid<-vector("list",7)
    ##### DOSE 1 checking #####
    current_baseline<-last_baseline<-nrow(antigens_count)
    counter<-1
    cycle<-0
    while(cycle==0|(current_baseline<last_baseline)){
      cycle<-cycle+1
      last_baseline<-nrow(antigens_count)
      # POLIO: No difference between OPV and IPV schedule currently. tOPV vs bOPV/mOPV should be addressed in antigen assignment
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="POLIO"&AGE_IMM_GIVEN>=wk_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8)),DELAYED=(AGE_IMM_GIVEN>mon_no_grace(4)))]
      #HIB: 
      #First Dose is FINAL DOSE HIB >15 months
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>=mon_no_grace(15),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(4)))]
      #otherwise additional doses are needed
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>=wk_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE
                                                                                                                   ,NEXT_DOSE_MIN=ifelse(AGE_IMM_GIVEN<yr_no_grace(1)
                                                                                                                                         ,(DATE_GIVEN+wk_no_grace(4))
                                                                                                                                         ,(DATE_GIVEN+wk_no_grace(8)))
                                                                                                                   ,NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8))
                                                                                                                   ,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(4)))]
      #HEP B
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPB"&AGE_IMM_GIVEN>=yr_no_grace(0),`:=`(DOSE_COUNTER=counter,VALID=TRUE,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8)),DELAYED=(AGE_IMM_GIVEN>mon_no_grace(4)))]
      #RSV - count any single dose of RSV monovalent Antibody for now
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="RSV"&AGE_IMM_GIVEN>=yr_no_grace(0),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(4)))]
      #TETANUS- ignore Td vs. DT
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=wk_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8)),DELAYED=(AGE_IMM_GIVEN>mon_no_grace(4)))]
      #ROTA - while maximum age to start is 14 weeks 6 days. Still count if late as series should be continued
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="ROTA"&AGE_IMM_GIVEN>=wk_with_grace(4)&AGE_IMM_GIVEN<mon_no_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE
                                                                                                                                                  ,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8)),DELAYED=(AGE_IMM_GIVEN>=wk_no_grace(15)))]
      #PCV: FINAL DOSE PCV after 2 YRS (not considering HR groups)
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=yr_no_grace(2),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(4)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=wk_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX
                                                                                                                   ,NEXT_DOSE_MIN=ifelse(AGE_IMM_GIVEN<yr_no_grace(1)
                                                                                                                                         ,(DATE_GIVEN+wk_no_grace(4))
                                                                                                                                         ,(DATE_GIVEN+wk_no_grace(8)))
                                                                                                                   ,NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8))
                                                                                                                   ,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(4)))]
      #MMR
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MMR"&AGE_IMM_GIVEN>=yr_with_grace(1),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(18))
                                                                                                                   ,NEXT_DOSE_RECOMMENDED=ifelse(AGE_IMM_GIVEN<yr_no_grace(4)
                                                                                                                                        ,(pmax(DOB+yr_no_grace(4),DATE_GIVEN+wk_no_grace(4)))
                                                                                                                                        ,(DATE_GIVEN+wk_no_grace(4)))
                                                                                                                   ,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4))
                                                                                                                   )]
      #VZV
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="VZV"&AGE_IMM_GIVEN>=yr_with_grace(1),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(18))
                                                                                                                   ,NEXT_DOSE_RECOMMENDED=ifelse(AGE_IMM_GIVEN<yr_no_grace(13)
                                                                                                                                          ,pmax(DOB+yr_no_grace(4),DATE_GIVEN+wk_no_grace(12))
                                                                                                                                          ,(DATE_GIVEN+wk_no_grace(4)))
                                                                                                                   ,NEXT_DOSE_MIN=ifelse(AGE_IMM_GIVEN<yr_no_grace(13)
                                                                                                                                         ,DATE_GIVEN+wk_no_grace(12)
                                                                                                                                         ,DATE_GIVEN+wk_no_grace(4))
                                                                                                                   )]
      #HEP A
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPA"&AGE_IMM_GIVEN>=yr_with_grace(1),`:=`(DOSE_COUNTER=counter,VALID=TRUE,NEXT_DOSE_MIN=(DATE_GIVEN+mon_no_grace(6)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+mon_no_grace(6)),DELAYED=(AGE_IMM_GIVEN>mon_no_grace(18)))]
      #HPV
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HPV"&AGE_IMM_GIVEN>=yr_with_grace(9),`:=`(DOSE_COUNTER=counter,VALID=TRUE,NEXT_DOSE_MIN=(DATE_GIVEN+mon_no_grace(6)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+mon_no_grace(6)),DELAYED=(AGE_IMM_GIVEN>yr_no_grace(13)))]
      #MCV: 
      # Adolescents who receive a first dose after their 16th birthday do not need a booster dose unless they become at increased risk for meningococcal disease. 
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MCV"&AGE_IMM_GIVEN>=yr_no_grace(16),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=TRUE)]
      # Guidance updated (2020): ACIP recommends a single dose of MenACWY at age 11 or 12 years followed by a booster dose administered at age 16 years (Table 2). Children who received MenACWY at age 10 years do not need an additional dose at age 11–12 years but should receive the booster dose at age 16 years. Children who received MenACWY before age 10 years and with no ongoing risk for meningococcal disease for which boosters are recommended should still receive MenACWY according to the recommended adolescent schedule, with the first dose at age 11–12 years and a booster dose at age 16 years. For example, a healthy child who received MenACWY at age 9 years because of short-term travel to a country where meningococcal disease is hyperendemic or epidemic and who is not otherwise at increased risk should receive the MenACWY at age 11–12 years according to the recommended ACIP adolescent vaccination schedule. Children who received MenACWY before age 10 years and for whom boosters are recommended because of an ongoing increased risk for meningococcal disease (e.g., those with complement deficiency, HIV infection, or asplenia) should follow the booster schedule for persons at increased risk.
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MCV"&AGE_IMM_GIVEN>=yr_no_grace(10),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(13))
                                                                                                                  ,NEXT_DOSE_RECOMMENDED=ifelse(AGE_IMM_GIVEN<yr_no_grace(16)
                                                                                                                                                ,pmax(DOB+yr_no_grace(16),DATE_GIVEN+wk_no_grace(8))
                                                                                                                                                ,(DATE_GIVEN+wk_no_grace(8)))
                                                                                                                  ,NEXT_DOSE_MIN=ifelse(AGE_IMM_GIVEN<yr_no_grace(16)
                                                                                                                                        ,pmax(DOB+yr_no_grace(16),DATE_GIVEN+wk_no_grace(8))
                                                                                                                                        ,DATE_GIVEN+wk_no_grace(8))
                                                                                                                  )]
      #MENB - SCDM vaccine cannot be delayed
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MENB"&AGE_IMM_GIVEN>=yr_with_grace(16),`:=`(DOSE_COUNTER=counter,VALID=TRUE,NEXT_DOSE_MIN=(DATE_GIVEN+mon_no_grace(6)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+mon_no_grace(6)),DELAYED=FALSE)]
      #INFLU - Yearly vaccine cannot be delayed as first dose
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="INFLUENZA"&AGE_IMM_GIVEN>=yr_no_grace(9),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(2))
                                                                                                                       ,NEXT_DOSE_RECOMMENDED=next_sept1(DATE_GIVEN))]
                                                                                                                                             
                                                                                                                                             
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="INFLUENZA"&AGE_IMM_GIVEN>=mon_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(2)),NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(4)))]
      #COIVD - Yearly vaccine cannot be delayed
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID"&AGE_IMM_GIVEN>=mon_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=FALSE)]
      #antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID2023"&AGE_IMM_GIVEN>=mon_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=FALSE)]
      #store skipped
      skipped_group<-antigens_count[ADMIN_COUNTER==counter&is.na(VALID)]
      if(nrow(skipped_group)> 0L){
        skipped_group[,`:=`(NOTES=(""),COUNTER=counter,CYCLE=cycle)]
        skipped_list[[length(skipped_list) + 1L]]<-skipped_group
        antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
        if (verbose) message("Invalid administrations identified: ", nrow(skipped_group))
      }
      #reset the ADMIN_COUNTER for first run where the first given was skipped
      antigens_count[STUDY_ID %in% skipped_group$STUDY_ID,`:=`(ADMIN_COUNTER=seq_len(.N)+counter-1),by=.(STUDY_ID,ANTIGEN)]
      # find only the combinations that need reset
      current_baseline<-nrow(antigens_count)
      if(verbose)message(paste0("Dose: ",counter," Cycle: ",cycle,". Validated total of ",antigens_count[VALID&DOSE_COUNTER==counter,.N]," antigens for this dose. ",antigens_count[is.na(VALID),.N]," antigen adminisrations still to be tested (total) @",Sys.time()))
    }
    #clean up for next dose
    antigens_count[,SERIES_COMPLETE:=max(DOSE_COMPLETES_SERIES,na.rm = FALSE),by=.(STUDY_ID, ANTIGEN)]#apply series completion to whole set
    #store skipped - series completed previously, no further dose evaluations needed
    skipped_group<-antigens_count[SERIES_COMPLETE==TRUE&DOSE_COUNTER==0&is.na(VALID)]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=paste0("Series completed previously with ",counter," valid doses."),COUNTER=counter)]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if (verbose) message("Not evaluated since series already completed: ", nrow(skipped_group))
    }
    # find only the combinations that need reset
    current_baseline<-nrow(antigens_count)
    antigens_count<-antigens_count[VALID==TRUE|SERIES_COMPLETE==FALSE]#only keep future doses or ADMIN_1s that are correct
    #Post Processing
    antigens_count[VALID==TRUE,FIRST_DOSE_CVX:=CVX]
    antigens_count[,AGE_FIRST_DOSE:=min(AGE_IMM_GIVEN,na.rm = FALSE),by=.(STUDY_ID, ANTIGEN)]#apply series completion to whole set
    antigens_count[,FIRST_DOSE_CVX:=max(FIRST_DOSE_CVX,na.rm = FALSE),by=.(STUDY_ID, ANTIGEN)]#apply first dose CVX
    #antigens_count[ANTIGEN=="PCV",PRIOR_PNEUMOCOCCAL:=max_na(LAST_PNEUMOCOCCAL_CVX,na.rm = TRUE),by=.(STUDY_ID)]#handling of PPV23 and PCV
    antigens_count[ANTIGEN=="HIB",PRIOR_HIB4:=max_na(FIRST_DOSE_CVX %in% cvx$HIB4,na.rm = TRUE),by=.(STUDY_ID)]
    antigens_count[ANTIGEN=="ROTA",PRIOR_ROTA2:=max_na(FIRST_DOSE_CVX %in% cvx$ROTA2,na.rm = TRUE),by=.(STUDY_ID)]
    antigens_count[ANTIGEN=="MENB",PRIOR_BEXSERO:=max_na(FIRST_DOSE_CVX %in% cvx$BEXSERO,na.rm = TRUE),by=.(STUDY_ID)]
    antigens_count[ANTIGEN=="MENB",PRIOR_TRUMENBA:=max_na(FIRST_DOSE_CVX %in% cvx$TRUMENBA,na.rm = TRUE),by=.(STUDY_ID)]
    antigens_count[ANTIGEN=="HEPB"&AGE_IMM_GIVEN>=yr_with_grace(11),PREVIOUS_HEPB_BRAND_ADOL := FIRST_DOSE_CVX %in% cvx$ADOL_HEPB2]
    antigens_count[ANTIGEN=="HEPB"&AGE_IMM_GIVEN>=yr_with_grace(18),PREVIOUS_HEPB_BRAND_ADULT2DOSE := FIRST_DOSE_CVX %in% cvx$ADULT_HEPB]
    #report out
    if(verbose){message(paste0("Administrations validated dose ",counter,": ",antigens_count[VALID==TRUE&DOSE_COUNTER==counter,.N]))}
    
    if(verbose){message(paste0("Total administrations remaining: ",antigens_count[is.na(VALID),.N]))}
    antigens_valid[[counter]]<-antigens_count[VALID==TRUE]
    antigens_count<-antigens_count[is.na(VALID)]
    if(verbose){message(paste0("Completed Dose ",counter," evaluation @ ",Sys.time()))}
    gc()#clean up
    ##### DOSE 2 checking #####
    #adjust the loop variables
    current_baseline<-last_baseline<-nrow(antigens_count)
    counter<-counter+1
    cycle<-0
    while(cycle==0|(current_baseline<last_baseline)){
      cycle<-cycle+1
      last_baseline<-nrow(antigens_count)
      #VALIDATE DOSE 2 
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="POLIO"&AGE_IMM_GIVEN>=wk_with_grace(10)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8)),VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      #HIB completion checks
      #1 ANY 2nd dose after 15 mo of age is completion (8wks since previous to ensure immune response)
      #2 2 doses after First year of life, 8 weeks between
      #3 anything else will require a third dose or 4th dose to complete series
      #ANY DOSE after 15 mo of age is completion (8wks since previous to ensure immune response)
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>mon_no_grace(15)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]  
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>=yr_no_grace(1)&AGE_FIRST_DOSE>=yr_no_grace(1)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]  
      #4 week If current age is younger than 12 months and first dose was administered at younger than age 7 months and at least 1 previous dose was PRP-T (ActHib, Pentacel, Hiberix), Vaxelis or unknown
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>=wk_with_grace(10)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE
                                                                                                                                               ,NEXT_DOSE_MIN=ifelse(AGE_IMM_GIVEN<yr_no_grace(1)&AGE_FIRST_DOSE<mon_no_grace(7)&(FIRST_DOSE_CVX %in% cvx$HIB4 | CVX %in% cvx$HIB4)
                                                                                                                                                                 ,(DATE_GIVEN+wk_no_grace(4))
                                                                                                                                                                 ,(DATE_GIVEN+wk_no_grace(8)))#must already be above age 12mo since starting age is given
                                                                                                                                               ,NEXT_DOSE_RECOMMENDED=ifelse(AGE_IMM_GIVEN<yr_no_grace(1)&AGE_FIRST_DOSE<mon_no_grace(7)&(FIRST_DOSE_CVX %in% cvx$HIB4|CVX %in% cvx$HIB4)
                                                                                                                                                                             ,ifelse(AGE_IMM_GIVEN<mon_no_grace(7),(DATE_GIVEN+wk_no_grace(8)),(DATE_GIVEN+wk_no_grace(4)))
                                                                                                                                                                             ,(DATE_GIVEN+wk_no_grace(8)))#must already be above age 12mo since starting age is given
                                                                                                                                               ,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      #hepb dose 2 can be 4 weeks after dose 1, no age minimum yet
      #HEPB 2 dose adol series (4 months)
      #hep B adult 2 dose series (1 month)
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPB"&AGE_IMM_GIVEN>=yr_with_grace(11)&AGE_IMM_GIVEN<=yr_with_grace(16)&INTERVAL>=mon_with_grace(4)&PREVIOUS_HEPB_BRAND_ADOL&CVX %in% cvx$ADOL_HEPB2&CVX %in% cvx$ADOL_HEPB2 &INTERVAL>=mon_with_grace(4),`:=`(HEPB_DOSE_2=AGE_IMM_GIVEN,DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPB"&AGE_IMM_GIVEN>=yr_with_grace(18)&INTERVAL>=wk_with_grace(4)&FIRST_DOSE_CVX %in% cvx$ADULT_HEPB2 &CVX %in% cvx$ADULT_HEPB2,`:=`(HEPB_DOSE_2=AGE_IMM_GIVEN,DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      #HEP B dose 3 next dose is min x weeks from dose 2 and 16 weeks from dose 1,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8))
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPB"&AGE_IMM_GIVEN>=yr_no_grace(0)&INTERVAL>=wk_with_grace(4),`:=`(HEPB_DOSE_2=AGE_IMM_GIVEN,DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4))
                                                                                                                                             ,NEXT_DOSE_MIN=pmax(DOB+wk_no_grace(24),DATE_GIVEN+wk_no_grace(8),DOB+AGE_FIRST_DOSE+wk_no_grace(16))#
                                                                                                                                             ,NEXT_DOSE_RECOMMENDED=pmax(DOB+wk_no_grace(24),DATE_GIVEN+wk_no_grace(8),DOB+AGE_FIRST_DOSE+wk_no_grace(16))#
                                                                                                                                             )]
      #TETANUS dose 1-2 is 4 week minimum, age is not a factor yet
      #Age is a factor in dose 2-3 if over age 7 for dose 2
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=wk_with_grace(10)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4))
                                                                                                                                                   ,NEXT_DOSE_MIN=pmax(DATE_GIVEN+wk_no_grace(4))#
                                                                                                                                                   ,NEXT_DOSE_RECOMMENDED=pmax(DOB+wk_no_grace(24),DATE_GIVEN+wk_no_grace(4))#
                                                                                                                                                   )]
      #ROTA series complete if both 2 dose series
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="ROTA"&AGE_IMM_GIVEN>=wk_with_grace(10)&AGE_IMM_GIVEN<mon_no_grace(8)&INTERVAL>=wk_with_grace(4)&IS_ROTA2==TRUE&PRIOR_ROTA2==TRUE,`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="ROTA"&AGE_IMM_GIVEN>=wk_with_grace(10)&AGE_IMM_GIVEN<mon_no_grace(8)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8)),DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      #PCV as final if dose 2 >24 mo OR first dose after first yr_no_grace(0)
      #PCV 3rd dose needed (at least) ,NEXT_DOSE_MIN=(DATE_GIVEN+wk_no_grace(4)),NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+wk_no_grace(8)) but if age 1yr weeks would have min and recommended at 8 weeks
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&INTERVAL>=wk_with_grace(8)&AGE_IMM_GIVEN>=yr_no_grace(2),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=yr_no_grace(1)&AGE_FIRST_DOSE>=yr_no_grace(1)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=wk_with_grace(10)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE
                                                                                                                                               ,NEXT_DOSE_MIN=ifelse(AGE_IMM_GIVEN<yr_no_grace(1)&AGE_FIRST_DOSE<mon_no_grace(7)
                                                                                                                                                                     ,(DATE_GIVEN+wk_no_grace(4))
                                                                                                                                                                     ,pmax(DOB+yr_no_grace(1),DATE_GIVEN+wk_no_grace(8)))
                                                                                                                                               ,NEXT_DOSE_RECOMMENDED=ifelse(AGE_IMM_GIVEN<yr_no_grace(1)&AGE_FIRST_DOSE<mon_no_grace(7)
                                                                                                                                                                             ,(DATE_GIVEN+wk_no_grace(8))
                                                                                                                                                                             ,pmax(DOB+yr_no_grace(1),DATE_GIVEN+wk_no_grace(8)))
                                                                                                                                               ,LAST_PNEUMOCOCCAL_CVX=CVX,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      #as per CDC/ACIP 2022-1-25 checked, 4 wks minimum dose 2, CAN BE GIVEN BEFORE AGE 4
      #TODO: MMRV has minimum timing of 3 months - THIS IS NOT CHECKED FOR YET
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MMR"&AGE_IMM_GIVEN>=yr_with_grace(1)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="VZV"&AGE_IMM_GIVEN>=yr_with_grace(13)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(6)&INTERVAL>mon_with_grace(4)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="VZV"&AGE_IMM_GIVEN>=yr_with_grace(1)&INTERVAL>=wk_with_grace(12),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(6)&INTERVAL>mon_with_grace(15)))]
      #HEPA A adult has 3 dose series 
      #TODO: ADULT HEPA ANTIGEN, 2 doses completes series (currently not implemented)
      #HEP A DOSE 2 (VALID) is 6 months from Dose 1 (VALID) - it ignores skipped interval doses
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPA"&AGE_IMM_GIVEN>=yr_with_grace(1)&AGE_IMM_GIVEN<yr_no_grace(19)&(AGE_IMM_GIVEN-AGE_FIRST_DOSE)>=mon_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(2)&INTERVAL>mon_with_grace(15)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPA"&IS_ADULT_HEPA==TRUE&AGE_IMM_GIVEN>yr_no_grace(19)&(AGE_IMM_GIVEN-AGE_FIRST_DOSE)>=mon_with_grace(6),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(2)&INTERVAL>mon_with_grace(15)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPA"&AGE_IMM_GIVEN>=yr_no_grace(19)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(2)&INTERVAL>mon_with_grace(15))
                                                                                                                                              ,NEXT_DOSE_RECOMMENDED=pmax(DATE_GIVEN+mon_no_grace(5),DOB+AGE_FIRST_DOSE+mon_no_grace(6))
                                                                                                                                              ,NEXT_DOSE_MIN=pmax(DATE_GIVEN+mon_no_grace(5),DOB+AGE_FIRST_DOSE+mon_no_grace(6))
                                                                                                                                              )]
      #HPV: Completion of series if first dose given before 15th birthday
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HPV"&AGE_IMM_GIVEN>=yr_with_grace(9)&AGE_FIRST_DOSE<yr_no_grace(15)&INTERVAL>=mon_with_grace(5),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(9)&INTERVAL>mon_with_grace(15)))]
      #HPV: otherwise count as second dose if 4 weeks/grave from dose 1
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HPV"&AGE_IMM_GIVEN>=yr_with_grace(9)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(9)&INTERVAL>mon_with_grace(15))
                                                                                                                                              ,NEXT_DOSE_RECOMMENDED=pmax(DATE_GIVEN+wk_no_grace(16),DOB+AGE_FIRST_DOSE+mon_no_grace(6))
                                                                                                                                              ,NEXT_DOSE_MIN=pmax(DATE_GIVEN+wk_no_grace(12),DOB+AGE_FIRST_DOSE+wk_no_grace(21))
                                                                                                                                              )]
      #MCV: 
      # Adolescents who receive their first dose at age 13–15 years should receive a booster dose at age 16–18 years; 
      # the booster dose can be administered at any time, as long as a minimum interval of 8 weeks between doses is maintained. 
      # Adolescents who receive a first dose after their 16th birthday do not need a booster dose unless they become at increased risk for meningococcal disease. 
      # Persons aged 19–21 years who have not received a dose after their 16th yr_no_grace(0)day can receive a single MenACWY dose as part of catch-up vaccination. MenACWY vaccines are interchangeable; the same vaccine product is recommended, but not required, for all doses. MenACWY vaccines can be administered simultaneously with other vaccines indicated for this age group, but at a different anatomic site, if feasible. MenACWY-TT, which is conjugated to tetanus toxoid, is only licensed for the prevention of meningococcal disease; use of this vaccine does not replace doses or affect the dosing intervals of routinely recommended tetanus toxoid–containing vaccines in any age group.
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MCV"&AGE_IMM_GIVEN>=yr_with_grace(16)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(17)&INTERVAL>mon_with_grace(6)))]
      #MENB:Note: MenB-FHbp and MenB-4C are not interchangeable
      #MENB doses given > 4 wks after dose 1 are valid historically
      # Otherwise 6 months recommended. If given 4 months after previous dose and multiple doses that is OK. no additional time from first dose. Only applies to first dose >16
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MENB"&AGE_IMM_GIVEN>=yr_with_grace(16)&INTERVAL>=wk_with_grace(4)&PRIOR_BEXSERO&IS_BEXSERO&DATE_GIVEN<MENB_DATE,`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(16)&INTERVAL>mon_with_grace(15)))]      
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MENB"&AGE_IMM_GIVEN>=yr_with_grace(16)&INTERVAL>=mon_with_grace(6)&PRIOR_BEXSERO&IS_BEXSERO,`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(16)&INTERVAL>mon_with_grace(15)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MENB"&AGE_IMM_GIVEN>=yr_with_grace(16)&INTERVAL>=mon_with_grace(4)&PRIOR_BEXSERO&IS_BEXSERO&ABS_ADMIN_COUNTER>=3&AGE_FIRST_DOSE>=yr_with_grace(16),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(16)&INTERVAL>mon_with_grace(15)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="MENB"&AGE_IMM_GIVEN>=yr_with_grace(16)&INTERVAL>=mon_with_grace(6)&PRIOR_TRUMENBA&IS_TRUMENBA,`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(16)&INTERVAL>mon_with_grace(15)))]
      #INFLUENZA BOOSTER AND YEARLY CHECK
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="INFLUENZA",`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(2)),NEXT_DOSE_RECOMMENDED=next_sept1(DATE_GIVEN))]
      #need data on first_COVID_brand
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID"&FIRST_DOSE_CVX %in% cvx$PFIZER&AGE_IMM_GIVEN>=mon_with_grace(6)&INTERVAL>=wk_with_grace(3),`:=`(DOSE_COUNTER=counter,VALID=TRUE)]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID"&!FIRST_DOSE_CVX %in% cvx$PFIZER&AGE_IMM_GIVEN>=mon_with_grace(6)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE)]
      #store skipped
      skipped_group<-antigens_count[ADMIN_COUNTER==counter&is.na(VALID)]
      if(nrow(skipped_group)> 0L){
        skipped_group[,`:=`(NOTES=(""),COUNTER=counter,CYCLE=cycle)]
        skipped_list[[length(skipped_list) + 1L]]<-skipped_group
        antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
        if (verbose) message("Invalid administrations identified: ", nrow(skipped_group))
      }
      #reset the ADMIN_COUNTER for first run where the first given was skipped
      antigens_count[STUDY_ID %in% skipped_group$STUDY_ID,`:=`(ADMIN_COUNTER=seq_len(.N)+counter-1),by=.(STUDY_ID,ANTIGEN)]
      current_baseline<-nrow(antigens_count)
      if(verbose)message(paste0("Dose: ",counter," Cycle: ",cycle,". Validated total of ",antigens_count[VALID&DOSE_COUNTER==counter,.N]," antigens for this dose. ",antigens_count[is.na(VALID),.N]," antigen adminisrations still to be tested (total) @",Sys.time()))
    }
    #clean up and prep for next dose
    antigens_count[,SERIES_COMPLETE:=max(DOSE_COMPLETES_SERIES,na.rm = FALSE),by=.(STUDY_ID, ANTIGEN)]#apply series completion to whole set
    #store skipped - series completed previously, no further dose evaluations needed
    skipped_group<-antigens_count[SERIES_COMPLETE==TRUE&DOSE_COUNTER==0&is.na(VALID)]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=paste0("Series completed previously with ",counter," valid doses."),COUNTER=counter)]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if (verbose) message("Not evaluated since series already completed: ", nrow(skipped_group))
    }
    #Dose post-processing - HEP B, HIB and PCV care about second dose status
    #antigens_count[ANTIGEN=="PCV"&ADMIN_COUNTER>=counter,PRIOR_PNEUMOCOCCAL:=max_na(LAST_PNEUMOCOCCAL_CVX,na.rm = TRUE),by=.(STUDY_ID)]#handling of PPV23 and PCV
    antigens_count[ANTIGEN=="HEPB"&ADMIN_COUNTER>=counter,SECOND_HEPB:=max_na(HEPB_DOSE_2,na.rm = TRUE),by=.(STUDY_ID)]
    antigens_count[ADMIN_COUNTER<=(counter+1)&ANTIGEN=="HIB",EXISTS_HIB4:=any(ifelse(is.na(IS_HIB4), F, IS_HIB4)),by=.(STUDY_ID)]#force HIB 4 calculations if dose 1/2/3 are HIB4
    #report out
    if(verbose){message(paste0("Administrations validated dose ",counter,": ",antigens_count[VALID==TRUE&DOSE_COUNTER==counter,.N]))}
    
    if(verbose){message(paste0("Total administrations remaining: ",antigens_count[is.na(VALID),.N]))}
    antigens_valid[[counter]]<-antigens_count[VALID==TRUE]
    antigens_count<-antigens_count[is.na(VALID)]
    if(verbose){message(paste0("Completed Dose ",counter," evaluation @ ",Sys.time()))}
    gc()#clean up
    ##### DOSE 3 checking #####
    current_baseline<-last_baseline<-nrow(antigens_count)
    counter<-counter+1
    cycle<-0
    while(cycle==0|(current_baseline<last_baseline)){
      cycle<-cycle+1#cycle counter
      last_baseline<-nrow(antigens_count)
      #VALIDATE DOSE 3
      #POLIO 3rd dose validate_immunizations pre/post 4th yr_no_grace(0)day
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="POLIO"&AGE_IMM_GIVEN>=yr_no_grace(4)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(12)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="POLIO"&AGE_IMM_GIVEN>=wk_with_grace(14)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6))
                                                                                                                                                 ,NEXT_DOSE_MIN=pmax(DOB+yr_no_grace(4),DATE_GIVEN+mon_no_grace(6))
                                                                                                                                                 ,NEXT_DOSE_RECOMMENDED=pmax(DOB+yr_no_grace(4),DATE_GIVEN+mon_no_grace(6))
                                                                                                                                                 )]
      #third dose after age 1 yr completes series if all first doses are HIB3
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&EXISTS_HIB4==FALSE&AGE_IMM_GIVEN>=yr_with_grace(1)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(2)&INTERVAL>mon_with_grace(6)))]
      #third dose completes series if 8 weeks since previous, is 12+ mo and first dose was administered at age 7 through 11 months;
      #third dose completes series if 8 weeks since previous, is 12+ mo and second dose was administered 12+ mo;
      #otherwise 4th dose is needed  
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>=yr_no_grace(1)&AGE_FIRST_DOSE<yr_no_grace(1)&AGE_FIRST_DOSE>=mon_no_grace(7)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>=yr_no_grace(1)&AGE_FIRST_DOSE<yr_no_grace(1)&(AGE_IMM_GIVEN-INTERVAL)>=yr_no_grace(1)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>=wk_with_grace(14)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE
                                                                                                                                               ,NEXT_DOSE_MIN=pmax(DOB+yr_no_grace(1),(DATE_GIVEN+wk_no_grace(8)))
                                                                                                                                               ,NEXT_DOSE_RECOMMENDED=pmax(DOB+mon_no_grace(15),(DATE_GIVEN+wk_no_grace(8)))
                                                                                                                                               ,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6)))]
      #hepb dose 3 is 8 weeks after a valid dose 2 and 16 weeks after dose 1 and given at >24 weeks
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPB"&AGE_IMM_GIVEN>=wk_no_grace(24)&(AGE_IMM_GIVEN-SECOND_HEPB)>=wk_with_grace(8)&(AGE_IMM_GIVEN-AGE_FIRST_DOSE)>=wk_with_grace(16),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6)))]
      
      #PRIMARY TETANUS  COMPLETE WITH DOSE 3 AT 7yrs if 2 doses after first yr_no_grace(0)day - 
      #HOWEVER, DOSE after AGE 10 is required, therefore
      #SERIES IS NEVER COMPLETE as needed q10 years
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(10)&AGE_FIRST_DOSE>=yr_no_grace(1)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(13)&INTERVAL>mon_with_grace(15))
                                                                                                                                                                                    ,NEXT_DOSE_MIN=DATE_GIVEN+yr_no_grace(10)
                                                                                                                                                                                    ,NEXT_DOSE_RECOMMENDED=DATE_GIVEN+yr_no_grace(10)
                                                                                                                                                                                    )]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=wk_with_grace(14)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&AGE_IMM_GIVEN<yr_no_grace(11)&INTERVAL>mon_with_grace(6))
                                                                                                                                                   ,NEXT_DOSE_MIN=pmax(DOB+yr_no_grace(4),DATE_GIVEN+mon_no_grace(6))
                                                                                                                                                   ,NEXT_DOSE_RECOMMENDED=pmax(DOB+yr_no_grace(4),DATE_GIVEN+mon_no_grace(6))
                                                                                                                                                   )]
      #ROTA 
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="ROTA"&AGE_IMM_GIVEN>=wk_with_grace(14)&AGE_IMM_GIVEN<mon_no_grace(8)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(7)&INTERVAL>mon_with_grace(4)))]
      #PCV dose 3 >24 months OR 
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=yr_no_grace(2)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=yr_no_grace(1)&AGE_FIRST_DOSE>=mon_no_grace(7)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=yr_no_grace(1)&(AGE_IMM_GIVEN-INTERVAL)>=yr_no_grace(1)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=wk_with_grace(14)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX
                                                                                                                                               ,NEXT_DOSE_MIN=pmax(DOB+yr_no_grace(1),(DATE_GIVEN+wk_no_grace(8)))
                                                                                                                                               ,NEXT_DOSE_RECOMMENDED=pmax(DOB+mon_no_grace(15),(DATE_GIVEN+wk_no_grace(8)))
                                                                                                                                               ,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(6)))]
      #as per CDC/ACIP 2022-1-25 checked, 4 wks minimum dose 2-3, CAN BE GIVEN BEFORE AGE 4
      #need_data_first HPV age for COMPLETION ONLY, not for validate_immunizations of first dose
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HPV"&AGE_IMM_GIVEN>=yr_with_grace(9)&(AGE_IMM_GIVEN-AGE_FIRST_DOSE)>=wk_with_grace(21)&INTERVAL>=wk_with_grace(12),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(13)&INTERVAL>mon_with_grace(15)))]
      #HEPA A dose 3 is 6 months from dose #1 and 5 months from dose #2
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HEPA"&AGE_IMM_GIVEN>=yr_no_grace(19)&(AGE_IMM_GIVEN-AGE_FIRST_DOSE)>=mon_with_grace(6)&INTERVAL>=mon_with_grace(5),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(15)))]
      #Influenza - Booster and Yearly
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="INFLUENZA",`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(2))
                                                                                         ,NEXT_DOSE_RECOMMENDED=next_sept1(DATE_GIVEN))]
                                                                                                                       
                                                                                                                       
      #need data on first_COVID_brand
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID"&FIRST_DOSE_CVX %in% cvx$PFIZER&AGE_IMM_GIVEN>=mon_with_grace(6)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE)]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID"&!FIRST_DOSE_CVX %in% cvx$PFIZER&AGE_IMM_GIVEN>=mon_with_grace(6)&INTERVAL>=wk_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE)]
      #store skipped
      skipped_group<-antigens_count[ADMIN_COUNTER==counter&is.na(VALID)]
      if(nrow(skipped_group)> 0L){
        skipped_group[,`:=`(NOTES=(""),COUNTER=counter,CYCLE=cycle)]
        skipped_list[[length(skipped_list) + 1L]]<-skipped_group
        antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
        if (verbose) message("Invalid administrations identified: ", nrow(skipped_group))
      }
      #reset the ADMIN_COUNTER for first run where the first given was skipped
      antigens_count[STUDY_ID %in% skipped_group$STUDY_ID,`:=`(ADMIN_COUNTER=seq_len(.N)+counter-1),by=.(STUDY_ID,ANTIGEN)]
      # find only the combinations that need reset
      current_baseline<-nrow(antigens_count)
      if(verbose)message(paste0("Dose ",counter," Cycle #",cycle," validated total of ",antigens_count[VALID&DOSE_COUNTER==counter,.N]," antigens ",antigens_count[is.na(VALID),.N]," remaining @",Sys.time()))
    }
    #clean up and prep for next dose
    antigens_count[,SERIES_COMPLETE:=max(DOSE_COMPLETES_SERIES,na.rm = FALSE),by=.(STUDY_ID, ANTIGEN)]#apply series completion to whole set
    #store skipped - series completed previously, no further dose evaluations needed
    skipped_group<-antigens_count[SERIES_COMPLETE==TRUE&DOSE_COUNTER==0&is.na(VALID)]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=paste0("Series completed previously with ",counter," valid doses."),COUNTER=counter)]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if (verbose) message("Not evaluated since series already completed: ", nrow(skipped_group))
    }
    #Dose post-processing
    #antigens_count[ANTIGEN=="PCV"&ADMIN_COUNTER>=counter,PRIOR_PNEUMOCOCCAL:=max_na(LAST_PNEUMOCOCCAL_CVX,na.rm = TRUE),by=.(STUDY_ID)]#handling of PPV23 and PCV
    #report out
    if(verbose){message(paste0("Administrations validated dose ",counter,": ",antigens_count[VALID==TRUE&DOSE_COUNTER==counter,.N]))}
    
    if(verbose){message(paste0("Total administrations remaining: ",antigens_count[is.na(VALID),.N]))}
    antigens_valid[[counter]]<-antigens_count[VALID==TRUE]
    antigens_count<-antigens_count[is.na(VALID)]
    if(verbose){message(paste0("Completed Dose ",counter," evaluation @ ",Sys.time()))}
    gc()#clean up
    ##### DOSE 4 checking ######
    current_baseline<-last_baseline<-nrow(antigens_count)
    counter<-counter+1
    cycle<-0
    while(cycle==0|(current_baseline<last_baseline)){
      cycle<-cycle+1#cycle counter
      last_baseline<-nrow(antigens_count)
      #IPV
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="POLIO"&AGE_IMM_GIVEN>=yr_no_grace(4)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(6)&INTERVAL>mon_with_grace(15)))]
      #HIB
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="HIB"&AGE_IMM_GIVEN>=yr_no_grace(1)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(18)&INTERVAL>mon_with_grace(6)))]
      #PCV
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&AGE_IMM_GIVEN>=yr_no_grace(1)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(18)&INTERVAL>mon_with_grace(6)))]
      #PPV23
      #antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PCV"&PRIOR_PNEUMOCOCCAL %in% cvx$PCV &CVX %in% cvx$PPV23 &AGE_IMM_GIVEN>=yr_no_grace(2)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(3)&INTERVAL>mon_with_grace(6)))]
      #TETANUS 4 doses counts if last is >7yr (operationalized as 10yr due to HEDIS) - COMPLETES SERIES UNTIL ADULTHOOD
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(10)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(13)&INTERVAL>mon_with_grace(15))
                                                                                                                                                                                              ,NEXT_DOSE_RECOMMENDED=(DATE_GIVEN+yr_no_grace(10))
                                                                                                                                                                                              ,NEXT_DOSE_MIN=(DATE_GIVEN+yr_no_grace(10))
                                                                                                                                                     )]
      #CDC/ACIP: Dose 5 (DTAP) is not necessary if dose 4 was administered at age 4 years or older and at least 6 months after dose 3.
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(4)&AGE_IMM_GIVEN<yr_no_grace(7)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&AGE_IMM_GIVEN<yr_no_grace(11)&INTERVAL>mon_with_grace(15))
                                                                                                                                                                                 ,NEXT_DOSE_MIN=pmax(DOB+yr_no_grace(10),DATE_GIVEN+mon_no_grace(6))
                                                                                                                                                                                 ,NEXT_DOSE_RECOMMENDED=pmax(DOB+yr_no_grace(11),DATE_GIVEN+mon_no_grace(6)))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(7)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(1)&INTERVAL>mon_with_grace(15)&AGE_IMM_GIVEN<yr_no_grace(11))
                                                                                                                                                    ,NEXT_DOSE_MIN=DATE_GIVEN+yr_no_grace(10)
                                                                                                                                                    ,NEXT_DOSE_RECOMMENDED=DATE_GIVEN+yr_no_grace(10))]
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(1)&AGE_IMM_GIVEN<yr_no_grace(7)&INTERVAL>=mon_with_grace(4),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>mon_no_grace(18)&AGE_IMM_GIVEN<yr_no_grace(11)&INTERVAL>mon_with_grace(6))
                                                                                                                                                                              ,NEXT_DOSE_MIN=pmax(DOB+yr_no_grace(4),DATE_GIVEN+mon_no_grace(6))
                                                                                                                                                                              ,NEXT_DOSE_RECOMMENDED=pmax(DOB+yr_no_grace(4),DATE_GIVEN+mon_no_grace(6)))]
      #Influenza - Yearly
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="INFLUENZA",`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(2))
                                                                                         ,NEXT_DOSE_RECOMMENDED=next_sept1(DATE_GIVEN))]
                                                                                                                       
                                                                                                                       
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID"&AGE_IMM_GIVEN>=mon_with_grace(6)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE)]
      #store skipped
      skipped_group<-antigens_count[ADMIN_COUNTER==counter&is.na(VALID)]
      if(nrow(skipped_group)> 0L){
        skipped_group[,`:=`(NOTES=(""),COUNTER=counter,CYCLE=cycle)]
        skipped_list[[length(skipped_list) + 1L]]<-skipped_group
        antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
        if (verbose) message("Invalid administrations identified: ", nrow(skipped_group))
      }
      #reset the ADMIN_COUNTER for first run where the first given was skipped
      antigens_count[STUDY_ID %in% skipped_group$STUDY_ID,`:=`(ADMIN_COUNTER=seq_len(.N)+counter-1),by=.(STUDY_ID,ANTIGEN)]
      # find only the combinations that need reset
      current_baseline<-nrow(antigens_count)
      if(verbose)message(paste0("Dose: ",counter," Cycle: ",cycle,". Validated total of ",antigens_count[VALID&DOSE_COUNTER==counter,.N]," antigens for this dose. ",antigens_count[is.na(VALID),.N]," antigen adminisrations still to be tested (total) @",Sys.time()))
    }
    #clean up and prep for next dose
    antigens_count[,SERIES_COMPLETE:=max(DOSE_COMPLETES_SERIES,na.rm = FALSE),by=.(STUDY_ID, ANTIGEN)]#apply series completion to whole set
    #store skipped - series completed previously, no further dose evaluations needed
    skipped_group<-antigens_count[SERIES_COMPLETE==TRUE&DOSE_COUNTER==0&is.na(VALID)]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=paste0("Series completed previously with ",counter," valid doses."),COUNTER=counter)]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if (verbose) message("Not evaluated since series already completed: ", nrow(skipped_group))
    }
    #Dose post-processing
    #antigens_count[ANTIGEN=="PCV"&ADMIN_COUNTER>=counter,PRIOR_PNEUMOCOCCAL:=max_na(LAST_PNEUMOCOCCAL_CVX,na.rm = TRUE),by=.(STUDY_ID)]#handling of PPV23 and PCV
    #report out
    if(verbose){message(paste0("Administrations validated dose ",counter,": ",antigens_count[VALID==TRUE&DOSE_COUNTER==counter,.N]))}
    
    if(verbose){message(paste0("Total administrations remaining: ",antigens_count[is.na(VALID),.N]))}
    antigens_valid[[counter]]<-antigens_count[VALID==TRUE]
    antigens_count<-antigens_count[is.na(VALID)]
    if(verbose){message(paste0("Completed Dose ",counter," evaluation @ ",Sys.time()))}
    gc()#clean up
    ##### DOSE 5 checking ######
    current_baseline<-last_baseline<-nrow(antigens_count)
    counter<-counter+1
    cycle<-0
    while(cycle==0|(current_baseline<last_baseline)){
      cycle<-cycle+1#cycle counter
      last_baseline<-nrow(antigens_count)
      #TETANUS
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(10)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(13)&INTERVAL>mon_with_grace(15))
                                                                                                                                                     ,NEXT_DOSE_MIN=DATE_GIVEN+yr_no_grace(10)
                                                                                                                                                     ,NEXT_DOSE_RECOMMENDED=DATE_GIVEN+yr_no_grace(10))]
      #CDC/ACIP: Dose 5 (DTAP) is not necessary if dose 4 was administered at age 4 years or older and at least 6 months after dose 3.
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(4)&AGE_IMM_GIVEN<yr_no_grace(7)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(6)&AGE_IMM_GIVEN<yr_no_grace(11)&INTERVAL>mon_with_grace(15))
                                                                                                                                                                                 ,NEXT_DOSE_MIN=pmax(DOB+yr_no_grace(10),DATE_GIVEN+mon_no_grace(6))
                                                                                                                                                                                 ,NEXT_DOSE_RECOMMENDED=pmax(DOB+yr_no_grace(11),DATE_GIVEN+mon_no_grace(6)))]
      #PPV23 - DOSE 1
      #antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PPV23"&PRIOR_PNEUMOCOCCAL %in% cvx$PCV &CVX %in% cvx$PPV23 &AGE_IMM_GIVEN>=yr_no_grace(2)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(3)&INTERVAL>mon_with_grace(6)))]
      #PPV23 - DOSE 2
      #antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="PPV23"&PRIOR_PNEUMOCOCCAL %in% cvx$PPV23&CVX %in% cvx$PPV23&AGE_IMM_GIVEN>=yr_no_grace(5)&INTERVAL>=yr_no_grace(5),`:=`(DOSE_COUNTER=counter,VALID=TRUE,LAST_PNEUMOCOCCAL_CVX=CVX,DOSE_COMPLETES_SERIES=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(3)&INTERVAL>mon_with_grace(6)))]
      #Influenza - Yearly
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="INFLUENZA",`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(2))
                                                                                         ,NEXT_DOSE_RECOMMENDED=next_sept1(DATE_GIVEN))]
                                                                                                                       
                                                                                                                       
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID"&AGE_IMM_GIVEN>=mon_with_grace(6)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE)]
      #store skipped
      skipped_group<-antigens_count[ADMIN_COUNTER==counter&is.na(VALID)]
      if(nrow(skipped_group)> 0L){
        skipped_group[,`:=`(NOTES=(""),COUNTER=counter,CYCLE=cycle)]
        skipped_list[[length(skipped_list) + 1L]]<-skipped_group
        antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
        if (verbose) message("Invalid administrations identified: ", nrow(skipped_group))
      }
      #reset the ADMIN_COUNTER for first run where the first given was skipped
      antigens_count[STUDY_ID %in% skipped_group$STUDY_ID,`:=`(ADMIN_COUNTER=seq_len(.N)+counter-1),by=.(STUDY_ID,ANTIGEN)]
      current_baseline<-nrow(antigens_count)
      if(verbose)message(paste0("Dose: ",counter," Cycle: ",cycle,". Validated total of ",antigens_count[VALID&DOSE_COUNTER==counter,.N]," antigens for this dose. ",antigens_count[is.na(VALID),.N]," antigen adminisrations still to be tested (total) @",Sys.time()))
    }
    #clean up and prep for next dose
    antigens_count[,SERIES_COMPLETE:=max(DOSE_COMPLETES_SERIES,na.rm = FALSE),by=.(STUDY_ID, ANTIGEN)]
    #store skipped - series completed previously, no further dose evaluations needed
    skipped_group<-antigens_count[SERIES_COMPLETE==TRUE&DOSE_COUNTER==0&is.na(VALID)]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=paste0("Series completed previously with ",counter," valid doses."),COUNTER=counter)]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if (verbose) message("Not evaluated since series already completed: ", nrow(skipped_group))
    }
    #Dose 5 post-processing
    #antigens_count[ANTIGEN=="PCV"&ADMIN_COUNTER>=counter,PRIOR_PNEUMOCOCCAL:=max_na(LAST_PNEUMOCOCCAL_CVX,na.rm = TRUE),by=.(STUDY_ID)]#handling of PPV23 and PCV
    #report out
    if(verbose){message(paste0("Administrations validated dose ",counter,": ",antigens_count[VALID==TRUE&DOSE_COUNTER==counter,.N]))}
    
    if(verbose){message(paste0("Total administrations remaining: ",antigens_count[is.na(VALID),.N]))}
    antigens_valid[[counter]]<-antigens_count[VALID==TRUE]
    antigens_count<-antigens_count[is.na(VALID)]
    if(verbose){message(paste0("Completed Dose ",counter," evaluation @ ",Sys.time()))}
    gc()#clean up
    ##### DOSE 6 checking ######
    current_baseline<-last_baseline<-nrow(antigens_count)
    counter<-counter+1
    cycle<-0
    while(cycle==0|(current_baseline<last_baseline)){
      cycle<-cycle+1#cycle counter
      last_baseline<-nrow(antigens_count)
      #TETANUS - Never complete
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(10)&INTERVAL>=days_with_grace(180),`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(AGE_IMM_GIVEN>yr_no_grace(13)&INTERVAL>mon_with_grace(15)),NEXT_DOSE_MIN=DATE_GIVEN+yr_no_grace(10)
                                                                                                                                                     ,NEXT_DOSE_RECOMMENDED=DATE_GIVEN+yr_no_grace(10))]
      #Influenza - Yearly
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="INFLUENZA",`:=`(DOSE_COUNTER=counter,VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(2))
                                                                                         ,NEXT_DOSE_RECOMMENDED=next_sept1(DATE_GIVEN))]
                                                                                                                       
                                                                                                                       
      antigens_count[is.na(VALID)==TRUE&ADMIN_COUNTER==counter&ANTIGEN=="COVID"&AGE_IMM_GIVEN>=mon_with_grace(6)&INTERVAL>=wk_with_grace(8),`:=`(DOSE_COUNTER=counter,VALID=TRUE)]
      #store skipped
      skipped_group<-antigens_count[ADMIN_COUNTER==counter&is.na(VALID)]
      if(nrow(skipped_group)> 0L){
        skipped_group[,`:=`(NOTES=(""),COUNTER=counter,CYCLE=cycle)]
        skipped_list[[length(skipped_list) + 1L]]<-skipped_group
        antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
        if (verbose) message("Invalid administrations identified: ", nrow(skipped_group))
      }
      #reset the ADMIN_COUNTER for first run where the first given was skipped
      antigens_count[STUDY_ID %in% skipped_group$STUDY_ID,`:=`(ADMIN_COUNTER=seq_len(.N)+counter-1),by=.(STUDY_ID,ANTIGEN)]
      current_baseline<-nrow(antigens_count)
      if(verbose)message(paste0("Dose: ",counter," Cycle: ",cycle,". Validated total of ",antigens_count[VALID&DOSE_COUNTER==counter,.N]," antigens for this dose. ",antigens_count[is.na(VALID),.N]," antigen adminisrations still to be tested (total) @",Sys.time()))
    }
    #clean up and prep for next dose
    antigens_count[,SERIES_COMPLETE:=max(DOSE_COMPLETES_SERIES,na.rm = FALSE),by=.(STUDY_ID, ANTIGEN)]
    #store skipped - series completed previously, no further dose evaluations needed
    skipped_group<-antigens_count[SERIES_COMPLETE==TRUE&DOSE_COUNTER==0&is.na(VALID)]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=paste0("Series completed previously with ",counter," valid doses."),COUNTER=counter)]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if (verbose) message("Not evaluated since series already completed: ", nrow(skipped_group))
    }
    #report out
    if(verbose){message(paste0("Administrations validated dose ",counter,": ",antigens_count[VALID==TRUE&DOSE_COUNTER==counter,.N]))}
    
    if(verbose){message(paste0("Total administrations remaining: ",antigens_count[is.na(VALID),.N]))}
    antigens_valid[[counter]]<-antigens_count[VALID==TRUE]
    antigens_count<-antigens_count[is.na(VALID)]
    if(verbose){message(paste0("Completed Dose ",counter," evaluation @ ",Sys.time()))}
    gc()#clean up
    ##### DOSE 7+ checking ######
    #TO get to this point all previous doses would be valid and only matters for tetanus for re-dosing q10 years and yearly flu
    counter<-counter+1
    antigens_count[ANTIGEN=="TETANUS",`:=`(DOSE_COUNTER=as.numeric(ADMIN_COUNTER),VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(11)),NEXT_DOSE_MIN=DATE_GIVEN+yr_no_grace(10)
                                           ,NEXT_DOSE_RECOMMENDED=DATE_GIVEN+yr_no_grace(10))]
    #FLU is yearly, delayed is if more than 2 years between doses
    #Influenza - Yearly
    antigens_count[ANTIGEN=="INFLUENZA",`:=`(DOSE_COUNTER=as.numeric(ADMIN_COUNTER),VALID=TRUE,DELAYED=(INTERVAL>yr_with_grace(2))
                                                                                       ,NEXT_DOSE_RECOMMENDED=next_sept1(DATE_GIVEN))]
                                                                                                                     
    #COVID
    antigens_count[ANTIGEN=="COVID",`:=`(DOSE_COUNTER=as.numeric(ADMIN_COUNTER),VALID=TRUE)]
    #SEAONALITY OF DOSE ADMINISTER (e.g which Flu season was this considered season X is from 8/1/x-7/31/x+1, therefore if <8 will -1 the year to get X)
    antigens_count[ANTIGEN %in% SEASONAL_ANTIGENS,`:=`(MONTH_GIVEN=month_from_date(DATE_GIVEN),YEAR_GIVEN=year_from_date(DATE_GIVEN))]
    antigens_count[ANTIGEN == "INFLUENZA",FLU_SEASON := ifelse(MONTH_GIVEN > 7L, YEAR_GIVEN, YEAR_GIVEN - 1L)]
    antigens_count[ANTIGEN == "RSV",RSV_SEASON := ifelse(MONTH_GIVEN > 7L, YEAR_GIVEN, YEAR_GIVEN - 1L)]
    antigens_count[ANTIGEN == "COVID",COVID_SEASON := ifelse(MONTH_GIVEN > 7L, YEAR_GIVEN, YEAR_GIVEN - 1L)]
    #store skipped - series completed previously, no further dose evaluations needed
    skipped_group<-antigens_count[is.na(VALID)==TRUE]
    if(nrow(skipped_group)> 0L){
      skipped_group[,`:=`(NOTES=paste0("Invalid at endstep."))]
      skipped_list[[length(skipped_list) + 1L]]<-skipped_group
      antigens_count<-antigens_count[!skipped_group,on="TABLE_INDEX"]
      if (verbose) message("Not evaluated since series already completed: ", nrow(skipped_group))
    }
    #report out
    if(verbose){message(paste0("Administrations validated dose ",7,"+: ",antigens_count[VALID==TRUE,.N]))}
    if(verbose){message(paste0("Total administrations discarded: ",sum(vapply(skipped_list, nrow, integer(1L)))))}
    if(verbose){message(paste0("Total administrations remaining: ",antigens_count[is.na(VALID),.N]))}
    antigens_valid[[counter]]<-antigens_count[VALID==TRUE]
    antigens_count<-antigens_count[is.na(VALID)]
    if(verbose){message(paste0("Completed evaluation @ ",Sys.time()))}
    gc()#clean up
    #build return lists
    antigens_complete<-data.table::rbindlist(lapply(antigens_valid, function(dt) if(nrow(dt) > 0) dt else NULL),fill = TRUE)
    skipped <- data.table::rbindlist(lapply(skipped_list, function(dt) if(nrow(dt) > 0) dt else NULL),fill = TRUE)
    #report out
    if(verbose){message(paste0("Total administrations discarded: ",skipped[,.N]))}
    if(verbose){message(paste0("Total administrations ignrored: ",antigens_count[is.na(VALID),.N]))}
    if(verbose){message(paste0("Total administrations validated: ",antigens_complete[,.N]))}
    gc()
    ##### HANDLE OF PPV2: NOTES - TO BE ADDED within Immunocompromised build#####
    # Among persons aged ≥2 years with medical indications to receive both PCV13 and PPSV23 in a series, 
    # including adults aged ≥65 years with immunocompromising conditions, functional or anatomic asplenia, cochlear implants, or cerebrospinal fluid leaks, 
    # a dose of PPSV23 should be given ≥8 weeks after a dose of PCV13.
    # If a dose of PPSV23 is inadvertently given earlier than the recommended interval, the dose need not be repeated.
    # 1. do not need to check intervals for PPV23
    # 2. do need to check minimum age (24 months) (Immunize.org - does not generate an immune response)
    # 3. Do not need to redose PCV13 after PPV23 (ACIP guidance for adults)
    # 4. SERIES COMPLETE is dependent on reasons for dosing immunocompromised

    ##### REBUILD IMMUNIZATIONS GIVEN #####
    immunizations <- unique(
      antigens_complete[, .(STUDY_ID, CVX, DATE_GIVEN, PRODUCT, AGE_IMM_GIVEN)],
      by = c("STUDY_ID", "CVX", "DATE_GIVEN")
    )
    
    # Add one column per antigen in antigens_list
    for (ant in antigens_list) {
      immunizations[, (ant) := CVX %in% cvx[[ant]] ]
    }
    #finalize the list for return
    valid_output <- list(
      immunizations    = immunizations,
      antigens         = antigens_complete[,..ANTIGEN_RETURN_COLUMNS],
      skipped_antigens = skipped[,..INVALID_DOSE_COLUMNS]
    )
    #output messaging
    if(verbose){message(paste0("Completed validate_immunizations of immunizations @ ",Sys.time()))}
    if(verbose){message(paste0("Total administrations discards: ",skipped[,.N]))}
    if(verbose)message(paste0("Output returned includes 3 large lists:"))
    if(verbose)message(paste0("--immunizations list contains product level data"))      
    if(verbose)message(paste0("--antigens list contains dose level data"))      
    if(verbose)message(paste0("--skipped_antigens list contains all administrations that were not given at valid times"))
    return(valid_output)
}
