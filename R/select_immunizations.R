.datatable.aware = TRUE

#' This function filters and processes immunization data by mapping antigens to CVX codes, ensuring that the correct immunization data is processed and returned.
#' @param immunization_data DataTable of immunizations with STUDY_ID, PRODUCT, and DATE_GIVEN or AGE_IMM_GIVEN.
#' @param lim_cvx_map Local-to-CVX mapping, default is NA. Used if 'CVX' is not present in the immunization data.
#' @param antigen_to_eval List of antigens to evaluate. Default is 'ALL', meaning all antigens are evaluated.
#' @param verbose Boolean indicating whether to print log messages. Default is TRUE.
#' @param date_of_birth_column_name Name of the column containing date of birth information. Default is 'DOB'.
#' @param study_id_column_name Name of the column containing study ID information. Default is 'STUDY_ID'.
#' @param immunization_code_column_name Name of the column containing CVX codes. Default is 'CVX'.
#' @param local_immunization_identifier_column_name Name of the column containing local immunization identifiers. Default is 'LIM'.
#' @param immunization_product_name_column_name Name of the column containing immunization product names. Default is 'PRODUCT'.
#' @param immunization_date_given_column_name Name of the column containing the date the immunization was given. Default is 'DATE_GIVEN'.
#' @param age_at_immunization_column_name Name of the column containing age at immunization. Default is 'AGE_IMM_GIVEN'.
#' @return Returns a refined/limited immunization data with a column PROCESSED added. 

select_immunizations<-function(immunization_data,lim_cvx_map=NA,antigen_to_eval=c('ALL'),verbose=TRUE
                               ,date_of_birth_column_name='DOB',study_id_column_name='STUDY_ID'
                               ,immunization_code_column_name='CVX',local_immunization_identifier_column_name='LIM'
                               ,immunization_product_name_column_name='PRODUCT'
                               ,immunization_date_given_column_name='DATE_GIVEN',age_at_immunization_column_name='AGE_IMM_GIVEN'
                               ){
  
    ##### Validate and rename column names####
  #create local copy
  immunization_data <- data.table::copy(immunization_data)
  # Standard column mapping
  standard_cols <- list(
    STUDY_ID = study_id_column_name,
    DOB = date_of_birth_column_name,
    PRODUCT = immunization_product_name_column_name
  )
  
  # Loop through standard columns and rename if needed
  for (std_name in names(standard_cols)) {
    col_name <- standard_cols[[std_name]]
    if (!col_name %in% colnames(immunization_data)) {
      stop(paste0("Column '", col_name, "' required but not found in immunization_data."))
    }
    if (col_name != std_name) {
      setnames(immunization_data, col_name, std_name)
    }
  }
  
  # Handle CVX / LIM logic
  if (!immunization_code_column_name %in% colnames(immunization_data)) {
    if (is.null(lim_cvx_map)) {
      stop("Column 'CVX' is required unless a LIM-to-CVX map is supplied.")
    } else if (!local_immunization_identifier_column_name %in% colnames(immunization_data)) {
      stop("Column 'LIM' is required if CVX is missing and a LIM-to-CVX map is supplied.")
    } else if (local_immunization_identifier_column_name != 'LIM') {
      setnames(immunization_data,local_immunization_identifier_column_name,"LIM")
    }
  } else if (immunization_code_column_name != 'CVX') {
    setnames(immunization_data,immunization_code_column_name,"CVX")
  }
  
  # Handle DATE_GIVEN / AGE_IMM_GIVEN
  if (!immunization_date_given_column_name %in% colnames(immunization_data) &&
      !age_at_immunization_column_name %in% colnames(immunization_data)) {
    stop("At least one of 'DATE_GIVEN' or 'AGE_IMM_GIVEN' must exist in immunization_data")
  }
  
  if (immunization_date_given_column_name %in% colnames(immunization_data) &&
      immunization_date_given_column_name != 'DATE_GIVEN') {
    setnames(immunization_data,immunization_date_given_column_name,'DATE_GIVEN')
  }
  
  if (age_at_immunization_column_name %in% colnames(immunization_data) &&
      age_at_immunization_column_name != 'AGE_IMM_GIVEN') {
    setnames(immunization_data,age_at_immunization_column_name,'AGE_IMM_GIVEN')
  }
  
  # Decide which to use for downstream processing
  date_switch <- if ("DATE_GIVEN" %in% colnames(immunization_data)) "DATE_GIVEN" else "AGE_IMM_GIVEN"
  
  #### Load CVX mapping resource file ####
  cvx_files <- list.files(system.file("extdata", package = "VaccineValidateR"), full.names = TRUE)
  cvx_map <- data.table::rbindlist(lapply(cvx_files, data.table::fread))
  if (verbose) message("CVX-antigen map loaded. Mapping immunization data...")
  #add the CVX to map to the data if not already present
  if(!"CVX" %in% colnames(immunization_data)){
    if (verbose) message("Using locally supplied Local Code to CVX map. This may take some time.")
    immunization_data[lim_cvx_map, CVX := i.CVX, on = "LIM"]#map on the CVX from the LIM_CVX_MAP
  }
  
 
  ##### APPLY ANTIGEN-CVX MAPS ####
  #base AG maps from incoming data is mapped antigen or CVX code - do not need special cases here as there is not the need for expansive data evaluation
  tetanus_cvx <-cvx_map[grepl('TD|DT',ANTIGEN,ignore.case=TRUE)]$CVX
  polio_cvx <- cvx_map[grepl('IPV|OPV|POLIO',ANTIGEN,ignore.case=TRUE)]$CVX
  hib_cvx  <- cvx_map[grepl('HIB',ANTIGEN,ignore.case=TRUE)]$CVX
  pcv_cvx  <- cvx_map[grepl('PREVNAR|PCV',ANTIGEN,ignore.case=TRUE)]$CVX
  ppv23_cvx  <- cvx_map[grepl('PPV23',ANTIGEN,ignore.case=TRUE)]$CVX
  hepb_cvx  <- cvx_map[grepl('HEPB',ANTIGEN,ignore.case=TRUE)]$CVX
  hepa_cvx  <- cvx_map[grepl('HEPA',ANTIGEN,ignore.case=TRUE)]$CVX
  mmr_cvx  <- cvx_map[grepl('MMR',ANTIGEN,ignore.case=TRUE)]$CVX
  vzv_cvx  <- cvx_map[grepl('VARICELLA',ANTIGEN,ignore.case=TRUE)]$CVX
  hpv_cvx  <- cvx_map[grepl('HPV',ANTIGEN,ignore.case=TRUE)]$CVX
  mcv_cvx  <- cvx_map[grepl('MCV|MENACPOLY',ANTIGEN,ignore.case=TRUE)]$CVX
  menb_cvx  <- cvx_map[grepl('MENB',ANTIGEN,ignore.case=TRUE)]$CVX
  flu_cvx  <- cvx_map[grepl('INFLUENZA',ANTIGEN,ignore.case=TRUE)]$CVX
  live_flu_cvx<-  cvx_map[grepl('NASALINFLUENZA',ANTIGEN,ignore.case=TRUE)]$CVX
  rota_cvx  <- cvx_map[grepl('ROTA',ANTIGEN,ignore.case=TRUE)]$CVX
  rsv_cvx <- cvx_map[grepl('RSV',ANTIGEN,ignore.case=TRUE)]$CVX
  covid_cvx  <- cvx_map[grepl('COVID19',ANTIGEN,ignore.case=TRUE)]$CVX
  oral_polio_cvx <- cvx_map[grepl('OPV',ANTIGEN,ignore.case=TRUE)]$CVX
  yellowfever_cvx<-cvx_map[grepl('YELLOWFEVER',ANTIGEN,ignore.case=TRUE)]$CVX
  smallpox_cvx<-cvx_map[grepl('SMALLPOX',ANTIGEN,ignore.case=TRUE)]$CVX
  live_non_enteral_cvx<-c(mmr_cvx,vzv_cvx,live_flu_cvx,yellowfever_cvx,smallpox_cvx)#group all the live virus CVX codes
  
  ##### MERGE DATASETS AND APPLY ANTIGEN MAP TO IMMUNIZATIONS####
  #CONTAINS ANTIGEN VARIABLE ON EACH VACCINE
  immunization_data[, `:=`(TETANUS=CVX %in% tetanus_cvx
                           , POLIO=CVX %in% polio_cvx
                           , HIB=CVX %in% hib_cvx
                           , PCV=CVX %in% pcv_cvx
                           , PPV23=CVX %in% ppv23_cvx
                           , HEPB=CVX %in% hepb_cvx
                           , HEPA=CVX %in% hepa_cvx
                           , MMR=CVX %in% mmr_cvx
                           , VZV=CVX %in% vzv_cvx
                           , HPV=CVX %in% hpv_cvx
                           , MCV=CVX %in% mcv_cvx
                           , MENB=CVX %in% menb_cvx
                           , INFLUENZA=CVX %in% flu_cvx
                           , ROTA=CVX %in% rota_cvx
                           , COVID=CVX %in% covid_cvx
                           , LAIV=CVX %in% live_flu_cvx
                           , OPV=CVX %in% oral_polio_cvx
                           , RSV=CVX %in% rsv_cvx
                           , IS_LIVE_NON_ENTERAL = CVX %in% live_non_enteral_cvx
                           
  )]
  if (verbose) message("CVX codes employed to determine antigens within immunizations. Merging data. Please Wait...")
  
  #isolate antigens of interest in order to discard extraneous data
  antigens_list<-c('POLIO','HIB','PCV','HEPA','HEPB','HPV','MCV','MENB','TETANUS','ROTA','COVID','RSV','MMR','VZV','INFLUENZA')
  #if there is a live vaccine in antigens_to_eval OR antigens_to_eval is ALL we need to keep all like vaccines
  if(!"ALL" %in% antigen_to_eval){
      antigens_list<-intersect(antigens_list, antigen_to_eval)
  }
  #ensure all live virus data is kept if any potential Live data in data set.
  if(any(c('MMR','VZV','INFLUENZA') %in% antigens_list)==TRUE){
    antigens_list <- c(antigens_list, 'MMR','VZV','INFLUENZA')
  }
  antigens_list<-unique(antigens_list)#force only 1/antigen
  
  #create a list of immunizations based on what we care about
  rows_to_keep <- rowSums(immunization_data[, ..antigens_list])>0
  immunization_data <- immunization_data[rows_to_keep]
  if(verbose) message("Content restricted to antigens of interest. Finalizing content.")
  #rebuild the immunization list 
  #add columns of interest
  if(!"GIVEN_STATUS" %in% colnames(immunization_data)){immunization_data[,GIVEN_STATUS:=1]}
  if(!"ADMIN_LOCATION" %in% colnames(immunization_data)){immunization_data[,ADMIN_LOCATION:=NA]}
  immunization_data[,PROCESSED:=TRUE]#set a marker to prevent reprocessing of data
  immunization_data<-data.table::as.data.table(unique(immunization_data,by=c('STUDY_ID','CVX',date_switch)))
  
  if(verbose) message("This table can be sent to validation function directly and will not be reprocessed.")
  if(verbose) message("Do not change 'Processed' column. N=",nrow(immunization_data))
  #add columns of interest
  return(immunization_data)
}