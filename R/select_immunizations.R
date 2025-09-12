.datatable.aware = TRUE
#' Filter and Process Immunization Data
#'
#' Filters and processes immunization data by mapping antigens to CVX codes, ensuring that the correct immunization data is processed and returned.
#'
#' @param immunization_data A \code{data.table} of immunizations containing at minimum \code{STUDY_ID}, \code{PRODUCT}, and either \code{DATE_GIVEN} or \code{AGE_IMM_GIVEN}.
#' @param lim_cvx_map Optional \code{data.table} for local-code→CVX mappings. Default is \code{NULL}. Used if 'CVX' is not present in the immunization data.
#' @param antigen_to_eval Character vector of antigens to evaluate (e.g., \code{c('POLIO','MMR')}) or 'ALL' (case-insensitive).
#' @param verbose Logical indicating whether to print log messages. Default is \code{TRUE}.
#' @param date_of_birth_column_name Character string specifying the column containing date of birth. Default is 'DOB'.
#' @param study_id_column_name Character string specifying the column containing study ID. Default is 'STUDY_ID'.
#' @param immunization_code_column_name Character string specifying the column containing CVX codes. Default is 'CVX'.
#' @param local_immunization_identifier_column_name Character string specifying the column containing local immunization identifiers. Default is 'LIM'.
#' @param immunization_product_name_column_name Character string specifying the column containing immunization product names. Default is 'PRODUCT'.
#' @param immunization_date_given_column_name Character string specifying the column containing the date the immunization was given. Default is 'DATE_GIVEN'.
#' @param age_at_immunization_column_name Character string specifying the column containing age at immunization. Default is 'AGE_IMM_GIVEN'.
#' @details
#' Antigen selection has special handling for live, non-enteral vaccines 
#' (e.g., \code{MMR}, \code{VZV}, intranasal \code{INFLUENZA}). 
#' If any of these are included in \code{antigen_to_eval} or present in the 
#' dataset, all such live antigens are retained together. This ensures that 
#' live-virus scheduling rules (e.g., spacing between doses) are preserved 
#' during downstream validate_immunizations.
#' @examples
#' library(data.table)
#'
#' # Example immunization data
#' imm_dt <- data.table(
#'   STUDY_ID   = c(1, 2, 3, 3),
#'   DOB        = as.Date(c("2020-01-15", "2019-07-30", "2018-03-15", "2018-03-15")),
#'   PRODUCT    = c("MMR Vax", "DTaP", "Polio", "Varicella"),
#'   CVX        = c(03, 20, 10, 21),           # simplified CVX codes
#'   DATE_GIVEN = as.Date(c("2021-02-01", "2020-09-01", "2020-05-01", "2021-06-01"))
#' )
#'
#' # Select all immunizations
#' selected_all <- select_immunizations(imm_dt)
#' print(selected_all)
#'
#' # Select only POLIO and MMR antigens
#' selected_subset <- select_immunizations(imm_dt, antigen_to_eval = c("POLIO", "MMR"))#will also return Varicella(VZV) since any live nonenteral antigen forces all to be kept
#' print(selected_subset)
#' @return A \code{data.table} containing the refined and filtered immunization data with an added attribute \code{processed} indicating the table is ready.
#'   An attribute \code{antigens_used} (character vector) is attached, recording
#'   which antigens were evaluated.
#'
#' @export

select_immunizations<-function(immunization_data,lim_cvx_map=NULL,antigen_to_eval=c('ALL'),verbose=TRUE
                               ,date_of_birth_column_name='DOB',study_id_column_name='STUDY_ID'
                               ,immunization_code_column_name='CVX',local_immunization_identifier_column_name='LIM'
                               ,immunization_product_name_column_name='PRODUCT'
                               ,immunization_date_given_column_name='DATE_GIVEN',age_at_immunization_column_name='AGE_IMM_GIVEN'
                               ){
  
    ##### Validate and rename column names####
  #create local copy
  immunization_data <- data.table::setDT(data.table::copy(immunization_data))
  antigen_to_eval <- toupper(antigen_to_eval)
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
  
  #add CVX to immunization data if not already present
  if(!"CVX" %in% colnames(immunization_data)){
    if (verbose) message("Using locally supplied Local Code to CVX map. This may take some time.")
    # Harmonize types to character to avoid join misses
    immunization_data[, LIM := as.character(LIM)]
    lim_cvx_map[, LIM := as.character(LIM)]
    immunization_data[lim_cvx_map, CVX := i.CVX, on = "LIM"]#map on the CVX from the LIM_CVX_MAP
  }
  
  if (verbose) message("CVX codes employed to determine antigens within immunizations. Merging data. Please Wait...")
  
  #isolate antigens of interest in order to discard extraneous data
  antigens_list<-SYSTEM_ANTIGENS
  #if there is a live vaccine in antigens_to_eval OR antigens_to_eval is ALL we need to keep all like vaccines
  if(!"ALL" %in% antigen_to_eval){
      antigens_list<-intersect(antigens_list, antigen_to_eval)
  }
  #Hanld Live Virus Antigens
  if(any(c('MMR','VZV','INFLUENZA') %in% antigens_list)==TRUE){
    antigens_list <- c(antigens_list, 'MMR','VZV','INFLUENZA')#ensure all live antigens kept (non - enteral)
  }
  antigens_list<-unique(antigens_list)#force only 1/antigen
  #confirm that at least 1 antigen is being evaluated fromthe approved list
  if (length(antigens_list) == 0L) {
    if (verbose) message("No matching antigens requested; returning 0 rows.")
    return(immunization_data[0])#no processed column - prevents continuation with no data
  }
  #add the antigens that matter to the immunization_data list
  for (antigen in antigens_list) {
    # CVX = column of codes, cvx = global mapping list
    immunization_data[, (antigen) := CVX %in% cvx[[antigen]] ]
  }
  
  
  #create a list of immunizations based on what we care about existing
  rows_to_keep <- rowSums(immunization_data[, ..antigens_list], na.rm = TRUE)>0
  immunization_data <- immunization_data[rows_to_keep]
  if(verbose) message("Content restricted to antigens of interest. Finalizing content.")
  #rebuild the immunization list 
  #add columns of interest
  if(!"GIVEN_STATUS" %in% colnames(immunization_data)){immunization_data[,GIVEN_STATUS:=1]}
  if(!"ADMIN_LOCATION" %in% colnames(immunization_data)){immunization_data[,ADMIN_LOCATION:=NA]}
  #create the processed column
  immunization_data<-data.table::as.data.table(unique(immunization_data,by=c('STUDY_ID','CVX',date_switch)))
  if(verbose) message("This table can be sent to validate_immunizations function directly and will not be reprocessed.")
  if(verbose) message("Do not change 'processed' attribute. N=",nrow(immunization_data))
  #set metadata attributes
  data.table::setattr(immunization_data, "antigens_used", antigens_list)
  data.table::setattr(immunization_data, "processed", TRUE)
  #add columns of interest
  return(immunization_data)
}