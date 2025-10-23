#' Evaluate Population-Level Immunization Metrics
#'
#' This function evaluates childhood and adolescent immunization coverage against HEDIS and Non-HEDIS metrics (UTD by Interval/Dose Calculations)
#' It calculates whether each patient meets specific vaccine coverage criteria (e.g., HPV, Tetanus, MCV, DTaP, PCV, Rota, Influenza) 
#' for both HEDIS and Non-Hedis (direct calculation) standards based on the provided antigen data at age 13 (HEDIS ADOL) and age 2 (HEDIS CIS).
#'
#' @param patients A \code{data.table} of patients with at least a study ID column.
#' @param antigens A \code{data.table} of validated antigen-level immunization records 
#'   (output from the validate_immunizations functions) including AGE_IMM_GIVEN, DOSE_COUNTER, and DOSE_COMPLETES_SERIES.
#' @param date_of_birth_column_name Name of the column in \code{patients} containing patient DOB. Default is 'DOB'.
#' @param study_id_column_name Name of the column in \code{patients} containing study ID. Default is 'STUDY_ID'.
#' @param verbose Logical indicating whether to print log messages. Default is \code{TRUE}.
#' @return The \code{patients} data.table with additional columns for each HEDIS and UTD metric, 
#'   including adolescent and young child measures (e.g., HEDIS_ADOL2, UTD_ADOL2, HEDIS_CIS10, UTD_CIS10).
#' @export

population_metrics<-function(patients,antigens,date_of_birth_column_name='DOB',study_id_column_name='STUDY_ID',verbose=TRUE){
  if (verbose) message("Starting with validation checks.")
  #validation step
  check_antigen_table(antigens)
  patients <- data.table::setDT(data.table::copy(patients))
  antigens <- data.table::setDT(data.table::copy(antigens))
  # Map study id columns to 'STUDY_ID' consistently
  if (!study_id_column_name %in% names(patients)) {
    stop("patients is missing study id column: ", study_id_column_name)
  }
  if (!study_id_column_name %in% names(antigens)) {
    stop("antigens is missing study id column: ", study_id_column_name)
  }
  if (study_id_column_name != "STUDY_ID") {
    data.table::setnames(patients, study_id_column_name, "STUDY_ID")
    data.table::setnames(antigens, study_id_column_name, "STUDY_ID")
  }
  #antigens table is never modified, no need to create local copy
  if (verbose) message("Validation complete. Starting population-level immunization metric evaluation at ", Sys.time())
  # --- Adolescent Metrics ---
  if (verbose) message("Calculating HEDIS and UTD Adolescent metrics...")
  
  #HEDIS ADOLESCENT & UTD ADOLESCENT
  patients[,HPV_HEDIS_ADOL := (STUDY_ID %in% antigens[ANTIGEN=="HPV"&DOSE_COMPLETES_SERIES&AGE_IMM_GIVEN<yr_no_grace(13)]$STUDY_ID)]
  patients[,TETANUS_HEDIS_ADOL := (STUDY_ID %in% antigens[ANTIGEN=="TETANUS"&AGE_IMM_GIVEN>=yr_no_grace(10)&AGE_IMM_GIVEN<yr_no_grace(13)]$STUDY_ID)]
  #TETANUS IS NEVER COMPLETE
  patients[,TETANUS_UTD_ADOL := (STUDY_ID %in% antigens[ANTIGEN=="TETANUS"&
                                                          ((AGE_FIRST_DOSE>yr_no_grace(1)&AGE_FIRST_DOSE<=yr_no_grace(7)&AGE_IMM_GIVEN>=yr_no_grace(10)&DOSE_COUNTER>=4)|
                                                             (AGE_FIRST_DOSE>yr_no_grace(7)&AGE_IMM_GIVEN>=yr_no_grace(10)&DOSE_COUNTER>=3)|
                                                             (AGE_IMM_GIVEN>=yr_no_grace(10)&DOSE_COUNTER>=5)
                                                           )
                                                        &AGE_IMM_GIVEN<yr_no_grace(13)]$STUDY_ID)]
  patients[,MCV_HEDIS_ADOL := (STUDY_ID %in% antigens[ANTIGEN=="MCV"&AGE_IMM_GIVEN>=yr_no_grace(11)&AGE_IMM_GIVEN<yr_no_grace(13)]$STUDY_ID)]
  patients[,MCV_UTD_ADOL := (STUDY_ID %in% antigens[ANTIGEN=="MCV"&AGE_IMM_GIVEN>=yr_no_grace(10)&AGE_IMM_GIVEN<yr_no_grace(13)]$STUDY_ID)]
  patients[,HEDIS_ADOL2:=(HPV_HEDIS_ADOL&TETANUS_HEDIS_ADOL&MCV_HEDIS_ADOL)]
  patients[,UTD_ADOL2:=(HPV_HEDIS_ADOL&TETANUS_UTD_ADOL&MCV_UTD_ADOL)]
  
  #HEDIS YOUNG CHILD
  #Childhood Immunization Status (CIS): Combination 10 is a measure that tracks the percentage of children who have received all 10 recommended vaccines by their second birthday. 
  #These vaccines include DTaP, IPV, MMR, HiB, HepB, VZV, PCV, HepA, RV, and Influenza. 
  #This measure is used to assess the vaccination coverage of children and identify areas where immunization efforts may need to be improved. 
  #Added PCV and HIB SERIES COMPLETE evaluations as an alternative measure since HEDIS is inexact
  # --- Childhood CIS Metrics ---
  if (verbose) message("Calculating HEDIS and UTD Childhood CIS metrics...")
  
  patients[,DTAP_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="TETANUS"&DOSE_COUNTER==4&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,DTAP_UTD_CIS := (STUDY_ID %in% antigens[ANTIGEN=="TETANUS"&((AGE_FIRST_DOSE>=yr_no_grace(1)&DOSE_COUNTER==3)|DOSE_COUNTER==4)&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,IPV_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="POLIO"&DOSE_COUNTER==3&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,MMR_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="MMR"&DOSE_COUNTER==1&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,VZV_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="VZV"&DOSE_COUNTER==1&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,PCV_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="PCV"&DOSE_COUNTER==4&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,PCV_UTD_CIS := (STUDY_ID %in% antigens[ANTIGEN=="PCV"&DOSE_COMPLETES_SERIES&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,HIB_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="HIB"&DOSE_COUNTER==3&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,HIB_UTD_CIS := (STUDY_ID %in% antigens[ANTIGEN=="HIB"&DOSE_COMPLETES_SERIES&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,HEPB_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="HEPB"&DOSE_COUNTER==3&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,HEPA_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="HEPA"&DOSE_COUNTER==2&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  patients[,ROTA_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="ROTA"&DOSE_COMPLETES_SERIES&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]#uses ANTIGEN LOGIC
  patients[,INFLUENZA_HEDIS_CIS := (STUDY_ID %in% antigens[ANTIGEN=="INFLUENZA"&DOSE_COUNTER==2&AGE_IMM_GIVEN<yr_no_grace(2)]$STUDY_ID)]
  
  #Build population metric statements
  patients[,HEDIS_CIS2:=(DTAP_HEDIS_CIS&IPV_HEDIS_CIS&MMR_HEDIS_CIS&VZV_HEDIS_CIS&HIB_HEDIS_CIS&HEPB_HEDIS_CIS)]
  patients[,UTD_CIS2:=(DTAP_UTD_CIS&IPV_HEDIS_CIS&MMR_HEDIS_CIS&VZV_HEDIS_CIS&HIB_UTD_CIS&HEPB_HEDIS_CIS)]
  patients[,HEDIS_CIS3:=(HEDIS_CIS2&PCV_HEDIS_CIS)]
  patients[,UTD_CIS3:=(UTD_CIS2&PCV_UTD_CIS)]
  patients[,HEDIS_CIS7:=(HEDIS_CIS3&ROTA_HEDIS_CIS)]
  patients[,UTD_CIS7:=(UTD_CIS3&ROTA_HEDIS_CIS)]
  patients[,HEDIS_CIS10:=(HEDIS_CIS7&INFLUENZA_HEDIS_CIS&HEPA_HEDIS_CIS)]
  patients[,UTD_CIS10:=(UTD_CIS7&INFLUENZA_HEDIS_CIS&HEPA_HEDIS_CIS)]
  if (verbose) {
    message("Population-level metric definitions (final outputs):")
    
    # Childhood measures
    message("  • HEDIS_CIS2: By age 2, child received 4 DTaP, 3 Polio, 1 MMR, 1 Varicella, 3 Hib, and 3 HepB (count-based).")
    message("  • UTD_CIS3: By age 2, child completed: DTaP, Polio, MMR, Varicella, Hib, and HepB series (using completion rules).")
    
    message("  • HEDIS_CIS3: By age 2, child received 4 DTaP, 3 Polio, 1 MMR, 1 Varicella, 3 Hib, 3 HepB, and 4 PCV (count-based).")
    message("  • UTD_CIS3: By age 2, child completed: DTaP, Polio, MMR, Varicella, Hib, HepB, and PCV series (using completion rules).")
    
    message("  • HEDIS_CIS7: By age 2, child received 4 DTaP, 3 Polio, 1 MMR, 1 Varicella, 3 Hib, 3 HepB, 4 PCV, and Rota series completion (count-based).")
    message("  • UTD_CIS7: By age 2, child completed: DTaP, Polio, MMR, Varicella, Hib, HepB, PCV, and Rota (all series-completion based).")
    
    message("  • HEDIS_CIS10: By age 2, child received 4 DTaP, 3 Polio, 1 MMR, 1 Varicella, 3 Hib, 3 HepB, 4 PCV, Rota series, 2 Influenza, and 2 HepA (count-based).")
    message("  • UTD_CIS10: By age 2, child completed: DTaP, Polio, MMR, Varicella, Hib, HepB, PCV, Rota, Influenza (2 doses), and HepA (2-dose series).")
    
    # Adolescent measures
    message("  • HEDIS_ADOL2: By age 13, adolescent received 2 HPV, ≥1 Tetanus dose age 10–13, and ≥1 MCV dose age 11–13 (count-based).")
    message("  • UTD_ADOL2: By age 13, adolescent completed: HPV (series-based), is up to date on Tetanus (evaluated by series completion rules), and received MCV (≥1 dose age 10–13).")
    
    message("Completed evaluation of HEDIS and UTD population metrics at ", Sys.time())
  }
  # Store definitions as an attribute of patients
  population_metric_definitions <- list(
    HEDIS_CIS2  = "By age 2, child received 4 DTaP, 3 Polio, 1 MMR, 1 Varicella, 3 Hib, and 3 HepB (count-based).",
    UTD_CIS2    = "By age 2, up to date as recommended for: DTaP, Polio, MMR, Varicella, Hib, and HepB series (using completion rules).",
    HEDIS_CIS3  = "By age 2, child received 4 DTaP, 3 Polio, 1 MMR, 1 Varicella, 3 Hib, 3 HepB, and 4 PCV (count-based).",
    UTD_CIS3    = "By age 2, up to date as recommended for: DTaP, Polio, MMR, Varicella, Hib, HepB, and PCV series (using completion rules).",
    HEDIS_CIS7  = "By age 2, child received 4 DTaP, 3 Polio, 1 MMR, 1 Varicella, 3 Hib, 3 HepB, 4 PCV, and Rota series completion (count-based).",
    UTD_CIS7    = "By age 2, up to date as recommended for: DTaP, Polio, MMR, Varicella, Hib, HepB, PCV, and Rota (all series-completion based).",
    HEDIS_CIS10 = "By age 2, child received 4 DTaP, 3 Polio, 1 MMR, 1 Varicella, 3 Hib, 3 HepB, 4 PCV, Rota series, 2 Influenza, and 2 HepA (count-based).",
    UTD_CIS10   = "By age 2, up to date as recommended for: DTaP, Polio, MMR, Varicella, Hib, HepB, PCV, Rota, Influenza (2 doses), and HepA (2-dose series).",
    HEDIS_ADOL2 = "By age 13, adolescent received 2 HPV, ≥1 Tetanus dose age 10–13, and ≥1 MCV dose age 11–13 (count-based).",
    UTD_ADOL2   = "By age 13, adolescent up to date as recommended for: HPV (series-based), is up to date on Tetanus (evaluated by series completion rules), and received MCV (≥1 dose age 10–13)."
  )
  data.table::setattr(patients, "population_metric_definitions",population_metric_definitions)
  data.table::setattr(patients, "population_metrics", names(population_metric_definitions))
  
  if (verbose) message("Attached metric list as attribute: population_metrics")
  return(patients)
}
