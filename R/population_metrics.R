#' R script to evaluate immunization reports against HEDIS and other measures.
#' 
#' @param x SOURCE_NAME.
#' @param y .
#' @return The sum of \code{x} and \code{y}.
#' @examples


population_metrics<-function(patients,antigens,date_of_birth_column_name='DOB',study_id_column_name='STUDY_ID'){
  #SUPPORT FUNCTIONS
  max0<-function(x,na.rm = TRUE){as.numeric(if(!is.infinite(suppressWarnings(max(x,na.rm=TRUE)))){max(x,na.rm=TRUE)}else {0})}
  yr_no_grace <-function(x) {yr_with_grace(x,grace=0)}
  yr_with_grace<-function(x,grace=4) {floor(365.25*x)-grace}
  
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
  return(patients)
}
