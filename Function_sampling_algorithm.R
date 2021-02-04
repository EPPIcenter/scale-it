##############################
# Sampling algorithm
###############################

library(lubridate)
library(stringr)
library(data.table)
library(tidyverse)

#' Sampling EMR data for SCALE-IT serosurveillance 
#'
#' @param labs Lab test electronic medical record dataframe 
#' @param encounters encounters/diagnoses electronic medical record dataframe 
#' @param thresh  threshold for maximum number of samples (n.b. NOT exact number
#'   samples returned)
#' @param proportion_longitudinal proportion of data which are from patients
#'   previously sampled
#' @param LISdat Laboratory information system dataframe
#' @param prev_sampled Dataframe of previously sampled data - first column is
#'   date of sample collection, second column is MRN
#'
#' @return
#' @export
#'
#' @examples
gen_list <-
  function(labs,
           encounters,
           thresh,
           proportion_longitudinal,
           LISdat,
           LIS_startdate,
           LIS_enddate,
           prev_sampled) {
    
    
    
    # internal function used later
    substrRight <- function(x, n) {
      substr(x, nchar(x) - n + 1, nchar(x))
    }
    
    
    
    age_bins <- c(0, 20,40, 60,80,100)
    
    
    labs <- labs %>% mutate(
      age_groups = case_when(
        PAT_AGE >= 0  & PAT_AGE < 20 ~ "0-19",
        PAT_AGE >= 20  &
          PAT_AGE < 40 ~ '20-39',
        PAT_AGE >= 40  &
          PAT_AGE < 60 ~ '40-59',
        PAT_AGE >= 60  &
          PAT_AGE < 80 ~ '60-79',
        PAT_AGE >= 80  &
          PAT_AGE < 110 ~ '80+'
      )
    )
    
    
    labs_merge<- merge(labs, encounters, by = "PAT_MRN_ID")
    
   
    # Identify sentinel population in data 
    sentinel <-
      labs_merge[grepl("Z34", labs_merge[["CURRENT_ICD10_LIST"]]),]
    
    
    labs_sf <-
      labs[labs$zcta >= 94100 & labs$zcta < 94199 | labs$zcta == 99997, ]
    
    # Subset to chemistries
    chemistries_sf <-
      labs_sf[grepl("METABOLIC", labs_sf[["DESCRIPTION"]]) |
                grepl("CREATININE", labs_sf[["DESCRIPTION"]]) |
                grepl("ELECTROLYTES", labs_sf[["DESCRIPTION"]]) |
                grepl("ANTIBODIES", labs_sf[["DESCRIPTION"]]) |
                grepl("ANTIBODY", labs_sf[["DESCRIPTION"]]) |
                grepl("SYPHILIS", labs_sf[["DESCRIPTION"]]) |
                grepl("RUBELLA", labs_sf[["DESCRIPTION"]]) |
                grepl("SERUM / PLASMA", labs_sf[["DESCRIPTION"]]) , ]
    
    
    sentinel_sf <-
      chemistries_sf[chemistries_sf$PAT_MRN_ID %in% sentinel$PAT_MRN_ID,]#sentinel pops living in SF
    
    # select samples from emergency department 
    ed_chemistries_sf<- chemistries_sf[grepl("EMERGENCY", chemistries_sf[["DEPT"]]),]
    
    # subset to outpatient
    op_chemistries_sf <- chemistries_sf %>% filter(PAT_TYPE == "OP")
    
    # Adult outpatient chemistries
    adult_op_chemistries_sf <-
      op_chemistries_sf %>% filter(PAT_AGE > 17)
    
    # ALL under 18 chemistries (IP and OP)
    u18_chemistries_sf <- chemistries_sf %>% filter(PAT_AGE < 18)
    
    # ensure no repeated sentinel populations
    sentinel_sf <- distinct(sentinel_sf,  PAT_MRN_ID, .keep_all = TRUE)
    
    if (length(unique(sentinel_sf$PAT_MRN_ID)) > thresh * proportion_sentinel) {
      sentinel_sf <-
        dplyr::sample_n(
          sentinel_sf,
          size = thresh * proportion_sentinel,
          replace = F
        )
    }
    # Generate sample list containing  outpatient adults, under 18s, and sentinel
    sample_list <-
      rbind(adult_op_chemistries_sf, ed_chemistries_sf,u18_chemistries_sf, sentinel_sf)
    
    
    # Identify MRNS sampled previously
    prev_sampled_match<-prev_sampled[prev_sampled$PAT_MRN_ID %in% sample_list$PAT_MRN_ID == TRUE, c("PAT_MRN_ID", "SPEC_COLL_DATE", "SPEC_COLL_TIME")]
    sample_match<-unique(sample_list[sample_list$PAT_MRN_ID %in% prev_sampled[, "PAT_MRN_ID"] == TRUE,],  by = c("PAT_MRN_ID"))[, c("PAT_MRN_ID", "SPEC_COLL_DTTM")]
    prev_MRNs<-merge(sample_match, prev_sampled_match, by = "PAT_MRN_ID")
    prev_MRNs$TIME_DIFF<- as.Date(prev_MRNs$SPEC_COLL_DTTM) - as.Date(prev_MRNs$SPEC_COLL_DATE, format = "%m/%d/$Y")
    
    # Identify MRNS tested within last 30 days and exclude - currently no longitudinal MRNs so commented out
    recent_MRNS <-  prev_MRNs %>% filter( TIME_DIFF < 30)
    
    #longitudinal MRNs
    longitudinal_MRNs  <- prev_MRNs %>% filter( TIME_DIFF >= 30)
    longitudinal_MRNs<-distinct(longitudinal_MRNs, PAT_MRN_ID, .keep_all = TRUE)
    
    
    # Check how many longitudinal samples we have and, if it exceeds desired proportion, sample to that proportion
    if (length(unique(longitudinal_MRNs$PAT_MRN_ID)) > thresh * proportion_longitudinal) {
      longitudinal_MRNs <-
        dplyr::sample_n(
          longitudinal_MRNs,
          size = thresh * proportion_longitudinal,
          replace = F
        )
    }
    
    longitudinal_MRNs<-sample_list %>% filter (PAT_MRN_ID %in% longitudinal_MRNs$PAT_MRN_ID ==TRUE)
    longitudinal_MRNs<-distinct(longitudinal_MRNs, PAT_MRN_ID, .keep_all = TRUE)
    
    sample_list <-
      rbind(sample_list[sample_list$PAT_MRN_ID %in% recent_MRNS$PAT_MRN_ID == FALSE,], sentinel_sf, longitudinal_MRNs)
    
    sample_list<-distinct(sample_list, PAT_MRN_ID, .keep_all = T)
    
    
    ##############################################
    ## Stratification and sampling by zip and age
    ##############################################
    
    
    #if total sampled are below threshold, save all, if not sample from full dataset by age and zip
    
    # count overall samples - remembering sentinel has already been merged with the sample_list above 
    overall_samples<-length(unique(sample_list$PAT_MRN_ID))
    
    if (overall_samples <= thresh) {
      
      
    } else{
      
      
      
      #stratify by zip
      zips <- as.data.frame(table(sample_list$zcta))
      
      colnames(zips)<-c("zip", "Freq")
      zip_info<-read.csv("C:/Users/ir515/Box/EMR data/zipcode_info.csv")
      zip_info$target_zip <-as.integer(thresh*zip_info$proportion_population)
      zip_info$zip<-zip_info$zcta
      
      zips<-right_join(zip_info, zips, by ="zip")
      zips$diff <- zips$Freq - zips$target_zip
      zips[which(zips$zip == 99997),"diff"]<- thresh/length(unique(zips$zip))
      
      
      
      #stratify by age
      as.data.frame(table(sample_list$age_groups))
      
      # set target sample size to have equal numbers in age bin
      target_agebin <- as.integer(thresh / (length(age_bins)-1))
      
      #see if any age groups are over the target
      ages <- as.data.frame(table(sample_list$age_groups))
      
      ages_thin <- ages[ages$Freq > target_agebin, "Var1"]
      ages$diff <- ages$Freq - target_agebin
      
      #if a group is over, sample the target number from that group, but only from
      #over-represented zipcodes
      zips_to_sample <- zips %>% filter(diff > 0)
      ages_to_sample <- ages %>% filter(diff > 0)
      
      
      sample_data <- function(samplezip, sampleage) {
        tmp_age <- ages %>% filter(Var1 == sampleage)
        age_weight <- tmp_age$diff / sum(ages_to_sample$diff)
        tmp_df <-
          sample_list %>% filter(zcta == samplezip, age_groups == sampleage)
        tmp_zip <- zips %>% filter(zip == samplezip)
        if (round((tmp_zip$diff) * age_weight) > nrow(tmp_df)) {
          size <- nrow(tmp_df)
        } else{
          size <- round((tmp_zip$diff) * age_weight)
        }
        
        exclude <- dplyr::sample_n(tmp_df, size = size, replace = F)
        
        return(exclude)
      }
      
      
      exclude<-list()
      
      for(i in 1:length(ages_thin)){
        res <-
          lapply(zips_to_sample$zip, sample_data, sampleage = as.character(ages_thin[i]))
        res <- bind_rows(res)
        exclude[[i]]<-res
        
      }
      exclude<-bind_rows(exclude)
      
      
      thinned_sample <-
        sample_list[sample_list$PAT_MRN_ID %in% exclude$PAT_MRN_ID == FALSE,]
      
      
      # check for numbers removed
      nrow(sample_list)
      nrow(thinned_sample)
      
      #check distribution - note age weighting is still off
      
      table(thinned_sample$zcta)
      table(thinned_sample$age_groups)
      
      
      ####THE SAMPLES FROM THIS SHOULD THEN BE combined WITH sentinel_sf AND longitudinal_MRNs AS THEY ARE BELOW)#####
      sample_list <- rbind(thinned_sample, sentinel_sf, longitudinal_MRNs)
      
      
    }
    #Remove any duplicate MRNs
    sample_list <- distinct(sample_list, PAT_MRN_ID, .keep_all = TRUE)
    
    
    sample_list$PAT_MRN_ID<-as.character(sample_list$PAT_MRN_ID)
    
    #Restrict to UCSF labs
    sample_list <- sample_list %>% filter(RESULTING_LAB == "UCSF LAB" |
                                            RESULTING_LAB == "UCSF MOUNT ZION CLINICAL LABS" )
    
    #Make sure that there are no indivduals also tested for COVID 
    LISdat<-LISdat %>% filter( as.Date(LISdat$DateTimeCollected) > LIS_startdate &  as.Date(LISdat$DateTimeCollected) < LIS_enddate )
    testedforcovid <- LISdat[grepl("COV19", LISdat[["TestName"]]) |
                               grepl("COVID-19", LISdat[["TestName"]])  , ]
    testedforcovid_admit <- testedforcovid[grepl("PENDING", testedforcovid[["TestName"]]),]
    testedforcovid_nonadmit <- testedforcovid[!(testedforcovid$PtNumber %in% testedforcovid_admit$PtNumber),]
    sample_list$PAT_MRN_ID<-as.character(sample_list$PAT_MRN_ID)
    
    sample_list <- sample_list[!(sample_list$PAT_MRN_ID %in% testedforcovid_nonadmit$PtNumber),]
    
    return(sample_list)
  }



    




