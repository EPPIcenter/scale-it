###############################################
##  Code to generate stratified seroprevalence 
###############################################

# Script to generate posterior estimates of seroprevalence  stratified by
# Age, Sex, Ethnicity, Month of collection, Insurance and Neighborhood


library(rstan)
library(rstantools)
library(bayesplot)
library(readr)
library(tidyverse)
set.seed(1234)


# Set working directory
setwd("INSERT WORKING DIRECTORY PATH")

# read in data
samples_collected<-read_csv("~/Supplementary Table 2.csv")[,1:4]

# add data filtered to remove homeless individuals and samples which couldn't be geolocated to rooftop level
neighborhood_data<-read_csv("~/Supplementary Table 3.csv")[,1:3]


# function to summarize data by category
#' Generate posteriors
#'
#' @param category_names string naming the category used for stratification
#' @param data.samps  data frame with 3 or 4 columns: Class column for multiple stratification classification categories such as Age and Sex, column with unique names of groups for stratification (named Group if Class column is present, or as the same as category_names if not ), number of positive samples and total number of samples
#'
#' @return
#' @export
#'
#' @examples

gen_posteriors <- function(category_names, data.samps) {

  if(is.null(data.samps$Class) == FALSE){
    data.samps<-data.samps %>% filter(Class == category_names)
    
   test_neg_samps_by_category<-  data.samps$Count - data.samps$Positive
    n_samps_by_category<-data.samps$Count
    test_pos_samps_by_category <-data.samps$Positive
    names_by_category<-data.samps$Group
    
    
  } else{
    
    #remove groups with less than 50 samples for map 
    data.samps %>%
      mutate(category = get(category_names)) %>%
      filter(Count >50)-> category
    target_category<-category$category
    
    test_neg_samps_by_category<-  category$Count - category$Positive
    n_samps_by_category<-category$Count
    test_pos_samps_by_category <-category$Positive
    names_by_category<-category$category
  }
  
  
  sample_size_tab <-
    as.data.frame(cbind(names_by_category, n_samps_by_category))
  
  
  ## Fit the model
  fit_Stan_by_category <- stan(
    model_code = fit_Stan_by_category,
    data = list(
      ## This is data from the GLM_SRN (3 antigen: Spike + RBD + N) model
      ## This is ordered as (ELISA, Luminex): (+,+), (+,-), (-,+), (-,-)
      y_positive_controls_tested_on_both = c(55, 1, 0, 0),
      y_negative_controls_tested_on_both = c(0, 1, 0, 87),
      ## This is ordered as: +, -
      y_positive_controls_ELISA_only = c(59, 2),
      y_negative_controls_ELISA_only = c(1, 4),
      ## This is ordered as: +, -
      y_positive_controls_Luminex_only = c(72, 5),
      y_negative_controls_Luminex_only = c(0, 26),
      ##
      n_category = length(test_pos_samps_by_category),
      y_sample = test_pos_samps_by_category,
      n_sample = n_samps_by_category
    ),
    chains = 4,
    iter = 15000,
    thin = 20,
    init = 0,
    refresh = 0,
    control = list(adapt_delta = 0.99)
  )
  
  index <- c(7:(6 + length(names_by_category)))
  
  #Pull out posteriors and rename columns
  posterior_by_category <- as.matrix(fit_Stan_by_category)[, index]
  
  colnames(posterior_by_category) <- c(names_by_category)
  
  # p1 <-
  #   mcmc_hist(posterior_by_category, facet_args = list(ncol = 2)) + theme_minimal() + ggtitle(
  #     "Model allowing for conditional dependence in Se and Sp between ELISA and Luminex + SCALE-IT results"
  #   )
  # p2 <- mcmc_trace(posterior_by_category)
  # 
  
  mean <- summary(fit_Stan_by_category)$summary[index, "mean"]
  lowci <- summary(fit_Stan_by_category)$summary[index, "2.5%"]
  hici <- summary(fit_Stan_by_category)$summary[index, "97.5%"]
  res <- as.data.frame(cbind(mean, lowci, hici))
  
   res$strat_cat <- as.vector(names_by_category)
  res$count<- n_samps_by_category
  res$raw_seroprev<- test_pos_samps_by_category/n_samps_by_category
  return(list(posterior_by_category,res))
  
}


## Model adapted from: https://pubmed.ncbi.nlm.nih.gov/10802336/
fit_Stan_by_category <- '
data {

  // This is ordered as (ELISA, Luminex): (+,+), (+,-), (-,+), (-,-)
  int y_positive_controls_tested_on_both[4];
  int y_negative_controls_tested_on_both[4];

  // This is ordered as: +, -
  int y_positive_controls_ELISA_only[2];
  int y_negative_controls_ELISA_only[2];

  // This is ordered as: +, -
  int y_positive_controls_Luminex_only[2];
  int y_negative_controls_Luminex_only[2];


  int n_category;
  int y_sample[n_category];
  int n_sample[n_category];

}
parameters {

  real logit_Se_ELISA;
  real logit_Sp_ELISA;

  real logit_Se_Luminex;
  real logit_Sp_Luminex;

  real<lower=fmax(-(1.0-inv_logit(logit_Se_ELISA))*(1.0-inv_logit(logit_Se_Luminex)),-inv_logit(logit_Se_ELISA)*inv_logit(logit_Se_Luminex)), upper=fmin(inv_logit(logit_Se_ELISA)*(1.0-inv_logit(logit_Se_Luminex)), inv_logit(logit_Se_Luminex)*(1.0-inv_logit(logit_Se_ELISA)))> covariance_Se;
  real<lower=fmax(-(1.0-inv_logit(logit_Sp_ELISA))*(1.0-inv_logit(logit_Sp_Luminex)),-inv_logit(logit_Sp_ELISA)*inv_logit(logit_Sp_Luminex)), upper=fmin(inv_logit(logit_Sp_ELISA)*(1.0-inv_logit(logit_Sp_Luminex)), inv_logit(logit_Sp_Luminex)*(1.0-inv_logit(logit_Sp_ELISA)))> covariance_Sp;

  real<lower=0.0, upper=1.0> prev[n_category];

}
transformed parameters {

  real<lower=0.0, upper=1.0> Se_ELISA;
  real<lower=0.0, upper=1.0> Sp_ELISA;

  real<lower=0.0, upper=1.0> Se_Luminex;
  real<lower=0.0, upper=1.0> Sp_Luminex;

  simplex[4] p_positive_controls_tested_on_both;
  simplex[4] p_negative_controls_tested_on_both;

  real<lower=0.0, upper=1.0> Se_serial_overall;
  real<lower=0.0, upper=1.0> Sp_serial_overall;

  Se_ELISA = inv_logit(logit_Se_ELISA);
  Sp_ELISA = inv_logit(logit_Sp_ELISA);

  Se_Luminex = inv_logit(logit_Se_Luminex);
  Sp_Luminex = inv_logit(logit_Sp_Luminex);

  p_positive_controls_tested_on_both[1] = (Se_ELISA*Se_Luminex) + covariance_Se;
  p_positive_controls_tested_on_both[2] = (Se_ELISA*(1.0-Se_Luminex)) - covariance_Se;
  p_positive_controls_tested_on_both[3] = ((1.0-Se_ELISA)*Se_Luminex) - covariance_Se;
  p_positive_controls_tested_on_both[4] = ((1.0-Se_ELISA)*(1.0-Se_Luminex)) + covariance_Se;

  p_negative_controls_tested_on_both[1] = ((1.0-Sp_ELISA)*(1.0-Sp_Luminex)) + covariance_Sp;
  p_negative_controls_tested_on_both[2] = ((1.0-Sp_ELISA)*Sp_Luminex) - covariance_Sp;
  p_negative_controls_tested_on_both[3] = (Sp_ELISA*(1.0-Sp_Luminex)) - covariance_Sp;
  p_negative_controls_tested_on_both[4] = (Sp_ELISA*Sp_Luminex) + covariance_Sp;

  Se_serial_overall = Se_ELISA*Se_Luminex + covariance_Se;
  Sp_serial_overall = 1.0 - ((1.0-Sp_ELISA)*(1.0-Sp_Luminex)) - covariance_Sp;

}
model {
  real p_sample[n_category];

  // Test performance characteristics
  target += multinomial_lpmf(y_positive_controls_tested_on_both | p_positive_controls_tested_on_both);
  target += multinomial_lpmf(y_negative_controls_tested_on_both | p_negative_controls_tested_on_both);

  target += binomial_lpmf(y_positive_controls_ELISA_only[1] | sum(y_positive_controls_ELISA_only), Se_ELISA);
  target += binomial_lpmf(y_negative_controls_ELISA_only[2] | sum(y_negative_controls_ELISA_only), Sp_ELISA);

  target += binomial_lpmf(y_positive_controls_Luminex_only[1] | sum(y_positive_controls_Luminex_only), Se_Luminex);
  target += binomial_lpmf(y_negative_controls_Luminex_only[2] | sum(y_negative_controls_Luminex_only), Sp_Luminex);



  for(i in 1:n_category) {

		p_sample[i] = prev[i]*Se_serial_overall + (1-prev[i])*(1-Sp_serial_overall);
  	y_sample[i] ~ binomial(n_sample[i], p_sample[i]);
  }


}
'


samples_collected$Class
#samples_collected$month<- as.character(samples_collected$month)
age <- gen_posteriors("Age Group", data.samps = samples_collected)
sex<-gen_posteriors("Sex", data.samps = samples_collected)
ethn<-gen_posteriors("Race/Ethnicity", data.samps = samples_collected)
fin<-gen_posteriors("Insurance Type", data.samps = samples_collected)
month<-gen_posteriors("Month", data.samps = samples_collected)
hood<-gen_posteriors("Neighborhood", data.samps = neighborhood_data)


saveRDS(list(hood[[2]], "~/neighborhood_posteriors.RData")
write.csv(hood[[2]], "~/neighborhood_posteriors.csv")

posteriors <-
  list(
    age = age,
    sex = sex,
    ethn = ethn,
    month = month[[2]],
    fin=fin
  )
saveRDS(posteriors, file = "~/posteriors.RData")
write.csv(rbind(posteriors$age[[2]],posteriors$sex[[2]], posteriors$hosp[[2]], posteriors$ethn[[2]], posteriors$fin[[2]])[,c("strat_cat","count","raw_seroprev", "mean", "lowci", "hici")], "~/posteriors.csv")




#####################################
# testing differences between groups
####################################

sex_diff<- quantile((sex[[1]][,2]-sex[[1]][,1]), probs = c(0.025,0.5,0.975))


hosp_diff<- quantile(hosp[[1]][,2]-hosp[[1]][,1], probs = c(0.025,0.5,0.975))



calc_diff<-function(file){
  diff<-matrix(ncol = ncol(file[[1]]), nrow =nrow(file[[1]]) )
  for(i in 1:nrow(file[[1]])){
    x<- file[[1]][i,]
    diff[i,]<-outer(x,x, `-`)[3,]
    
  }
  return(diff)
}
age_diff<-calc_diff(age)
ethn_diff<-calc_diff(ethn)
fin_diff<-calc_diff(fin)

quantile(age_diff[,1], probs = c(0.025,0.5,0.975))
quantile(age_diff[,2], probs = c(0.025,0.5,0.975))
quantile(age_diff[,4], probs = c(0.025,0.5,0.975))
quantile(age_diff[,5], probs = c(0.025,0.5,0.975))


quantile(ethn_diff[,1], probs = c(0.025,0.5,0.975))
quantile(ethn_diff[,2], probs = c(0.025,0.5,0.975))
quantile(ethn_diff[,4], probs = c(0.025,0.5,0.975))
quantile(ethn_diff[,5], probs = c(0.025,0.5,0.975))


quantile(fin_diff[,1], probs = c(0.025,0.5,0.975))
quantile(fin_diff[,2], probs = c(0.025,0.5,0.975))
quantile(fin_diff[,3], probs = c(0.025,0.5,0.975))
quantile(fin_diff[,4], probs = c(0.025,0.5,0.975))
