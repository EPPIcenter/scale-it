###################################################################
## Script to estimate seroprevalence, weighted by the distribution
## of age and sex in sample.
## adjusted for sensitivty and specificity of test performance
###################################################################


library(rstan)
library(bayesplot)

set.seed(123)

## Model adapted from: https://pubmed.ncbi.nlm.nih.gov/10802336/
Stan_model_cond_covariance <- '
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
  
  //
  int n_categories;
  int y_sample[n_categories];
  int n_sample[n_categories];
  real wt[n_categories];

}
parameters {
  
  real logit_Se_ELISA;
  real logit_Sp_ELISA;
  
  real logit_Se_Luminex;
  real logit_Sp_Luminex;
  
  real<lower=fmax(-(1.0-inv_logit(logit_Se_ELISA))*(1.0-inv_logit(logit_Se_Luminex)),-inv_logit(logit_Se_ELISA)*inv_logit(logit_Se_Luminex)), upper=fmin(inv_logit(logit_Se_ELISA)*(1.0-inv_logit(logit_Se_Luminex)), inv_logit(logit_Se_Luminex)*(1.0-inv_logit(logit_Se_ELISA)))> covariance_Se;
  real<lower=fmax(-(1.0-inv_logit(logit_Sp_ELISA))*(1.0-inv_logit(logit_Sp_Luminex)),-inv_logit(logit_Sp_ELISA)*inv_logit(logit_Sp_Luminex)), upper=fmin(inv_logit(logit_Sp_ELISA)*(1.0-inv_logit(logit_Sp_Luminex)), inv_logit(logit_Sp_Luminex)*(1.0-inv_logit(logit_Sp_ELISA)))> covariance_Sp;
  
  real<lower=0.0, upper=1.0> prev;
  
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
 real p_sample[n_categories];
  
  // Test performance characteristics
  target += multinomial_lpmf(y_positive_controls_tested_on_both | p_positive_controls_tested_on_both);
  target += multinomial_lpmf(y_negative_controls_tested_on_both | p_negative_controls_tested_on_both);
  
  target += binomial_lpmf(y_positive_controls_ELISA_only[1] | sum(y_positive_controls_ELISA_only), Se_ELISA);
  target += binomial_lpmf(y_negative_controls_ELISA_only[2] | sum(y_negative_controls_ELISA_only), Sp_ELISA);
  
  target += binomial_lpmf(y_positive_controls_Luminex_only[1] | sum(y_positive_controls_Luminex_only), Se_Luminex);
  target += binomial_lpmf(y_negative_controls_Luminex_only[2] | sum(y_negative_controls_Luminex_only), Sp_Luminex);
  
  // Sero-prevalence
  
  for(i in 1:n_categories){
  
		p_sample[i] = prev*Se_serial_overall + (1-prev)*(1-Sp_serial_overall);
    target += binomial_lpmf(y_sample[i] | n_sample[i], 	p_sample[i])*wt[i];
  }

  
}
'

## make weighting matrix

sample_mat<-matrix(nrow=10,ncol=3)
#colnames(sample_mat)<-c("pos","total","weight")
# cat 1 - 0-19 M
sample_mat[1,]<-c( 7,120 , 0.508*0.15)
# cat 2 - 0-19 F
sample_mat[2,]<-c( 4,159 , 0.493*0.15)
# cat 3 - 20-39 M
sample_mat[3,]<-c(28 ,456 , 0.508*0.38)
# cat 4 - 20-39 F
sample_mat[4,]<-c(22, 801 , 0.493*0.38)
# cat 5 - 40-59 M
sample_mat[5,]<-c(47 ,711 , 0.508*0.253)
# cat 6 - 40-59 F
sample_mat[6,]<-c(25 ,610, 0.493*0.253)
# cat 7 - 60-79 M
sample_mat[7,]<-c( 26,712,0.58*0.173)
# cat 8 - 60-79 F
sample_mat[8,]<-c( 15,673 , 0.493*0.173)
# cat 9 - 80+ M
sample_mat[9,]<-c( 11,232, 0.58*0.043)
# cat 10 - 80+ F
sample_mat[10,]<-c( 6,248 , 0.493*0.043)

## Fit the model
fit_Stan_cond_covariance <- stan(
  model_code=Stan_model_cond_covariance,
  data=list(
    ## This is data from the GLM_SRN (3 antigen: Spike + RBD + N) model
    ## This is ordered as (ELISA, Luminex): (+,+), (+,-), (-,+), (-,-)
    y_positive_controls_tested_on_both=c(55,1,0,0),
    y_negative_controls_tested_on_both=c(0,1,0,87),
    ## This is ordered as: +, -
    y_positive_controls_ELISA_only=c(59,2),
    y_negative_controls_ELISA_only=c(1,4),
    ## This is ordered as: +, -
    y_positive_controls_Luminex_only=c(72,5),
    y_negative_controls_Luminex_only=c(0,26),
    ## This is ordered as: 10 x 3: +, SUM, weight
    ## TO UPDATE
    y_sample= c(sample_mat[,1]), n_sample= c(sample_mat[,2]), wt = c(sample_mat[,3]), n_categories = nrow(sample_mat),
    chains=4, iter=50000, thin=20, init=0, refresh=0, control=list(adapt_delta=0.9999)))

print(fit_Stan_cond_covariance, digits_summary=4)


mcmc_hist(fit_Stan_cond_covariance, pars=c("Se_ELISA", "Sp_ELISA", "Se_Luminex", "Sp_Luminex", "covariance_Se", "covariance_Sp", "Se_serial_overall", "Sp_serial_overall", "logit_Sp_Luminex", "prev"), facet_args=list(ncol=2)) + theme_minimal() + ggtitle("Model allowing for conditional dependence in Se and Sp between ELISA and Luminex + SCALE-IT results")
