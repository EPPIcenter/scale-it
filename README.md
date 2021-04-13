# Serosurveillance for Continuous, ActionabLe Epidemiologic Intelligence of Transmission (SCALE-IT)

Code to reproduce results described in "Citywide serosurveillance of the initial SARS-CoV-2 outbreak in San Francisco" by Routledge, Epstein, Takahashi et al. https://www.researchsquare.com/article/rs-180966/v1.

https://github.com/EPPIcenter/scale-it/blob/main/Function_sampling_algorithm.R contains the function used to select samples for analysis, stratifying by age and zipcode. 
   
https://github.com/EPPIcenter/scale-it/blob/main/Generate_fig_3.R generates Figure 3 in the manuscript, using Supplementary Table 2.
 
https://github.com/EPPIcenter/scale-it/blob/main/Generate_seroprev_category.R generates posterior seroprevalence estimates, adjusted by test performance characteristics, by demographic group and neighborhood , using Supplementary Table 2 and Supplementary Table 3.
 
https://github.com/EPPIcenter/scale-it/blob/main/Generate_seroprevalence_weighted.R generates posterior seroprevalence estimates overall, weighted by age and sex and adjusted by test performance characteristics.
