########################### 
# Code generating figure 3 
###########################

library(plotly)
library(tidyverse)
library(htmltools)


# Set working directory
setwd("INSERT WORKING DIRECTORY PATH")

# read in data
samples_collected<-read_csv("~/Supplementary Table 2.csv")

# redefine high and low 95% credible intervals as the difference from median for    plotting purposes
samples_collected$hici2<- samples_collected$`97.5% Credible Interval`-   samples_collected$`Adjusted seroprevalence`
samples_collected$lowci2<-samples_collected$`Adjusted seroprevalence`- samples_collected$`2.5% Credible Interval`


insurance_sero<-  samples_collected %>% filter(Class == "Insurance Type")
ethnicity_sero<- samples_collected %>% filter(Class == "Race/Ethnicity")
age_sero<-samples_collected %>% filter(Class == "Age Group")
sex_sero<-samples_collected %>% filter(Class == "Sex")

plot_figs<-function(dat){
  
  # plot parameters
  f <- list(
    size = 22
  )
  
  
  p<-plot_ly(
    data = dat,
    x = ~Group,
    y = ~`Adjusted seroprevalence`,
    text = ~Count,
    color= ~Group,
    error_y = list(
      type = "data",
      symmetric = FALSE,
      color = '#000000',
      arrayminus = ~lowci2,
      array = ~hici2)
  )%>%
    layout(xaxis=list( title = dat$Class[1], tickfont = list(size = 18),  titlefont = f),
           yaxis = list( title = "Adjusted seroprevalence", tickfont = list(size = 18),  titlefont = f), legend= list(font=list(size=18)))
  
  p
  return(p)
}

dat_list<-list(insurance_sero,age_sero,sex_sero,ethnicity_sero)
lapply(dat_list, plot_figs)