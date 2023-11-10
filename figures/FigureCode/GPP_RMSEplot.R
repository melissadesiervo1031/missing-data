# Load packages
library(here)
library(tidyverse)
library(rstan)
library(Amelia)
library(forecast)
library(xts)
library(nlme)
library(tidyverse)
library(lubridate)


## read in the cleaned beartooth ARIMA and STAN outputs together##

allbeartooth_gaussresults <- read_csv("data/model_results/allbeartooth_gaussresults.csv", show_col_types = FALSE)

##need to round##

allbeartooth_gaussresults2<-allbeartooth_gaussresults %>% mutate(missingprop_round=round(missingprop, 1), autocorr_round=round(autocorr, 1))

##summarize across the simulations ##

allsummary<-allbeartooth_gaussresults2 %>% group_by(parameter, missingprop_round, autocorr_round, type, missingness) %>% dplyr::summarise(sumdiffsquared=sum(diffsquared), n=n())%>% mutate(MSE=sumdiffsquared/(n-1))%>% mutate(RMSE=sqrt(MSE))


#Plot ###


allmethodsGPP<-ggplot(data=allsummary, aes(x=as.numeric(missingprop_round), y=RMSE, color=type))+
  facet_grid(~factor(parameter, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ missingness, scales="free_y")+
  geom_point(size=0.75, position = position_dodge(width=0.03))+
  theme_bw()+
  xlab("Proportion of missing data")+ theme(legend.position="top")+theme(legend.title=element_blank())+
  ylab("RMSE")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))

pdf(file = "Figures/GPPsimulationsfirstpass.pdf",width = 8, height = 5)
allmethodsGPP
dev.off()
