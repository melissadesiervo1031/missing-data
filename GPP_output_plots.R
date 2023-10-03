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




allmethodsGPP<-ggplot(data=paramallARIMASTAN3, aes(x=as.numeric(missingprop), y=value, color=type))+
  facet_grid(~factor(param, levels=c("intercept", "phi", "light", "discharge"),exclude = NA)~ missingness, scales="free_y")+
  geom_hline(data=trueestdf2, aes(yintercept=value), colour="gray")+
  #geom_hline(data=arimaestdf, aes(yintercept=value), colour="light blue")+
  geom_point(size=0.75, position = position_dodge(width=0.03))+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), size=0.3, width=0, position = position_dodge(width=0.03))+
  theme_bw()+
  xlab("Proportion of missing data")+ theme(legend.position="top")+theme(legend.title=element_blank())+
  ylab("Parameter estimate")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8))





####

