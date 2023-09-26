
#### 8/ 23/ 23 ###
## Melissa DeSiervo ###

### Clean beartooth output #####

library(tidyr)
library(dplyr)
library(stringr)

##  Read in the real simulation parameters ###

gauss_sim_params<-readRDS("data/missingDatasets/forBeartooth/gauss_sim_params.rds")

gauss_sim_param_long<-gather(gauss_sim_params, parameter, realestimate, phi:beta3, factor_key=TRUE)

gauss_sim_param_long$parameter <- str_replace(gauss_sim_param_long$parameter, "beta1", "intercept")
gauss_sim_param_long$parameter <- str_replace(gauss_sim_param_long$parameter, "beta2", "light")
gauss_sim_param_long$parameter <- str_replace(gauss_sim_param_long$parameter, "beta3", "discharge")


### ARIMA OUTPUT ####

Arima_MAR_A <- readRDS("data/BeartoothOutputData/gaussSim_MAR_A_arimaOut.rds")

### STAN BRMS output###

STAN_MAR_A<- read_csv("data/model_results/00_combined_gauss_sim_randMiss_A.csv", show_col_types = FALSE) #379999 X  10#
STAN_MAR_B<- read_csv("data/model_results/00_combined_gauss_sim_randMiss_B.csv", show_col_types = FALSE) # 379999 X 10#
STAN_MNAR<- read_csv("data/model_results/00_combined_gauss_sim_minMaxMiss.csv", show_col_types = FALSE) # 80999 X 10 # 

        ### remove the parameter rows ###
        
        STAN_MAR_A<-subset(STAN_MAR_A, parameter!="parameter") #375000 X  10#
        STAN_MAR_B<-subset(STAN_MAR_B, parameter!="parameter") #375000 X  10#
        STAN_MNAR<-subset(STAN_MNAR, parameter!="parameter") # 80000 X   10#
        
        ### change the names of params to match the ARIMA output ##
        
        STAN_MAR_A$parameter <- str_replace(STAN_MAR_A$parameter, "b_Intercept", "intercept")
        STAN_MAR_A$parameter <- str_replace(STAN_MAR_A$parameter, "b_light", "light")
        STAN_MAR_A$parameter <- str_replace(STAN_MAR_A$parameter, "b_discharge", "discharge")
        
        
        STAN_MAR_B$parameter <- str_replace(STAN_MAR_B$parameter, "b_Intercept", "intercept")
        STAN_MAR_B$parameter <- str_replace(STAN_MAR_B$parameter, "b_light", "light")
        STAN_MAR_B$parameter <- str_replace(STAN_MAR_B$parameter, "b_discharge", "discharge")
        
        STAN_MNAR$parameter <- str_replace(STAN_MNAR$parameter, "b_Intercept", "intercept")
        STAN_MNAR$parameter <- str_replace(STAN_MNAR$parameter, "b_light", "light")
        STAN_MNAR$parameter <- str_replace(STAN_MNAR$parameter, "b_discharge", "discharge")
        
        ### MNAR incorrectly labeled, fix ###
        
        STAN_MNAR$missingness <- str_replace(STAN_MNAR$missingness, "MAR", "MNAR")
        
        ### rbind the 3 STAN outputs ###
        
        STAN_MAR_A<-STAN_MAR_A %>% rename("2.5%" = "2.50%","97.5%" = "97.50%" )## first rename thess columns ##
        
        STAN_all<-rbind(STAN_MAR_A, STAN_MAR_B, STAN_MNAR)


      ### ASSIGN THE CORRECT PARAMETERS FROM THE SIMULATIONS ###

      STAN_all$run_no<-as.numeric(STAN_all$run_no)

      unique(STAN_all$run_no)

      ## Run number needs to only be 1-1000# 

      STAN_all_2<-STAN_all %>% mutate(SimNumber=ifelse(run_no > 1000, run_no-1000, run_no),SimNumber=ifelse(SimNumber> 1000, SimNumber-1000, SimNumber),SimNumber=ifelse(SimNumber > 1000, SimNumber-1000, SimNumber),SimNumber=ifelse(SimNumber > 1000, SimNumber-1000, SimNumber))

      STAN_all_withest<-merge(STAN_all_2, gauss_sim_param_long, by=c("SimNumber", "parameter"))
      
      STAN_all_withest1<-STAN_all_withest %>% select(run_no, SimNumber,missingness, type, missingprop_autocor, parameter, value=mean, error=sd, realestimate)
      
      STAN_all_withest1$value<-as.numeric(STAN_all_withest1$value)
      STAN_all_withest1$error<-as.numeric(STAN_all_withest1$error)
      

### ARIMA OUTPUT ####
      
      Arima_MAR_A <- readRDS("data/BeartoothOutputData/gaussSim_MAR_A_arimaOut.rds")
      
      Arima_MAR_A1<-Arima_MAR_A %>% select (run_no=CurSim, SimNumber=simName, missingness, type, parameter=param, missingprop_autocor, value, error=SE) 
      
      Arima_MNAR <- readRDS("data/BeartoothOutputData/gaussSim_MNAR_arimaOut.rds")
      
      Arima_MNAR1<-Arima_MNAR %>% select (run_no=CurSim, SimNumber=CurSim, missingness, type, parameter=param, missingprop_autocor, value, error=SE) 
      
      Arima_all<-rbind(Arima_MAR_A1, Arima_MNAR1)
      
      Arima_all$parameter <- str_replace(Arima_all$parameter, "ar1", "phi")
      
      
      ### NEED TO GET THE MNAR ONES TOO...BUT JUST THIS FOR NOW IS OKAY ###
      
      # real param #      
      Arima_realest<-Arima_MAR_A %>% select(SimNumber=simName, phi_sim, beta1_sim, beta2_sim, beta3_sim) 
      
      Arima_realest_long<-gather(Arima_realest, parameter, realestimate, phi_sim:beta3_sim, factor_key=TRUE)
            
      Arima_realest_long$parameter <- str_replace(Arima_realest_long$parameter, "phi_sim", "phi")
      Arima_realest_long$parameter <- str_replace(Arima_realest_long$parameter, "beta1_sim", "intercept")
      Arima_realest_long$parameter <- str_replace(Arima_realest_long$parameter, "beta2_sim", "light")
      Arima_realest_long$parameter <- str_replace(Arima_realest_long$parameter, "beta3_sim", "discharge")
      
      Arima_realest_long<-distinct(Arima_realest_long)
      
      ## merge the output with the real est ##
      
      ARIMA_all_withest<-merge(Arima_all, Arima_realest_long, by=c("SimNumber", "parameter"))

      ARIMA_all_withest2<-ARIMA_all_withest %>% select(run_no, SimNumber, missingness, type, missingprop_autocor, parameter, value, error, realestimate)
      
      
      
      
### MERGE THE STAN AND ARIMA OUTPUT!!! #####
      
Arimastan_output_all<-rbind(ARIMA_all_withest2, STAN_all_withest1)


      
### split the missingprop_autocor ##

Arimastan_output_all2<-separate(Arimastan_output_all, missingprop_autocor, into = c("missingprop1", "autocorr1"), sep = 16, remove = FALSE) 

Arimastan_output_all2$missingprop<-as.numeric(str_extract(Arimastan_output_all2$missingprop1, "\\d+\\.*\\d*"))
Arimastan_output_all2$autocorr<-as.numeric(str_extract(Arimastan_output_all2$autocorr1, "\\d+\\.*\\d*"))

Arimastan_output_all3<-Arimastan_output_all2 %>% mutate(diff=realestimate-value, diffsquared=diff^2) %>% select(run_no, SimNumber, missingness, type, missingprop, autocorr,parameter, value, error, realestimate, diff, diffsquared)

## write a csv for making the figure ###

write.csv(Arimastan_output_all3, "data/model_results/allbeartooth_gaussresults.csv")
