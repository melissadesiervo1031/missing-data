library(tidyverse)
library(brms)

# read in data
gauss_sim_randMiss_autoCorr_01 <- readRDS("data/missingDatasets/gauss_sim_randMiss_A.rds")

## subset the data 
# get simulated GPP data
sim_list <- gauss_sim_randMiss_autoCorr_01[[c(1)]]$y
# get simulated covariates (light and discharge)
sim_params <- gauss_sim_randMiss_autoCorr_01[[1]]$sim_params

# combine into a list to use in the brms function
# each element of the list contains the GPP values w/ missingness,  as well as light and discharge
# in this example, all of these datasets were made by using the same simulation data
simmissingdf <-lapply(X = sim_list, 
                      FUN = function(X) cbind.data.frame(GPP = X, 
                                                         light = sim_pars$X[,2], 
                                                         discharge = sim_pars$X[,3]))
# Make the model formula and priors
bform <- brms::bf(GPP | mi() ~ light + discharge + ar(p = 1))
bprior <- c(prior(normal(0,1), class = 'ar', ub = 1, lb = 0),
            prior(normal(0,5), class = 'b'))

# fit model to list of datasets
bmod <- brms::brm_multiple(bform, data = simmissingdf, 
                           prior = bprior, iter = 4000, 
                           combine = FALSE)

# a function to extract the brms derived parameters
extract_brms_pars <- function(bfit, include_missing = FALSE){
  bsum <- brms::posterior_summary(bfit, probs = c(0.025, 0.5, 0.975))
  bsum <- as.data.frame(bsum) %>%
    mutate(parameter = row.names(bsum)) %>%
    filter(!(parameter %in% c('lprior', 'lp__'))) %>%
    mutate(parameter = case_when(parameter == 'ar[1]' ~ 'phi',
                                 TRUE ~ parameter)) %>%
    select(parameter, mean = Estimate, sd = Est.Error, 
           '2.5%' = Q2.5, '50%' = Q50, '97.5%' = Q97.5)
  
  if(!include_missing){
    bsum <- bsum[grep('^Ymi', bsum$parameter, invert = TRUE),]
  }
  
  row.names(bsum) <- NULL
  
  return(bsum)
}
# actually extract the coefficients
bpars <- lapply(bmod, extract_brms_pars, include_missing = FALSE)
names(bpars) <- names(simmissingdf)
# return a list of the brms parameters, paired with the simulation parameters
outDat <- list(brms_pars = bpars,
            sim_params = sim_pars)


outDat_normPrior_bounds <- outDat

outDat_normPrior_bounds_phi <- map_df(c(1:length(outDat_normPrior_bounds$brms_pars)), function(x) 
    data.frame("name" = names(outDat_normPrior_bounds$brms_pars)[x],
               outDat_normPrior_bounds$brms_pars[[x]]
    )
    ) %>% 
  filter(parameter == "phi")

outDat_normPrior_noBounds_phi$prior <- "normPrior_noBounds"
outDat_uniformPrior_bounds_phi$prior <- "uniformPrior_bounds"
outDat_normPrior_bounds_phi$prior <- "normPrior_bounds"

outDat_all <- rbind(
  outDat_normPrior_noBounds_phi, outDat_uniformPrior_bounds_phi, outDat_normPrior_bounds_phi)

ggplot(data = outDat_all) +
  geom_density(aes(mean, col = prior)) + 
  geom_vline(aes(xintercept = 0.786986918))

## get the STAN code used in the mode
make_stancode(bform, data = simmissingdf[[1]], 
              prior = bprior, iter = 4000)
