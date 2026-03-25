#/////////////////////
# A figure to compare different missingness types
# Alice Stears
# updated 19 July 2023
#/////////////////////

# load packages -----------------------------------------------------------
library(tidyverse)
library(ggpubr)

# function ----------------------------------------------------------------
# source the makeMissing() function
source("./Functions/missing_data_functions.R")

# Making Figure ------------------------------------------------

# use the simulated gaussian AR1 data
gausSim <- readRDS("./data/gauss_ar1_0miss_datasets.rds")
# for now, get just the first dataset
gausSim <- data.frame("y" = gausSim[[1]]$y,
                      "time" = 1:length(gausSim[[1]]$y))

# get 40% missing completely at random
gausSim$y_lowAuto_40Miss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .4, autoCorr = 0)))

# get 40% missing in highly autocorrelated chunks
gausSim$y_highAuto_40miss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .4, autoCorr = .9)))

# get 10% missing completely at random
gausSim$y_lowAuto_10Miss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .1, autoCorr = 0)))

# get 10% missing in highly autocorrelated chunks
gausSim$y_highAuto_10miss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .1, autoCorr = .9)))

# get 40% missing in minmax of data
gausSim$y_minMaxMiss_40 <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "minMax", propMiss = .4)))

# get 10% missing in minmax of data
gausSim$y_minMaxMiss_10 <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "minMax", propMiss = .1)))

gausSim$density_x <- seq(from = -5, to = 7, length.out = nrow(gausSim))
gausSim$density_y <- dnorm(gausSim$density_x, mean = mean(gausSim$y), sd = sd(gausSim$y))

## put data into a long data.frame
gausSim <- gausSim %>% 
  pivot_longer(cols = c(y, y_lowAuto_40Miss, y_lowAuto_10Miss, y_highAuto_10miss, y_highAuto_40miss, y_minMaxMiss_40, y_minMaxMiss_10),
               names_to = "MissingnessType")
# fix names
gausSim[gausSim$MissingnessType == "y", "MissingnessType"] <- "A: No Missing Data"
gausSim[gausSim$MissingnessType == "y_lowAuto_40Miss", "MissingnessType"] <- "D: MCAR: 40% missing, 0 autocorr."
gausSim[gausSim$MissingnessType == "y_lowAuto_10Miss", "MissingnessType"] <- "B: MCAR: 10% missing, 0 autocorr."
gausSim[gausSim$MissingnessType == "y_highAuto_10miss", "MissingnessType"] <- "C: MCAR: 10% missing, .9 autocorr."
gausSim[gausSim$MissingnessType == "y_highAuto_40miss", "MissingnessType"] <- "E: MCAR: 40% missing, .9 autocorr."
gausSim[gausSim$MissingnessType == "y_minMaxMiss_10", "MissingnessType"] <- "F: MNAR: 10% missing"
gausSim[gausSim$MissingnessType == "y_minMaxMiss_40", "MissingnessType"] <- "G: MNAR: 40% missing"

gausSim$MissingnessType <- factor(gausSim$MissingnessType, 
                                      levels = c("A: No Missing Data", "B: MCAR: 10% missing, 0 autocorr.", "C: MCAR: 10% missing, .9 autocorr.",
                                                 "D: MCAR: 40% missing, 0 autocorr.", "E: MCAR: 40% missing, .9 autocorr.", "F: MNAR: 10% missing", "G: MNAR: 40% missing"),
                                      ordered = TRUE)
# make histograms like Dusty made (looking at histogram of values in the missing
# data time series vs. expected distribution from the AR1 process)
# no missing data

ts_lines <- ggplot(data = gausSim, aes(x = time, y = value)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~MissingnessType, ncol = 1 ) + 
  ylab("y") +
  theme_classic() +
  theme(strip.text = element_text(hjust = 0))

ts_hist <- ggplot(data = gausSim) +
  geom_histogram(aes(value, after_stat(density)), fill = "grey", color = "darkgrey") +
  geom_line(aes(y = density_y, x = density_x), color = "blue") +
  facet_wrap(~MissingnessType, ncol = 1) + 
  theme_classic() +
  xlab("y") + 
  theme(strip.text = element_blank())


## save the figure to file
png(filename = "./figures/CompareMissingnessTypes_fig.png", width = 7.5, height = 7.5, units = "in", res = 700)
ggarrange(ts_lines,
          ts_hist,
          widths = c(.75,.25),
          ncol = 2,
          nrow = 1)
dev.off()
