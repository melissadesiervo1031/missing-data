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
gausSim$y_lowAuto_40Miss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .4, autoCorr = .01)))

# get 40% missing in highly autocorrelated chunks
gausSim$y_highAuto_40miss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .4, autoCorr = .9)))

# get 10% missing completely at random
gausSim$y_lowAuto_10Miss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .1, autoCorr = .01)))

# get 10% missing in highly autocorrelated chunks
gausSim$y_highAuto_10miss <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "random", propMiss = .1, autoCorr = .9)))

# get 40% missing in minmax of data
gausSim$y_minMaxMiss_40 <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "minMax", propMiss = .4)))

# get 10% missing in minmax of data
gausSim$y_minMaxMiss_10 <- as.vector(unlist(makeMissing(timeSeries = gausSim$y, typeMissing = "minMax", propMiss = .1)))

# make histograms like Dusty made (looking at histogram of values in the missing
# data time series vs. expected distribution from the AR1 process)
# no missing data
gauss_noMiss_line <- ggplot(data = gausSim, aes(x = time, y = y)) +
  geom_line() +
  geom_point(size = 1) +
  ggtitle("No Missing Data") +
  ylab("y") + 
  theme_classic()

gauss_noMiss_hist <- ggplot() +
  geom_histogram(data = gausSim, aes(y, after_stat(density)), fill = "grey", color = "darkgrey") +
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)),
                x = seq(-5,7,.1)), color = "blue") +
  xlab("y") +
  theme_classic()

# randomly missing (low autocorrellation, 10% missing)
gauss_lowAuto_10Miss_line <- ggplot(data = gausSim, aes(x = time, y = y_lowAuto_10Miss)) +
  geom_line() +
  geom_point(size = 1) +
  ggtitle("Missing at Random: 10% missing, autoCorr = .01") +
  ylab("y") +
  theme_classic()

gauss_lowAuto_10Miss_hist <- ggplot() +
  geom_histogram(data = gausSim, aes(y_lowAuto_10Miss, after_stat(density)), fill = "grey", color = "darkgrey") +
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)),
                x = seq(-5,7,.1)), color = "blue") +
  xlab("y") +
  theme_classic()

# randomly missing (low autocorrellation, 40% missing)
gauss_lowAuto_40Miss_line <- ggplot(data = gausSim, aes(x = time, y = y_lowAuto_40Miss)) +
  geom_line() +
  geom_point(size = 1) +
  ggtitle("Missing at Random: 40% missing, autoCorr = .01") +
  ylab("y") +
  theme_classic()

gauss_lowAuto_40Miss_hist <- ggplot() +
  geom_histogram(data = gausSim, aes(y_lowAuto_40Miss, after_stat(density)), fill = "grey", color = "darkgrey") +
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)),
                x = seq(-5,7,.1)), color = "blue") +
  xlab("y") +
  theme_classic()

# missing randomly spaced chunks (high autocorr, 10% missing)
gauss_highAuto_10Miss_line <- ggplot(data = gausSim, aes(x = time, y = y_highAuto_10miss)) +
  geom_line() +
  geom_point(size = 1) +
  ggtitle("Missing at Random: 10% missing, autoCorr = .9") +
  ylab("y") +
  theme_classic()

gauss_highAuto_10Miss_hist <- ggplot() +
  geom_histogram(data = gausSim, aes(y_lowAuto_10Miss, after_stat(density)), fill = "grey", color = "darkgrey") +
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)),
                x = seq(-5,7,.1)), color = "blue") +
  xlab("y") +
  theme_classic()

# missing randomly spaced chunks (high autocorr, 40% missing)
gauss_highAuto_40Miss_line <- ggplot(data = gausSim, aes(x = time, y = y_highAuto_40miss)) +
  geom_line() +
  geom_point(size = 1) +
  ggtitle("Missing at Random: 40% missing, autoCorr = .9") +
  ylab("y") +
  theme_classic()

gauss_highAuto_40Miss_hist <- ggplot() +
  geom_histogram(data = gausSim, aes(y_lowAuto_40Miss, after_stat(density)), fill = "grey", color = "darkgrey") +
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)),
                x = seq(-5,7,.1)), color = "blue") +
  xlab("y") +
  theme_classic()

# missing at min and max (10%)
gauss_minMaxMiss_10miss_line <- ggplot(data = gausSim, aes(x = time, y = y_minMaxMiss_10)) +
  geom_line() +
  geom_point(size = 1) +
  ggtitle("Missing Not at Random: 10% missing") +
  ylab("y") +
  ylim(c(-2.5,5)) +
  theme_classic()

gauss_minMaxMiss_10miss_hist <- ggplot() +
  geom_histogram(data = gausSim, aes(y_minMaxMiss_10, after_stat(density)), fill = "grey", color = "darkgrey") +
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)),
                x = seq(-5,7,.1)), color = "blue") +
  xlab("y") +
  theme_classic()

# missing at min and max (40%)
gauss_minMaxMiss_40miss_line <- ggplot(data = gausSim, aes(x = time, y = y_minMaxMiss_40)) +
  geom_line() +
  geom_point(size = 1) +
  ggtitle("Missing Not at Random: 40% missing") +
  ylab("y") +
  ylim(c(-2.5,5)) +
  theme_classic()

gauss_minMaxMiss_40miss_hist <- ggplot() +
  geom_histogram(data = gausSim, aes(y_minMaxMiss_40, after_stat(density)), fill = "grey", color = "darkgrey") +
  geom_line(aes(y = dnorm(seq(-5,7,.1), mean = mean(gausSim$y), sd = sd(gausSim$y)),
                x = seq(-5,7,.1)), color = "blue") +
  xlab("y") +
  theme_classic()

## save the figure to file
png(filename = "./figures/CompareMissingnessTypes_fig.png", width = 700, height = 600)
ggarrange(gauss_noMiss_line, gauss_noMiss_hist,
          gauss_lowAuto_10Miss_line, gauss_lowAuto_10Miss_hist,
          gauss_lowAuto_40Miss_line, gauss_lowAuto_40Miss_hist,
          gauss_highAuto_10Miss_line, gauss_highAuto_10Miss_hist,
          gauss_highAuto_40Miss_line, gauss_highAuto_40Miss_hist,
          gauss_minMaxMiss_10miss_line, gauss_minMaxMiss_10miss_hist,
          gauss_minMaxMiss_40miss_line, gauss_minMaxMiss_40miss_hist,
          widths = c(.75,.25, .75,.25, .75,.25, .75,.25, .75,.25,.75,.25,.75,.25),
          ncol = 2,
          nrow = 7)
dev.off()
