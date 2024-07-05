# File descriptions

* `file 1`: An explanation of what this file does, and how it should be used in the analysis


* `file 2 `: An explanation of what this file does, and how it should be used in the analysis


## Josh's questions when going through code

* still going through `EM_fits2sims.R`
* `fit_Ricker_greattit.R` looks good
* `ricker3_.stan` looks good

### simulate_pop_Ricker.R

* Why are there two different Ricker models defined?
* How was the alpha value (0.0009997187) chosen for the first two plots? Do we actually use these plots
* for `missingProbs <- seq(0,1, by = .05)`, why is the max 1 when we only go up to 0.75? It is weird because after the for loop, the max missingness that is output is 0.75, presumably because of the error given by the nls function in the for loop. This is very strange to me

`Error in nls(popplus1 ~ Broods * exp(r - alpha * Broods), start = list(r = 0.1,  : 
  number of iterations exceeded maximum of 50`

