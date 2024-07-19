# File descriptions

* `gauss_Real_mods/Beartooth_scripts/modelruns_arima_MAR.sh`: Shell script for running modelruns_arima_MAR.R on HPC, very basic- I don't think it needs commenting

* `gauss_Real_mods/Beartooth_scripts/modelruns_brms_MAR.sh`: Shell script for running modelruns_brms_MAR.R on HPC analysis, very basic- I don't think it needs commenting

* `gauss_Real_mods/modelruns_arima_MAR.R`: R script prepared to run 3 ARIMA methods on the Gaussian real data with random missingness: drop missing, Kalman, multiple imputations, commenting is mostly clear, but some notes such as "data still doesn't work!" or "initially a lot of errors" should be cleaned up before publishing

* `gauss_Real_mods/modelruns_arima_MNAR.R `: R script prepared to run 3 ARIMA methods on the Gaussian real data with non-random missingness: drop missing, Kalman, multiple imputations, script is very similar to previous, comments look fine

* `gauss_Real_mods/modelruns_brms_MAR.R`: R script prepared to run brms methods on the Gaussian real data with random missingness, however the comments are talking about ARIMA methods, and were probably copied and pasted from some ARIMA codes... new commenting needed

* `gauss_Real_mods/modelruns_brms_MNAR.R`: R script prepared to run brms methods on the Gaussian real data with non-random missingness, one comment refers to "MODEL RUN ARIMA DROP" assuming this is not correct?

* `gauss_Sim_mods/Beartooth_scripts/modelruns_arima_MAR_A.sh`: Shell script for running modelruns_arima_MAR_A.R on HPC

* `gauss_Sim_mods/Beartooth_scripts/modelruns_arima_MAR_B.sh`: Shell script for running modelruns_arima_MAR_B.R on HPC

* `gauss_Sim_mods/Beartooth_scripts/modelruns_arima_MNAR.sh`: Unsure if this shell script is labeled right? It appears to run a different R script than indicated by it's name, a script that is now defunct? (modelruns_arima_MARautocor01.R)

* `gauss_Sim_mods/Beartooth_scripts/modelruns_brms_MAR_A.sh`: Shell script for running modelruns_brms_MAR_A.R on HPC

* `gauss_Sim_mods/Beartooth_scripts/modelruns_brms_MAR_B.sh`: Shell script for running modelruns_brms_MAR_B.R on HPC

* `gauss_Sim_mods/Beartooth_scripts/modelruns_brms_MNAR.sh`: Shell script for running modelruns_brms_MNAR.R on HPC

* `gauss_Sim_mods/modelruns_arima_MAR_A_Amelia_testing_05_14_24.R`: Appears similar to other R files that run ARIMA models, if it is for testing, I'm unsure that we still need it?

* `gauss_Sim_mods/modelruns_arima_MAR_A.R`: R script to run ARIMA models (drop missing, Kalman, multiple imputations) on the first half of the Gaussian simulated data with random missingness, seems straightforward

* `gauss_Sim_mods/modelruns_arima_MAR_B.R`: R script to run ARIMA models (drop missing, Kalman, multiple imputations) on the second half of the Gaussian simulated data with random missingness, almost exact same as for first half 

* `gauss_Sim_mods/modelruns_arima_MNAR.R`: R script to run ARIMA models (drop missing, Kalman, multiple imputations) on the Gaussian simulated data with non-random missingness, almost exact same as for first half of ARIMA MAR 

* `gauss_Sim_mods/modelruns_brms_MAR_A.R`: R script to run brms models on the first half of the Gaussian simulated data with random missingness, needs re-commenting to remove ARIMA labels

* `gauss_Sim_mods/modelruns_brms_MAR_B.R`: R script to run brms models on the second half of the Gaussian simulated data with random missingness, almost exact same as for first half, and like first half needs re-commenting to remove ARIMA labels

* `gauss_Sim_mods/modelruns_brms_MNAR.R`: R script to run brms models on the second half of the Gaussian simulated data with non-random missingness, similar to MAR first half, and also needs re-commenting to remove ARIMA labels
