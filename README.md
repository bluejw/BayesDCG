# BayesDCG: Bayesian Structural Learning of Directed Cyclic Graphs for Time-Series Causal Discovery

R data and code for the paper:
"Directed Cyclic Graphs for Simultaneous Discovery of Time-Lagged and Instantaneous Causality from Time-Series Data".

## Instructions for Use

This repository contains the R code generating simulated datasets and implementing the Bayesian structural learning algorithm proposed in the paper "Directed Cyclic Graphs for Simultaneous Discovery of Time-Lagged and Instantaneous Causality from Time-Series Data".

In the folder "Data_and_Code":

* The R script "Data_Generate_Scenario_I.R" generates a simulated dataset (simulation scenario I); The R data file "Treatment_History_Data.Rdata" contains the treatment history data for n=200 individuals randomly sampled from the Women's Interagency HIV Study (WIHS) dataset, which will be used to generate a simulated dataset (simulation scenario II) in the R script "Data_Generate_Scenario_II.R";

* The R script "MCMC_R_Functions.R" provides R functions used for MCMC, and the Rcpp script "MCMC_Rcpp_Functions.cpp" provides Rcpp functions used for MCMC;

* The R script "MCMC_Main.R" runs the MCMC algorithm for the proposed Bayesian structural learning. 
