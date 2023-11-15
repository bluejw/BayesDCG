# BayesDCG

**BayesDCG**: **Bayes**ian Structural Learning of **D**irected **C**yclic **G**raphs

R data and code for the paper:
"Directed Cyclic Graphs for Simultaneous Discovery of Time-Lagged and Instantaneous Causality from Time-Series Data".

## Instructions for Use

This repository contains the R code generating simulated datasets and implementing the Bayesian structural learning algorithm proposed in the paper "Directed Cyclic Graphs for Simultaneous Discovery of Time-Lagged and Instantaneous Causality from Time-Series Data".

In the folder "Data_and_Code":

* The R script "Data_Generate_Scenario_I.R" generates a simulated dataset (i.e., the simulation scenario I of the paper); The R data file "Treatment_History_Data.Rdata" contains the treatment history data for n=200 individuals randomly sampled from the Women's Interagency HIV Study (WIHS) dataset, which will be used to generate a simulated dataset (i.e., the simulation scenario II of the paper) in the R script "Data_Generate_Scenario_II.R";
