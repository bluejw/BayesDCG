##########################################################################################
#                                   MCMC Main                                            #
##########################################################################################

library(Rcpp)
library(RcppArmadillo)
library(Matrix)
library(MCMCpack)
library(LaplacesDemon)

Rcpp::sourceCpp('MCMC_Rcpp_Functions.cpp')
source('MCMC_Functions.R')

# Generate Simulated Data: Scenario I
source('Data_Generate_Scenario_I.R') 

# Generate Simulated Data: Scenario II
load("Treatment_History_Data.Rdata")
source('Data_Generate_Scenario_II.R') 


# MCMC Setup
Nit <- 5000 # number of MCMC iterations
burn.in <- 2500 # burn-in period
thin.fac <- 5 # thinning factor 
post_index <- seq(burn.in+1, Nit, by=thin.fac) # index of posterior samples
post_num <- (Nit-burn.in)/thin.fac # number of posterior samples

mcmc <- NULL # list of MCMC samples
mcmc$beta <- array(NA, dim=c(Nit, L+1, Q, Q))
mcmc$gamma_beta <- array(NA, dim=c(Nit, L+1, Q, Q))
mcmc$nu_beta <- array(NA, dim=c(Nit, L+1, Q, Q))
mcmc$rho_beta <- array(NA, dim=c(Nit))
mcmc$alpha <- array(NA, dim=c(Nit, L+1, Q, S))
mcmc$gamma_alpha <- array(NA, dim=c(Nit, L+1, Q, S))
mcmc$nu_alpha <- array(NA, dim=c(Nit, L+1, Q, S))
mcmc$rho_alpha <- array(NA, dim=c(Nit))
mcmc$mu <- array(NA, dim=c(Nit, Q))
mcmc$sigma2 <- array(NA, dim=c(Nit, Q))


# Initialize
seed_index <- 1
set.seed(seed_index)

initial <- init()
mcmc$beta[1,,,] <- initial$beta
mcmc$gamma_beta[1,,,] <- initial$gamma_beta
mcmc$nu_beta[1,,,] <- initial$nu_beta
mcmc$rho_beta[1] <- initial$rho_beta
mcmc$alpha[1,,,] <- initial$alpha
mcmc$gamma_alpha[1,,,] <- initial$gamma_alpha
mcmc$nu_alpha[1,,,] <- initial$nu_alpha
mcmc$rho_alpha[1] <- initial$rho_alpha
mcmc$Z <- initial$Z
mcmc$mu[1,] <- initial$mu
mcmc$sigma2[1,] <- initial$sigma2


# Start of the Chain
start.time = proc.time()

for (nit in 2:Nit){
  
  print(nit)
  
  # Update beta
  mcmc$beta[nit,1,,] <- update_beta1(mcmc$beta[nit-1,1,,], mcmc$beta[nit-1,2,,], mcmc$beta[nit-1,3,,], 
                                     mcmc$gamma_beta[nit-1,,,], mcmc$nu_beta[nit-1,,,], 
                                     mcmc$alpha[nit-1,1,,], mcmc$alpha[nit-1,2,,], mcmc$alpha[nit-1,3,,], 
                                     mcmc$Z, mcmc$mu[nit-1,], mcmc$sigma2[nit-1,])
  mcmc$beta[nit,2,,] <- update_beta2_rcpp(mcmc$beta[nit,1,,], mcmc$beta[nit-1,3,,],
                                          mcmc$alpha[nit-1,1,,], mcmc$alpha[nit-1,2,,], mcmc$alpha[nit-1,3,,],
                                          mcmc$gamma_beta[nit-1,2,,], mcmc$nu_beta[nit-1,2,,],
                                          mcmc$Z, mcmc$mu[nit-1,], mcmc$sigma2[nit-1,], n, J, Q, S, Y, X)
  mcmc$beta[nit,3,,] <- update_beta3_rcpp(mcmc$beta[nit,1,,], mcmc$beta[nit,2,,],
                                          mcmc$alpha[nit-1,1,,], mcmc$alpha[nit-1,2,,], mcmc$alpha[nit-1,3,,],
                                          mcmc$gamma_beta[nit-1,3,,], mcmc$nu_beta[nit-1,3,,],
                                          mcmc$Z, mcmc$mu[nit-1,], mcmc$sigma2[nit-1,], n, J, Q, S, Y, X)
  
  # Update gamma_beta, nu_beta, rho_beta
  mcmc$gamma_beta[nit,,,] <- update_gamma_beta(mcmc$beta[nit,,,], mcmc$nu_beta[nit-1,,,], mcmc$rho_beta[nit-1])
  mcmc$nu_beta[nit,,,] <- update_nu_beta(mcmc$beta[nit,,,], mcmc$gamma_beta[nit,,,])
  mcmc$rho_beta[nit] <- update_rho_beta(mcmc$gamma_beta[nit,,,])
  
  # Update alpha
  mcmc$alpha[nit,1,,] <- update_alpha1_rcpp(mcmc$beta[nit,1,,], mcmc$beta[nit,2,,], mcmc$beta[nit,3,,],
                                            mcmc$alpha[nit-1,2,,], mcmc$alpha[nit-1,3,,],
                                            mcmc$gamma_alpha[nit-1,1,,], mcmc$nu_alpha[nit-1,1,,],
                                            mcmc$Z, mcmc$mu[nit-1,], mcmc$sigma2[nit-1,], n, J, Q, S, Y, X)
  mcmc$alpha[nit,2,,] <- update_alpha2_rcpp(mcmc$beta[nit,1,,], mcmc$beta[nit,2,,], mcmc$beta[nit,3,,],
                                            mcmc$alpha[nit,1,,], mcmc$alpha[nit-1,3,,],
                                            mcmc$gamma_alpha[nit-1,2,,], mcmc$nu_alpha[nit-1,2,,],
                                            mcmc$Z, mcmc$mu[nit-1,], mcmc$sigma2[nit-1,], n, J, Q, S, Y, X)
  mcmc$alpha[nit,3,,] <- update_alpha3_rcpp(mcmc$beta[nit,1,,], mcmc$beta[nit,2,,], mcmc$beta[nit,3,,],
                                            mcmc$alpha[nit,1,,], mcmc$alpha[nit,2,,],
                                            mcmc$gamma_alpha[nit-1,3,,], mcmc$nu_alpha[nit-1,3,,],
                                            mcmc$Z, mcmc$mu[nit-1,], mcmc$sigma2[nit-1,], n, J, Q, S, Y, X)
  
  # Update gamma_alpha, nu_alpha, rho_alpha
  mcmc$gamma_alpha[nit,,,] <- update_gamma_alpha(mcmc$alpha[nit,,,], mcmc$nu_alpha[nit-1,,,], mcmc$rho_alpha[nit-1])
  mcmc$nu_alpha[nit,,,] <- update_nu_alpha(mcmc$alpha[nit,,,], mcmc$gamma_alpha[nit,,,])
  mcmc$rho_alpha[nit] <- update_rho_alpha(mcmc$gamma_alpha[nit,,,])
  
  # Update Z 
  mcmc$Z <- update_Z_rcpp(mcmc$beta[nit,1,,], mcmc$beta[nit,2,,], mcmc$beta[nit,3,,], 
                          mcmc$alpha[nit,1,,], mcmc$alpha[nit,2,,], mcmc$alpha[nit,3,,], 
                          mcmc$mu[nit-1,], mcmc$sigma2[nit-1,], n, J, Q, S, Y, X)
  
  # Update mu 
  mcmc$mu[nit,] <- update_mu_rcpp(mcmc$beta[nit,1,,], mcmc$beta[nit,2,,], mcmc$beta[nit,3,,],
                                  mcmc$alpha[nit,1,,], mcmc$alpha[nit,2,,], mcmc$alpha[nit,3,,],
                                  mcmc$Z, mcmc$sigma2[nit-1,], n, J, Q, S, Y, X)
  
  # Update sigma2
  mcmc$sigma2[nit,] <- update_sigma2_rcpp(mcmc$beta[nit,1,,], mcmc$beta[nit,2,,], mcmc$beta[nit,3,,],
                                          mcmc$alpha[nit,1,,], mcmc$alpha[nit,2,,], mcmc$alpha[nit,3,,],
                                          mcmc$Z, mcmc$mu[nit,], n, J, Q, S, Y, X, a_sigma, b_sigma)
}

duration = proc.time()-start.time
# End of the Chain


# Posterior Inference
post <- NULL # posterior samples
post$beta <- mcmc$beta[post_index,,,]
post$gamma_beta <- mcmc$gamma_beta[post_index,,,]
post$alpha <- mcmc$alpha[post_index,,,]
post$gamma_alpha <- mcmc$gamma_alpha[post_index,,,]
post$mu <- mcmc$mu[post_index,]
post$sigma2 <- mcmc$sigma2[post_index,]

graph <- NULL # causal graph 
graph$beta1 <- (colMeans(post$gamma_beta[,1,,]) > 0.5) * 1.0
graph$beta2 <- (colMeans(post$gamma_beta[,2,,]) > 0.5) * 1.0
graph$beta3 <- (colMeans(post$gamma_beta[,3,,]) > 0.5) * 1.0
graph$alpha1 <- (colMeans(post$gamma_alpha[,1,,]) > 0.5) * 1.0
graph$alpha2 <- (colMeans(post$gamma_alpha[,2,,]) > 0.5) * 1.0
graph$alpha3 <- (colMeans(post$gamma_alpha[,3,,]) > 0.5) * 1.0

coef <- NULL # causal effect 
coef$beta1 <- colMeans(post$beta[,1,,]) * graph$beta1
coef$beta2 <- colMeans(post$beta[,2,,]) * graph$beta2
coef$beta3 <- colMeans(post$beta[,3,,]) * graph$beta3
coef$alpha1 <- colMeans(post$alpha[,1,,]) * graph$alpha1
coef$alpha2 <- colMeans(post$alpha[,2,,]) * graph$alpha2
coef$alpha3 <- colMeans(post$alpha[,3,,]) * graph$alpha3
