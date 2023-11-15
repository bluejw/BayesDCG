############################################### Simulated Data: Scenario I #################################################

seed_index <- 1
set.seed(seed_index)

n <- 1000 # number of subjects
J <- rep(5, n) # number of visits
Q <- 4 # number of outcomes 
S <- 3 # number of covariates
L <- 2 # number of lags

X <- array(NA, dim=c(n, max(J), S)) # continuous covariates 
for (i in 1:n){
  for (s in 1:S){
    X[i,1:J[i],s] <- rnorm(J[i], 0, 1) 
  }
}

mu <- rep(0, Q) # global intercept
sigma2 <- rep(0.125, Q) # variance of i.i.d non-Gaussain errors
EY <- array(NA, dim=c(n, max(J), Q)) # outcome errors
for (i in 1:n){
  for (j in 1:J[i]){
    for (q in 1:Q){
      EY[i,j,q] <- rlaplace(1, location=0, scale=2*sqrt(sigma2[q]))
    }
  }
}

# simulation truth
beta1 <- rbind(c(0,0,-0.95,0),c(1.05,0,0,0),c(0,1,0,1),c(-0.1,0,0,0))
alpha1 <- rbind(c(0.5,0,0),c(0,0,0),c(0,0,0),c(0,0,0))

# generate data
data <- NULL; Nit <- 1000
data$Y <- array(NA, dim=c(Nit, n, max(J), Q))
for (i in 1:n){
  for (j in 1:J[i]){
    # initial values
    data$Y[1,i,j,] <- rmvn_rcpp(1, rep(0,Q), diag(1,Q))
  }
}

generate_data <- function(Y){
  Y_update <- array(NA, dim=c(n, max(J), Q))
  for (i in 1:n){
    for (j in 1:J[i]){
      Y_update[i,j,1] <- EY[i,j,1] - 0.95*Y[i,j,3] + 0.5*X[i,j,1] + mu[1]
      Y_update[i,j,2] <- EY[i,j,2] + 1.05*Y[i,j,1] + mu[2]
      Y_update[i,j,3] <- EY[i,j,3] + 1*Y[i,j,2] + 1*Y[i,j,4] + mu[3]
      Y_update[i,j,4] <- EY[i,j,4] - 0.1*Y[i,j,1] + mu[4]
    }
  }
  return(Y_update)
}

for (nit in 2:Nit){
  # iterative process
  data$Y[nit,,,] <- generate_data(data$Y[nit-1,,,])
}
Y <- data$Y[Nit,,,]

# hyper-parameters
nu_0 <- 2.5e-4 
a_nu <- 5; b_nu <- 50
a_rho <- 0.5; b_rho <- 0.5 
a_sigma <- 1; b_sigma <- 1
