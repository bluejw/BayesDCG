############################################### Simulated Data: Scenario II #################################################

seed_index <- 1
set.seed(seed_index)

# load the treatment history data
n <- data$n # number of individuals
J <- data$J # number of visits

eta <- 1 # effect size
Q <- 3 # number of outcomes 
S <- 3 # number of covariates
L <- 2 # lags 

X <- array(NA, dim=c(n, max(J), S)) # continuous covariates 
for (i in 1:n){
  X[i,1:J[i],1] <- sample(c(0,1), J[i], prob=c(0.6,0.4), replace=T) # time-variant discrete covariate
  X[i,1:J[i],2] <- rnorm(J[i], 0, 1) # time-variant continuous covariate
  X[i,1:J[i],3] <- rnorm(J[i], 0, 1) # time-variant continuous covariate
}

mu <- c(1, -1, 0) # global intercept
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
beta1 <- rbind(c(0,0.5,0),c(0,0,0.25),c(0.1,0,0))*eta
beta2 <- rbind(c(0.5,0.25,0),c(0,0.5,0.125),c(0,0,0.5))*eta
alpha1 <- rbind(c(0.75,0,0),c(0,-0.5,0),c(0,0,0.25))*eta

# generate outcome data
data <- NULL; Nit <- 100
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
    Y_update[i,1,1] <- EY[i,1,1] + eta*0.5*Y[i,1,2] + eta*0.75*X[i,1,1] + mu[1]
    Y_update[i,1,2] <- EY[i,1,2] + eta*0.25*Y[i,1,3] - eta*0.5*X[i,1,2] + mu[2]
    Y_update[i,1,3] <- EY[i,1,3] + eta*0.1*Y[i,1,1] + eta*0.25*X[i,1,3] + mu[3]
    for (j in 2:J[i]){
      Y_update[i,j,1] <- EY[i,j,1] + eta*0.5*Y[i,j,2] + eta*0.25*Y[i,j-1,2] + eta*0.5*Y[i,j-1,1] + eta*0.75*X[i,j,1] + mu[1]
      Y_update[i,j,2] <- EY[i,j,2] + eta*0.25*Y[i,j,3] + eta*0.125*Y[i,j-1,3] + eta*0.5*Y[i,j-1,2] - eta*0.5*X[i,j,2] + mu[2]
      Y_update[i,j,3] <- EY[i,j,3] + eta*0.1*Y[i,j,1] + eta*0*Y[i,j-1,1] + eta*0.5*Y[i,j-1,3] + eta*0.25*X[i,j,3] + mu[3]
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
nu_0 <- 5e-5*eta
a_nu <- 5; b_nu <- 50
a_rho <- 0.5; b_rho <- 0.5 
a_sigma <- 1; b_sigma <- 1
