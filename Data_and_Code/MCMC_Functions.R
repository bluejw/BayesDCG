################################################## Simulation Functions ###########################################



init <- function(){
  
  # MCMC initialization 
  # initialize parameters need to be estimated in MCMC
  
  alpha_init <- array(NA, dim=c(L+1, Q, S))
  gamma_alpha_init <- array(NA, dim=c(L+1, Q, S))
  nu_alpha_init <- array(NA, dim=c(L+1, Q, S))
  rho_alpha_init <- 0.5
  
  beta_init <- array(NA, dim=c(L+1, Q, Q))
  gamma_beta_init <- array(NA, dim=c(L+1, Q, Q))
  nu_beta_init <- array(NA, dim=c(L+1, Q, Q))
  rho_beta_init <- 0.5
  
  Z_init <- array(NA, dim=c(n, max(J), Q))
  mu_init <- rep(0, Q)
  sigma2_init <- rep(1, Q)
  
  # initialize alpha
  for (l in 1:(L+1)){
    for (q in 1:Q){
      for (s in 1:S){
        gamma_alpha_init[l,q,s] <- sample(c(nu_0,1), 1, prob=c(0.5,0.5)); nu_alpha_init[l,q,s] <- 0.01
        alpha_init[l,q,s] <- rnorm(1, 0, sqrt(gamma_alpha_init[l,q,s]*nu_alpha_init[l,q,s]))
      }
    }
  }
  
  # initialize beta
  for (l in 1:(L+1)){
    for (q in 1:Q){
      for (p in 1:Q){
        gamma_beta_init[l,q,p] <- sample(c(nu_0,1), 1, prob=c(0.5,0.5)); nu_beta_init[l,q,p] <- 0.01
        beta_init[l,q,p] <- rnorm(1, 0, sqrt(gamma_beta_init[l,q,p]*nu_beta_init[l,q,p]))
      }
    }
  }
  
  # initialize Z
  for (i in 1:n){
    for (j in 1:J[i]){
      for (q in 1:Q){
        Y_tilde_ijq <- Y[i,j,q] - t(Y[i,j,])%*%beta_init[1,q,] - t(X[i,j,])%*%alpha_init[1,q,]
        if (j > 1) { Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-1,])%*%beta_init[2,q,] - t(X[i,j-1,])%*%alpha_init[2,q,] }
        if (j > 2) { Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-2,])%*%beta_init[3,q,] - t(X[i,j-2,])%*%alpha_init[3,q,] }
        Z_init[i,j,q] <- rinvgaussian(1, mu=sqrt(sigma2_init[q])/(2*abs(Y_tilde_ijq)), lambda=1/4)
      }
    }
  }
  
  init_list <- list(alpha=alpha_init, gamma_alpha=gamma_alpha_init, nu_alpha=nu_alpha_init, rho_alpha=rho_alpha_init,
                    beta=beta_init, gamma_beta=gamma_beta_init, nu_beta=nu_beta_init, rho_beta=rho_beta_init,
                    Z=Z_init, mu=mu_init, sigma2=sigma2_init)
  return(init_list)            
}



logll_beta1 <- function(q, beta1, beta2, beta3, alpha1, alpha2, alpha3, Z, mu, sigma2){
  
  logll <- 0
  jacob_factor <- abs(det(diag(1, Q) - beta1))
  
  for (i in 1:n){
    # extra factor: absolute determinant of jacobian matrix 
    logll <- logll + J[i]*log(jacob_factor)
    for (j in 1:J[i]){
      Y_tilde_ijq <- Y[i,j,q] - mu[q] - t(Y[i,j,])%*%beta1[q,] - t(X[i,j,])%*%alpha1[q,]
      if (j > 1){ Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-1,])%*%beta2[q,] - t(X[i,j-1,])%*%alpha2[q,] }
      if (j > 2){ Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-2,])%*%beta3[q,] - t(X[i,j-2,])%*%alpha3[q,] }
      logll <- logll + dnorm(Y_tilde_ijq, 0, sqrt(sigma2[q]/Z[i,j,q]), log=T)
    }
  }
  
  return(logll)
}



update_beta1 <- function(beta1, beta2, beta3, gamma_beta, nu_beta, alpha1, alpha2, alpha3, Z, mu, sigma2){
  
  beta1_update <- beta1
  step <- 0.02 # step for normal random walk
  
  for (q in 1:Q){
    for (p in 1:Q){
      if (p == q){
        beta1_update[q,p] <- 0
      }else{
        # propose the new beta_1qp
        beta1_new <- beta1_update
        beta1_new[q,p] <- rnorm(1, beta1_update[q,p], step)
        
        # jacobian matrix norm constaint
        if (max(Mod(eigen(beta1_new)$values)) < 1){
          # calculate the acceptance ratio 
          ratio <- logll_beta1_rcpp(q-1, beta1_new, beta2, beta3, alpha1, alpha2, alpha3, Z, mu, sigma2, n, J, Q, S, Y, X) - 
            logll_beta1_rcpp(q-1, beta1_update, beta2, beta3, alpha1, alpha2, alpha3, Z, mu, sigma2, n, J, Q, S, Y, X) + 
            dnorm(x=beta1_new[q,p], 0, sqrt(gamma_beta[1,q,p]*nu_beta[1,q,p]), log=T) -
            dnorm(x=beta1_update[q,p], 0, sqrt(gamma_beta[1,q,p]*nu_beta[1,q,p]), log=T)
          # accept or reject
          test <- log(runif(1))
          if (test<ratio){ beta1_update <- beta1_new }
        }
      }
    }
  }
  
  return(beta1_update)
}



update_beta2 <- function(beta1, beta3, gamma_beta, nu_beta, alpha, Z, mu, sigma2){
  
  beta2_update <- matrix(NA, nrow=Q, ncol=Q)
  
  for (q in 1:Q){
    # calculate Y_sum and Yy_sum
    Y_sum <- matrix(0, nrow=Q, ncol=Q); Yy_sum <- rep(0, Q)
    for (i in 1:n){
      for (j in 2:J[i]){
        Y_sum <- Y_sum + Y[i,j-1,]%*%t(Y[i,j-1,])*Z[i,j,q]
        Y_tilde_ijq <- Y[i,j,q] - mu[q] - t(Y[i,j,])%*%beta1[q,] - t(X[i,j,])%*%alpha[1,q,] - t(X[i,j-1,])%*%alpha[2,q,]
        if (j > 2){ Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-2,])%*%beta3[q,] - t(X[i,j-2,])%*%alpha[3,q,] }
        Yy_sum <- Yy_sum + Y[i,j-1,]*c(Y_tilde_ijq)*Z[i,j,q]
      }
    }
    # Gaussian posterior distribution
    V_n <- chol2inv(chol(Y_sum/sigma2[q] + diag(1/(gamma_beta[2,q,]*nu_beta[2,q,]))))
    mu_n <- V_n %*% (Yy_sum/sigma2[q])
    beta2_update[q,] <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  }
  
  return(beta2_update)
}



update_beta3 <- function(beta1, beta2, gamma_beta, nu_beta, alpha, Z, mu, sigma2){
  
  beta3_update <- matrix(NA, nrow=Q, ncol=Q)
  
  for (q in 1:Q){
    # calculate Y_sum and Yy_sum
    Y_sum <- matrix(0, nrow=Q, ncol=Q); Yy_sum <- rep(0, Q)
    for (i in 1:n){
      for (j in 3:J[i]){
        Y_sum <- Y_sum + Y[i,j-2,]%*%t(Y[i,j-2,])*Z[i,j,q]
        Y_tilde_ijq <- Y[i,j,q] - mu[q] - t(Y[i,j,])%*%beta1[q,] - t(X[i,j,])%*%alpha[1,q,] 
          - t(Y[i,j-1,])%*%beta2[q,] - t(X[i,j-1,])%*%alpha[2,q,]- t(X[i,j-2,])%*%alpha[3,q,]
        Yy_sum <- Yy_sum + Y[i,j-2,]*c(Y_tilde_ijq)*Z[i,j,q]
      }
    }
    # Gaussian posterior distribution
    V_n <- chol2inv(chol(Y_sum/sigma2[q] + diag(1/(gamma_beta[3,q,]*nu_beta[3,q,]))))
    mu_n <- V_n %*% (Yy_sum/sigma2[q])
    beta3_update[q,] <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  }
  
  return(beta3_update)
}



update_gamma_beta <- function(beta, nu_beta, rho_beta){
  
  gamma_beta_update <- array(NA, dim=c(L+1, Q, Q))
  
  for (l in 1:(L+1)){
    for (q in 1:Q){
      for (p in 1:Q){
        if (l==1){
          # prob = pr(gamma_beta_lqp=1) / pr(gamma_beta_lqp=nu_0)
          prob <- (sqrt(1e-6)*rho_beta)/(1-rho_beta)*exp((1-1e-6)*beta[l,q,p]^2/(2*1e-6*nu_beta[l,q,p]))
          if (prob == Inf){ gamma_beta_update[l,q,p] <- 1
          }else{ gamma_beta_update[l,q,p] <- sample(c(1,1e-6), 1, prob=c(prob, 1)) }
        }else{
          # prob = pr(gamma_beta_lqp=1) / pr(gamma_beta_lqp=nu_0)
          prob <- (sqrt(nu_0)*rho_beta)/(1-rho_beta)*exp((1-nu_0)*beta[l,q,p]^2/(2*nu_0*nu_beta[l,q,p]))
          if (prob == Inf){ gamma_beta_update[l,q,p] <- 1
          }else{ gamma_beta_update[l,q,p] <- sample(c(1,nu_0), 1, prob=c(prob, 1)) }
        }
      }
    }
  }
  
  return(gamma_beta_update)
}



update_nu_beta <- function(beta, gamma_beta){
  
  nu_beta_update <- array(NA, dim=c(L+1, Q, Q))
  
  for (l in 1:(L+1)){
    for (q in 1:Q){
      for (p in 1:Q){
        a_nu_star <- a_nu + 1/2
        b_nu_star <- b_nu + beta[l,q,p]^2/(2*gamma_beta[l,q,p])
        nu_beta_update[l,q,p] <- rinvgamma(1, a_nu_star, b_nu_star)
      }
    }
  }
  
  return(nu_beta_update)
}



update_rho_beta <- function(gamma_beta){
  
  a_rho_star <- a_rho + sum(gamma_beta == 1)
  b_rho_star <- b_rho + sum(gamma_beta == nu_0)
  rho_beta_update <- rbeta(1, a_rho_star, b_rho_star)
  
  return(rho_beta_update)
}



update_alpha1 <- function(beta, alpha2, alpha3, gamma_alpha, nu_alpha, Z, mu, sigma2){
  
  alpha1_update <- matrix(NA, nrow=Q, ncol=S)
  
  for (q in 1:Q){
    # calculate X_sum and Xy_sum
    X_sum <- matrix(0, nrow=S, ncol=S); Xy_sum <- rep(0, S)
    for (i in 1:n){
      for (j in 1:J[i]){
        X_sum <- X_sum + X[i,j,]%*%t(X[i,j,])*Z[i,j,q]
        Y_tilde_ijq <- Y[i,j,q] - mu[q] - t(Y[i,j,])%*%beta[1,q,]
        if (j > 1) { Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-1,])%*%beta[2,q,] - t(X[i,j-1,])%*%alpha2[q,] }
        if (j > 2) { Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-2,])%*%beta[3,q,]- t(X[i,j-2,])%*%alpha3[q,] }
        Xy_sum <- Xy_sum + X[i,j,]*c(Y_tilde_ijq)*Z[i,j,q]
      }
    }
    # Gaussian posterior distribution
    V_n <- chol2inv(chol(X_sum/sigma2[q] + diag(1/(gamma_alpha[1,q,]*nu_alpha[1,q,]))))
    mu_n <- V_n %*% (Xy_sum/sigma2[q])
    alpha1_update[q,] <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  }
  
  return(alpha1_update)
}



update_alpha2 <- function(beta, alpha1, alpha3, gamma_alpha, nu_alpha, Z, mu, sigma2){
  
  alpha2_update <- matrix(NA, nrow=Q, ncol=S)
  
  for (q in 1:Q){
    # calculate X_sum and Xy_sum
    X_sum <- matrix(0, nrow=S, ncol=S); Xy_sum <- rep(0, S)
    for (i in 1:n){
      for (j in 2:J[i]){
        X_sum <- X_sum + X[i,j-1,]%*%t(X[i,j-1,])*Z[i,j,q]
        Y_tilde_ijq <- Y[i,j,q] - mu[q] - t(Y[i,j,])%*%beta[1,q,] - t(X[i,j,])%*%alpha1[q,] - t(Y[i,j-1,])%*%beta[2,q,]
        if (j > 2) { Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-2,])%*%beta[3,q,] - t(X[i,j-2,])%*%alpha3[q,] }
        Xy_sum <- Xy_sum + X[i,j-1,]*c(Y_tilde_ijq)*Z[i,j,q]
      }
    }
    # Gaussian posterior distribution
    V_n <- chol2inv(chol(X_sum/sigma2[q] + diag(1/(gamma_alpha[2,q,]*nu_alpha[2,q,]))))
    mu_n <- V_n %*% (Xy_sum/sigma2[q])
    alpha2_update[q,] <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  }
  
  return(alpha2_update)
}



update_alpha3 <- function(beta, alpha1, alpha2, gamma_alpha, nu_alpha, Z, mu, sigma2){
  
  alpha3_update <- matrix(NA, nrow=Q, ncol=S)
  
  for (q in 1:Q){
    # calculate X_sum and Xy_sum
    X_sum <- matrix(0, nrow=S, ncol=S); Xy_sum <- rep(0, S)
    for (i in 1:n){
      for (j in 3:J[i]){
        X_sum <- X_sum + X[i,j-2,]%*%t(X[i,j-2,])*Z[i,j,q]
        Y_tilde_ijq <- Y[i,j,q] - mu[q] - t(Y[i,j,])%*%beta[1,q,] - t(X[i,j,])%*%alpha1[q,] 
        - t(Y[i,j-1,])%*%beta[2,q,] - t(X[i,j-1,])%*%alpha2[q,] - t(Y[i,j-2,])%*%beta[3,q,]
        Xy_sum <- Xy_sum + X[i,j-2,]*c(Y_tilde_ijq)*Z[i,j,q]
      }
    }
    # Gaussian posterior distribution
    V_n <- chol2inv(chol(X_sum/sigma2[q] + diag(1/(gamma_alpha[3,q,]*nu_alpha[3,q,]))))
    mu_n <- V_n %*% (Xy_sum/sigma2[q])
    alpha3_update[q,] <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  }
  
  return(alpha3_update)
}



update_gamma_alpha <- function(alpha, nu_alpha, rho_alpha){
  
  gamma_alpha_update <- array(NA, dim=c(L+1, Q, S))
  
  for (l in 1:(L+1)){
    for (q in 1:Q){
      for (s in 1:S){
        # prob = pr(gamma_alpha_lqs=1) / pr(gamma_alpha_lqs=nu_0)
        prob <- (sqrt(nu_0)*rho_alpha)/(1-rho_alpha)*exp((1-nu_0)*alpha[l,q,s]^2/(2*nu_0*nu_alpha[l,q,s]))
        if (prob == Inf){ gamma_alpha_update[l,q,s] <- 1
        }else{ gamma_alpha_update[l,q,s] <- sample(c(1,nu_0), 1, prob=c(prob, 1)) }
      }
    }
  }
  
  return(gamma_alpha_update)
}



update_nu_alpha <- function(alpha, gamma_alpha){
  
  nu_alpha_update <- array(NA, dim=c(L+1, Q, S))
  
  for (l in 1:(L+1)){
    for (q in 1:Q){
      for (s in 1:S){
        a_nu_star <- a_nu + 1/2
        b_nu_star <- b_nu + alpha[l,q,s]^2/(2*gamma_alpha[l,q,s])
        nu_alpha_update[l,q,s] <- rinvgamma(1, a_nu_star, b_nu_star)
      }
    }
  }
  
  return(nu_alpha_update)
}



update_rho_alpha <- function(gamma_alpha){
  
  a_rho_star <- a_rho + sum(gamma_alpha == 1)
  b_rho_star <- b_rho + sum(gamma_alpha == nu_0)
  rho_alpha_update <- rbeta(1, a_rho_star, b_rho_star)
  
  return(rho_alpha_update)
}



update_Z <- function(beta, alpha, mu, sigma2){
  
  Z_update <- array(NA, dim=c(n, max(J), Q))
  eps <- 1e-3
  
  for (i in 1:n){
    for (j in 1:J[i]){
      for (q in 1:Q){
        Y_tilde_ijq <- Y[i,j,q] - mu[q] - t(Y[i,j,])%*%beta[1,q,] - t(X[i,j,])%*%alpha[1,q,]
        if (j > 1){ Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-1,])%*%beta[2,q,] - t(X[i,j-1,])%*%alpha[2,q,] }
        if (j > 2){ Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-2,])%*%beta[3,q,] - t(X[i,j-2,])%*%alpha[3,q,] }
        # Invere-Gaussian posterior distribution
        Z_update[i,j,q] <- max(rinvgaussian(1, mu=sqrt(sigma2[q])/(2*abs(Y_tilde_ijq)), lambda=1/4), eps)
      }
    }
  }
  
  return(Z_update)
}



update_mu <- function(beta, alpha, Z, sigma2){
  
  mu_update <- rep(NA, Q)
  
  for (q in 1:Q){
    
    Y_tilde_sum <- 0
    sigma2_tilde_sum <- 0
    
    for (i in 1:n){
      for (j in 1:J[i]){
        Y_tilde_ijq <- Y[i,j,q] - t(Y[i,j,])%*%beta[1,q,] - t(X[i,j,])%*%alpha[1,q,]
        if (j > 1) { Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-1,])%*%beta[2,q,]- t(X[i,j-1,])%*%alpha[2,q,] }
        if (j > 2) { Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-2,])%*%beta[3,q,]- t(X[i,j-2,])%*%alpha[3,q,] }
        Y_tilde_sum <- Y_tilde_sum + Y_tilde_ijq*Z[i,j,q]
        sigma2_tilde_sum <- sigma2_tilde_sum + Z[i,j,q]/sigma2[q]
      }
    }
    # Gaussian posterior distribution
    V_n <- 1/(sigma2_tilde_sum + 1/100)
    mu_n <- V_n*Y_tilde_sum/sigma2[q]
    mu_update[q] <- rnorm(1, mu_n, sqrt(V_n))
  }
  
  return(mu_update)
}



update_sigma2 <- function(beta, alpha, Z, mu){
  
  sigma2_update <- rep(NA, Q)
  
  for (q in 1:Q){
    Y_tilde_sum <- 0
    for (i in 1:n){
      for (j in 1:J[i]){
        Y_tilde_ijq <- Y[i,j,q] - mu[q] - t(Y[i,j,])%*%beta[1,q,] - t(X[i,j,])%*%alpha[1,q,]
        if (j > 1){ Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-1,])%*%beta[2,q,] - t(X[i,j-1,])%*%alpha[2,q,] }
        if (j > 2){ Y_tilde_ijq <- Y_tilde_ijq - t(Y[i,j-2,])%*%beta[3,q,] - t(X[i,j-2,])%*%alpha[3,q,] }
        Y_tilde_sum <- Y_tilde_sum + Y_tilde_ijq^2*Z[i,j,q]
      }
    }
    # Inverse-Gamma posterior distribution
    a_star <- a_sigma + sum(J)/2
    b_star <- b_sigma + Y_tilde_sum/2
    sigma2_update[q] <- rinvgamma(1, a_star, b_star)
  }
  
  return(sigma2_update)
}