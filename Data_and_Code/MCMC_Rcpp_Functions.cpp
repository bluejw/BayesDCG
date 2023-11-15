//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp; 


const double log2pi = log(2.0*M_PI);


// [[Rcpp::export]]
double dmvn_rcpp(rowvec& x, rowvec& mean, mat& sigma, bool logd = false){ 
  
  // calculate density of multivariate normal distribution
  // args: x: row vector data
  //      mean: row vector mean, sigma: covariance matrix  
  //      logd: true for taking log
  // returns: out: pdf (or log pdf) of multivariate normal distribution
  
  int xdim = x.size(); 
  mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
  
  vec z = rooti*trans(x-mean);
  double out = constants-0.5*sum(z%z)+rootisum;
  
  if (logd == false){ out = exp(out); }
  return(out);
}


// [[Rcpp::export]]
mat rmvn_rcpp(const int n, vec& mean, mat& sigma){
  
  // randomly generate samples from multivariate normal distribution
  // args: n: number of data 
  //      mean: row vector mean, sigma: covariance matrix  
  // returns: out: random samples from multivariate normal distribution
  
  int k = sigma.n_cols; // dimension of the multivariate normal distribution
  mat z = randn(n, k);
  mat out = repmat(mean,1,n).t()+z*chol(sigma);
  return(out);
}


// [[Rcpp::export]]
double rinvgamma_rcpp(const double a, const double b){
  
  // generate random samples from inverse-gamma distribution
  // args: inverse-gamma(a, b)
  // returns: random sample from inverse-gamma distribution
  
  return(1/R::rgamma(a, 1/b));
}


// [[Rcpp::export]]
double rinvgaussian_rcpp(const double mu, const double lambda){
  
  // generate random samples from inverse-gaussian distribution
  // args: inverse-gaussian(mu, lambda)
  // returns: random sample from inverse-gaussian distribution
  
  double out;
  double nu = R::rnorm(0, 1);
  double y = pow(nu, 2);
  double x = mu + ((pow(mu,2)*y)/(2*lambda)) - (mu/(2*lambda))*
    sqrt(4*mu*lambda*y + pow(mu,2)*pow(y,2));
  double z = R::runif(0, 1);
  if (z > (mu/(mu + x))){ out = mu*mu/x; 
  }else{ out = x; }
  return(out);
}


// [[Rcpp::export]]
double logll_beta1_rcpp(const int q, mat& beta1, mat& beta2, mat& beta3, mat& alpha1, mat& alpha2, mat& alpha3,
                        cube& Z, vec& mu, vec& sigma2, const int n, vec& J, const int Q, const int S, cube& Y, cube& X){
  
  double logll = 0;
  double jacob_factor, Y_tilde_ijq;
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;   
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  
  // calculate the jacobian matrix 
  mat jacob; jacob.eye(Q,Q); jacob -= beta1;
  jacob_factor = abs(det(jacob));
  
  for (int i=0; i<n; i++){
    
    // extra factor: absolute determinant of jacobian matrix 
    logll += J(i)*log(jacob_factor);
    
    for (int j=0; j<J(i); j++){
      Y_tilde_ijq = Y(i,j,q) - mu(q);
      Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
      beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
      Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
      X_ij1 = X.subcube(i, j, 0, i, j, S-1); 
      alpha_1q = conv_to<vec>::from(alpha1.submat(q, 0, q, S-1));
      Y_tilde_ijq -= conv_to<double>::from(X_ij1 * alpha_1q);
      if (j > 0){
        Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1); 
        beta_2q = conv_to<vec>::from(beta2.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij2 * beta_2q);
        X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1); 
        alpha_2q = conv_to<vec>::from(alpha2.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij2 * alpha_2q);
      }
      if (j > 1){
        Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1); 
        beta_3q = conv_to<vec>::from(beta3.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij3 * beta_3q);
        X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1); 
        alpha_3q = conv_to<vec>::from(alpha3.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij3 * alpha_3q);
      }
      // calculate the log-likelihood of beta1q
      logll += R::dnorm(Y_tilde_ijq, 0, sqrt(sigma2(q)/Z(i,j,q)), true);
    }
  }
  
  return(logll);
}


// [[Rcpp::export]]
mat update_beta2_rcpp(mat& beta1, mat& beta3, mat& alpha1, mat& alpha2, mat& alpha3, 
                      mat& gamma_beta2, mat& nu_beta2, cube& Z, vec& mu, vec& sigma2,
                      const int n, vec& J, const int Q, const int S, cube& Y, cube& X){
  
  mat beta2_update(Q, Q);
  vec Yy_sum(Q); mat Y_sum(Q,Q);
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  double Y_tilde_ijq;
  vec mu_n(Q); mat V_n(Q,Q);
  
  for (int q=0; q<Q; q++){
    
    Y_sum.fill(0); Yy_sum.fill(0);
    
    for (int i=0; i<n; i++){
      for (int j=1; j<J(i); j++){
        // calculate Y_sum
        Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1);
        Y_sum += Y_ij2.t()*Y_ij2*Z(i,j,q);
        // calculate Yy_sum
        Y_tilde_ijq = Y(i,j,q) - mu(q);
        Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
        beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
        X_ij1 = X.subcube(i, j, 0, i, j, S-1); 
        alpha_1q = conv_to<vec>::from(alpha1.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij1 * alpha_1q);
        X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1); 
        alpha_2q = conv_to<vec>::from(alpha2.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij2 * alpha_2q);
        if (j > 1){
          Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1); 
          beta_3q = conv_to<vec>::from(beta3.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij3 * beta_3q);
          X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1); 
          alpha_3q = conv_to<vec>::from(alpha3.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij3 * alpha_3q);
        }
        Yy_sum += conv_to<vec>::from(Y_ij2*Y_tilde_ijq*Z(i,j,q));
      }
    }
    // Gaussian posterior distribution
    mat V_beta; V_beta.eye(Q,Q);
    for (int q=0; q<Q; q++){ V_beta(q,q) *= 1/(gamma_beta2(q,q)*nu_beta2(q,q)); }
    V_n = inv_sympd(Y_sum/sigma2(q) + V_beta);
    mu_n = V_n*(Yy_sum/sigma2(q));
    beta2_update.row(q) = conv_to<rowvec>::from(rmvn_rcpp(1, mu_n, V_n));
  }
  
  return(beta2_update);
}


// [[Rcpp::export]]
mat update_beta3_rcpp(mat& beta1, mat& beta2, mat& alpha1, mat& alpha2, mat& alpha3, 
                      mat& gamma_beta3, mat& nu_beta3, cube& Z, vec& mu, vec& sigma2,
                      const int n, vec& J, const int Q, const int S, cube& Y, cube& X){
  
  mat beta3_update(Q, Q);
  vec Yy_sum(Q); mat Y_sum(Q,Q);
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  double Y_tilde_ijq;
  vec mu_n(Q); mat V_n(Q,Q);
  
  for (int q=0; q<Q; q++){
    
    Y_sum.fill(0); Yy_sum.fill(0);
    
    for (int i=0; i<n; i++){
      for (int j=2; j<J(i); j++){
        // calculate Y_sum
        Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1);
        Y_sum += Y_ij3.t()*Y_ij3*Z(i,j,q);
        // calculate Yy_sum
        Y_tilde_ijq = Y(i,j,q) - mu(q);
        Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
        beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
        X_ij1 = X.subcube(i, j, 0, i, j, S-1); 
        alpha_1q = conv_to<vec>::from(alpha1.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij1 * alpha_1q);
        Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1); 
        beta_2q = conv_to<vec>::from(beta2.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij2 * beta_2q);
        X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1); 
        alpha_2q = conv_to<vec>::from(alpha2.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij2 * alpha_2q);
        X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1); 
        alpha_3q = conv_to<vec>::from(alpha3.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij3 * alpha_3q);
        Yy_sum += conv_to<vec>::from(Y_ij3*Y_tilde_ijq*Z(i,j,q));
      }
    }
    // Gaussian posterior distribution
    mat V_beta; V_beta.eye(Q,Q);
    for (int q=0; q<Q; q++){ V_beta(q,q) *= 1/(gamma_beta3(q,q)*nu_beta3(q,q)); }
    V_n = inv_sympd(Y_sum/sigma2(q) + V_beta);
    mu_n = V_n*(Yy_sum/sigma2(q));
    beta3_update.row(q) = conv_to<rowvec>::from(rmvn_rcpp(1, mu_n, V_n));
  }
  
  return(beta3_update);
}


// [[Rcpp::export]]
mat update_alpha1_rcpp(mat& beta1, mat& beta2, mat& beta3, mat& alpha2, mat& alpha3, 
                       mat& gamma_alpha1, mat& nu_alpha1, cube& Z, vec& mu, vec& sigma2,
                       const int n, vec& J, const int Q, const int S, cube& Y, cube& X){
  
  mat alpha1_update(Q, S);
  vec Xy_sum(S); mat X_sum(S,S);
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  double Y_tilde_ijq;
  vec mu_n(S); mat V_n(S,S);
  
  for (int q=0; q<Q; q++){
    
    X_sum.fill(0); Xy_sum.fill(0);
    
    for (int i=0; i<n; i++){
      for (int j=0; j<J(i); j++){
        // calculate X_sum
        X_ij1 = X.subcube(i, j, 0, i, j, S-1);
        X_sum += X_ij1.t()*X_ij1*Z(i,j,q);
        // calculate Xy_sum
        Y_tilde_ijq = Y(i,j,q) - mu(q);
        Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
        beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
        if (j > 0){
          Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1); 
          beta_2q = conv_to<vec>::from(beta2.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij2 * beta_2q);
          X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1); 
          alpha_2q = conv_to<vec>::from(alpha2.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij2 * alpha_2q);
        }
        if (j > 1){
          Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1); 
          beta_3q = conv_to<vec>::from(beta3.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij3 * beta_3q);
          X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1); 
          alpha_3q = conv_to<vec>::from(alpha3.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij3 * alpha_3q);
        }
        Xy_sum += conv_to<vec>::from(X_ij1*Y_tilde_ijq*Z(i,j,q));
      }
    }
    // Gaussian posterior distribution
    mat V_alpha; V_alpha.eye(S,S);
    for (int s=0; s<S; s++){ V_alpha(s,s) *= 1/(gamma_alpha1(q,s)*nu_alpha1(q,s)); }
    V_n = inv_sympd(X_sum/sigma2(q) + V_alpha);
    mu_n = V_n*(Xy_sum/sigma2(q));
    alpha1_update.row(q) = conv_to<rowvec>::from(rmvn_rcpp(1, mu_n, V_n));
  }
  
  return(alpha1_update);
}


// [[Rcpp::export]]
mat update_alpha2_rcpp(mat& beta1, mat& beta2, mat& beta3, mat& alpha1, mat& alpha3, 
                       mat& gamma_alpha2, mat& nu_alpha2, cube& Z, vec& mu, vec& sigma2,
                       const int n, vec& J, const int Q, const int S, cube& Y, cube& X){
  
  mat alpha2_update(Q, S);
  vec Xy_sum(S); mat X_sum(S,S);
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  double Y_tilde_ijq;
  vec mu_n(S); mat V_n(S,S);
  
  for (int q=0; q<Q; q++){
    
    X_sum.fill(0); Xy_sum.fill(0);
    
    for (int i=0; i<n; i++){
      for (int j=1; j<J(i); j++){
        // calculate X_sum
        X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1);
        X_sum += X_ij2.t()*X_ij2*Z(i,j,q);
        // calculate Xy_sum
        Y_tilde_ijq = Y(i,j,q) - mu(q);
        Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
        beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
        X_ij1 = X.subcube(i, j, 0, i, j, S-1); 
        alpha_1q = conv_to<vec>::from(alpha1.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij1 * alpha_1q);
        Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1); 
        beta_2q = conv_to<vec>::from(beta2.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij2 * beta_2q);
        if (j > 1){
          Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1); 
          beta_3q = conv_to<vec>::from(beta3.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij3 * beta_3q);
          X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1); 
          alpha_3q = conv_to<vec>::from(alpha3.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij3 * alpha_3q);
        }
        Xy_sum += conv_to<vec>::from(X_ij2*Y_tilde_ijq*Z(i,j,q));
      }
    }
    // Gaussian posterior distribution
    mat V_alpha; V_alpha.eye(S,S);
    for (int s=0; s<S; s++){ V_alpha(s,s) *= 1/(gamma_alpha2(q,s)*nu_alpha2(q,s)); }
    V_n = inv_sympd(X_sum/sigma2(q) + V_alpha);
    mu_n = V_n*(Xy_sum/sigma2(q));
    alpha2_update.row(q) = conv_to<rowvec>::from(rmvn_rcpp(1, mu_n, V_n));
  }
  
  return(alpha2_update);
}


// [[Rcpp::export]]
mat update_alpha3_rcpp(mat& beta1, mat& beta2, mat& beta3, mat& alpha1, mat& alpha2, 
                       mat& gamma_alpha3, mat& nu_alpha3, cube& Z, vec& mu, vec& sigma2,
                       const int n, vec& J, const int Q, const int S, cube& Y, cube& X){
  
  mat alpha3_update(Q, S);
  vec Xy_sum(S); mat X_sum(S,S);
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  double Y_tilde_ijq;
  vec mu_n(S); mat V_n(S,S);
  
  for (int q=0; q<Q; q++){
    
    X_sum.fill(0); Xy_sum.fill(0);
    
    for (int i=0; i<n; i++){
      for (int j=2; j<J(i); j++){
        // calculate X_sum
        X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1);
        X_sum += X_ij3.t()*X_ij3*Z(i,j,q);
        // calculate Xy_sum
        Y_tilde_ijq = Y(i,j,q) - mu(q);
        Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
        beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
        X_ij1 = X.subcube(i, j, 0, i, j, S-1); 
        alpha_1q = conv_to<vec>::from(alpha1.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij1 * alpha_1q);
        Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1); 
        beta_2q = conv_to<vec>::from(beta2.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij2 * beta_2q);
        X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1); 
        alpha_2q = conv_to<vec>::from(alpha2.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij2 * alpha_2q);
        Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1); 
        beta_3q = conv_to<vec>::from(beta3.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij3 * beta_3q);
        Xy_sum += conv_to<vec>::from(X_ij2*Y_tilde_ijq*Z(i,j,q));
      }
    }
    // Gaussian posterior distribution
    mat V_alpha; V_alpha.eye(S,S);
    for (int s=0; s<S; s++){ V_alpha(s,s) *= 1/(gamma_alpha3(q,s)*nu_alpha3(q,s)); }
    V_n = inv_sympd(X_sum/sigma2(q) + V_alpha);
    mu_n = V_n*(Xy_sum/sigma2(q));
    alpha3_update.row(q) = conv_to<rowvec>::from(rmvn_rcpp(1, mu_n, V_n));
  }
  
  return(alpha3_update);
}


// [[Rcpp::export]]
cube update_Z_rcpp(mat& beta1, mat& beta2, mat& beta3, mat& alpha1, mat& alpha2, mat& alpha3, 
                   vec& mu, vec& sigma2, const int n, vec& J, const int Q, const int S, cube& Y, cube& X){
    
  cube Z_update(n, max(J), Q); Z_update.fill(datum::nan);
  double eps = 1e-3; double Y_tilde_ijq;
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  
  for (int i=0; i<n; i++){
    for (int j=0; j<J(i); j++){
      for (int q=0; q<Q; q++){
        Y_tilde_ijq = Y(i,j,q) - mu(q);
        Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
        beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
        X_ij1 = X.subcube(i, j, 0, i, j, S-1); 
        alpha_1q = conv_to<vec>::from(alpha1.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij1 * alpha_1q);
        if (j > 0){
          Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1); 
          beta_2q = conv_to<vec>::from(beta2.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij2 * beta_2q);
          X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1); 
          alpha_2q = conv_to<vec>::from(alpha2.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij2 * alpha_2q);
        }
        if (j > 1){
          Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1); 
          beta_3q = conv_to<vec>::from(beta3.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij3 * beta_3q);
          X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1); 
          alpha_3q = conv_to<vec>::from(alpha3.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij3 * alpha_3q);
        }
        // Inverse-Gaussian posterior distribution
        double mu_ijq = sqrt(sigma2(q))/(2*abs(Y_tilde_ijq));
        double lambda_ijq = 1.0/4;
        double Z_update_ijq = rinvgaussian_rcpp(mu_ijq, lambda_ijq);
        Z_update(i,j,q) = max(Z_update_ijq, eps);
      }
    }
  }
  
  return(Z_update);
}


// [[Rcpp::export]]
vec update_mu_rcpp(mat& beta1, mat& beta2, mat& beta3, mat& alpha1, mat& alpha2, mat& alpha3, 
                   cube& Z, vec& sigma2, const int n, vec& J, const int Q, const int S, cube& Y, cube& X){
  
  vec mu_update(Q);
  double Y_tilde_sum, sigma2_tilde_sum;
  double Y_tilde_ijq, mu_n, V_n;
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  
  for (int q=0; q<Q; q++){
    
    Y_tilde_sum = 0; sigma2_tilde_sum = 0;
    
    for (int i=0; i<n; i++){
      for (int j=0; j<J(i); j++){
        // calcualte Y_tilde_sum and sigma2_tilde_sum 
        Y_tilde_ijq = Y(i,j,q);
        Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
        beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
        X_ij1 = X.subcube(i, j, 0, i, j, S-1); 
        alpha_1q = conv_to<vec>::from(alpha1.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij1 * alpha_1q);
        if (j > 0){
          Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1); 
          beta_2q = conv_to<vec>::from(beta2.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij2 * beta_2q);
          X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1); 
          alpha_2q = conv_to<vec>::from(alpha2.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij2 * alpha_2q);
        }
        if (j > 1){
          Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1); 
          beta_3q = conv_to<vec>::from(beta3.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij3 * beta_3q);
          X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1); 
          alpha_3q = conv_to<vec>::from(alpha3.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij3 * alpha_3q);
        }
        Y_tilde_sum += Y_tilde_ijq*Z(i,j,q);
        sigma2_tilde_sum += Z(i,j,q)/sigma2(q);
      }
    }
    // Gaussian posterior distribution
    V_n = 1/(sigma2_tilde_sum + 1/100);
    mu_n = V_n*Y_tilde_sum/sigma2(q);
    mu_update(q) = R::rnorm(mu_n, sqrt(V_n));
  }
  
  return(mu_update);
}


// [[Rcpp::export]]
vec update_sigma2_rcpp(mat& beta1, mat& beta2, mat& beta3, mat& alpha1, mat& alpha2, mat& alpha3, 
                       cube& Z, vec& mu, const int n, vec& J, const int Q, const int S, cube& Y, cube& X,
                       const double a_sigma, const double b_sigma){
  
  vec sigma2_update(Q);
  double Y_tilde_sum, Y_tilde_ijq;
  double a_sigma_star, b_sigma_star;
  rowvec Y_ij1, Y_ij2, Y_ij3, X_ij1, X_ij2, X_ij3;
  vec beta_1q, beta_2q, beta_3q, alpha_1q, alpha_2q, alpha_3q;
  
  for (int q=0; q<Q; q++){
    
    Y_tilde_sum = 0;
    
    for (int i=0; i<n; i++){
      for (int j=0; j<J(i); j++){
        // calcualte Y_tilde_sum 
        Y_tilde_ijq = Y(i,j,q) - mu(q);
        Y_ij1 = Y.subcube(i, j, 0, i, j, Q-1); 
        beta_1q = conv_to<vec>::from(beta1.submat(q, 0, q, Q-1));
        Y_tilde_ijq -= conv_to<double>::from(Y_ij1 * beta_1q);
        X_ij1 = X.subcube(i, j, 0, i, j, S-1); 
        alpha_1q = conv_to<vec>::from(alpha1.submat(q, 0, q, S-1));
        Y_tilde_ijq -= conv_to<double>::from(X_ij1 * alpha_1q);
        if (j > 0){
          Y_ij2 = Y.subcube(i, j-1, 0, i, j-1, Q-1); 
          beta_2q = conv_to<vec>::from(beta2.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij2 * beta_2q);
          X_ij2 = X.subcube(i, j-1, 0, i, j-1, S-1); 
          alpha_2q = conv_to<vec>::from(alpha2.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij2 * alpha_2q);
        }
        if (j > 1){
          Y_ij3 = Y.subcube(i, j-2, 0, i, j-2, Q-1); 
          beta_3q = conv_to<vec>::from(beta3.submat(q, 0, q, Q-1));
          Y_tilde_ijq -= conv_to<double>::from(Y_ij3 * beta_3q);
          X_ij3 = X.subcube(i, j-2, 0, i, j-2, S-1); 
          alpha_3q = conv_to<vec>::from(alpha3.submat(q, 0, q, S-1));
          Y_tilde_ijq -= conv_to<double>::from(X_ij3 * alpha_3q);
        }
        Y_tilde_sum += pow(Y_tilde_ijq,2)*Z(i,j,q);
      }
    }
    // Inverse-Gamma posterior distribution
    a_sigma_star = a_sigma + accu(J)/2;
    b_sigma_star = b_sigma + Y_tilde_sum/2;
    sigma2_update(q) = rinvgamma_rcpp(a_sigma_star, b_sigma_star);
  }
  
  return(sigma2_update);
}