// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List lasso_shooting_cpp(mat X, vec y, vec lambda, double tolerance = 1e-5, int max_iteration = 10000){
  int p = X.n_cols;
  // Precompute the Hessian diagonals, since they do not change
  mat XX = X.t() * X;
  vec Xy = X.t() * y;
  vec Xy2 = 2 * Xy;
  mat XX2 = 2 * XX;
  
  // Initialize beta with the Ridge solution as suggested by Murphy
  vec beta = solve(XX + diagmat(lambda), Xy);
  
  // Iterations
  int iteration = 0;
  bool converged = false;
  
  vec beta_prev, s;
  
  while (!converged && (iteration < max_iteration)){
    beta_prev = beta;
    for (int j = 0; j < p; j++){
      s = Xy2(j) - dot(XX2.row(j), beta) + beta(j) * XX2(j,j);
      // soft thresholding
      if ( abs(as_scalar(s)) <= lambda(j) ) {
        beta(j) = 0;
      } else if ( as_scalar(s) > lambda(j) ) {
        beta(j) = as_scalar((s - lambda(j))/XX2(j,j));
      } else {
        beta(j) = as_scalar((s + lambda(j))/XX2(j,j));
      }
    }
    iteration = iteration + 1;
    converged =  norm(beta_prev - beta, 1) < tolerance;  
  }
  return List::create(Named("b") = beta,
                      Named("num_iter") = iteration,
                      Named("converged") = converged);
}