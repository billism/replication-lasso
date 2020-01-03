# R Implementation of the Shooting Algorithm for LASSO
# Reference:
# 1. Fu (1998) Penalized Regressions: The Bridge Versus the Lasso
# 2. Murphy (2012) Machine Learning: A Probablilistic Perspective. MIT Press
# 3. Belloni et al (2014) High-Dimensional Methods and Inference on Treatment and Structural Effects in Economics

lasso_shooting <- function(X, y, lambda, tolerance=1e-5, max_iteration=10000) {
  p <- ncol(X)
  ## Precompute the Hessian diagonals, since they do not change
  XX <- crossprod(X)
  Xy <- crossprod(X,y)
  XX2 <- 2 * XX
  Xy2 <- 2 * Xy
  
  ## Belloni et al (2014) employed penalty loadings
  ## Ensure lambda is a vector of length p
  if (length(lambda) == 1) {
    lambda <- rep(lambda, p)
  }
  
  ## Initialize beta with the Ridge solution as suggested by Murphy
  beta <- solve(XX + diag(lambda, p), Xy)
  
  ## Iterations
  iteration <- 0
  converged <- FALSE
  while ( !converged & iteration < max_iteration) {
    beta_prev <- beta
    # optimize by coordinate
    for ( j in 1:p ) {
      s <- Xy2[j] - sum( XX2[j,] %*% beta ) + beta[j] * XX2[j,j]
      ## soft thresholding
      if ( abs(s) <= lambda[j] ) {
        beta[j] <- 0
      } else if ( s > lambda[j] ) {
        beta[j] <- (s - lambda[j])/XX2[j,j]
      } else {
        beta[j] <- (s + lambda[j])/XX2[j,j]
      }
    }
    iteration <- iteration + 1
    converged <- sum(abs(beta - beta_prev)) < tolerance
  }
  result <- list( b=beta, num_iter=iteration, converged=converged)
  return(result)
}
