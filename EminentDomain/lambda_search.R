#####################################################################
## This function iteratively adjusts penalty
## level to achieve desired sparsity
#####################################################################
## Input:
## - K: desired number of instruments
## Outpu:
## - L: penalty value found
## - IND: indices of the selected instruments (a logical vector)
#####################################################################

lambda_search <- function(X, y, K, Lambda0=0) {
  L   <- Lambda0
  PI  <- lasso_shooting(X, y, L)$b;
  IND <- abs(PI) > 1.0e-4
  k   <- sum(IND)
  
  iter <- 0
  swap <- 0
  
  while ( (k != K) & iter <= 1000) {
    direct <- sign(k-K)
    if (direct > 0) {
      L <- L + 5 * 10^(1-swap)
    } else {
      L <- L - 5 * 10^(1-swap)
    }
    PI <- lasso_shooting(X, y, L)$b
    IND <- abs(PI) > 1.0e-4
    k <- sum(IND)
    if (direct != sign(k-K)) {
      swap <- swap + 1
    }
    iter <- iter + 1
  }
  return(list(Lambda=L, Index=IND))
}