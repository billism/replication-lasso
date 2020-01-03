#############################################
# 2 Stage Least Squares Estimator
#############################################
## y: response variable
## d: endogenous variable
## x: included exogenous variables
## z: excluded exogenous variables
#############################################

tsls <- function(y, d, x, z) {
  n  <- length(y)
  k1 <- dim(d)[2]
  k2 <- dim(x)[2]
  if (is.vector(d)) {k1 <- 1}
  if (is.null(x))   {k2 <- 0}
  if (is.vector(x)) {k2 <- 1}
  k <- k1 + k2
  
  X  <- cbind(d, x)
  Z  <- cbind(z, x)

  XZ    <- crossprod(X, Z)
  ZZinv <- solve(crossprod(Z))
  M   <- solve(XZ %*% ZZinv %*% t(XZ))
  
  b   <- M %*% XZ %*% ZZinv %*% (t(Z) %*% y)
  e   <- y - X %*% b
  VC1 <- (crossprod(e) / (n-k)) * M
  
  return(list(b=b, VC1=VC1))
}