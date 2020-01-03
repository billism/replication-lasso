# The Fuller (1977) (FULL) Estimator
## Reference:
#### 1. Hansen, Hausman, Newey (2008) Estimation with Many Instrumental Variables
#### 2. Anatolyev (2018) Many Instruments and/or Regressors: A Friendly Guide
################################################
# y: outcome variable in structural equation
# d: endogenous variable
# x: included exogenous variable
# z: excluded exogenous variable
# c: constant for Fuller estimator
################################################

fuller <- function(y, d, x=NULL, z, c=1) {
  n  <- length(y)
  k1 <- dim(d)[2]
  k2 <- dim(x)[2]
  if (is.vector(d)) {k1 <- 1}
  if (is.null(x))   {k2 <- 0}
  if (is.vector(x)) {k2 <- 1}
  k <- k1 + k2
  
  if (!is.matrix(z)) {z <- as.matrix(z)}
  if (!is.null(x)) {
    XXinv <- solve(crossprod(x))
    My <- y - x %*% XXinv %*% t(x) %*% y
    Md <- d - x %*% XXinv %*% t(x) %*% d
    Mz <- solve(crossprod(z) - crossprod(z, x) %*% XXinv %*% crossprod(x, z))
  } else {
    My <- y
    Md <- d
    Mz <- solve(crossprod(z))
  }
  
  ## Compute alpha
  Y <- cbind(My, Md)
  M <- solve(crossprod(Y)) %*% crossprod(Y, z) %*% Mz %*% crossprod(z, Y)
  #### alpha for LIML Estimator
  alpha  <- min(eigen(M)$values)
  #### alpha for Fuller Estimator
  alpha1 <- (alpha - (1 - alpha) * c / (n - k - dim(z)[2])) / (1 - (1 - alpha) * c / (n - k - dim(z)[2]))
  
  ## Compute Fuller Estimator
  X    <- cbind(d, x)
  Z    <- cbind(z, x)
  Mxy  <- crossprod(X, y)
  Mzy  <- crossprod(Z, y)
  Mzx  <- crossprod(Z, X)
  Mxx  <- crossprod(X)
  Mzzinv    <- solve(crossprod(Z))
  Mxzzzzx   <- t(Mzx) %*% Mzzinv %*% Mzx
  
  H  <- Mxzzzzx - alpha1 * Mxx
  b  <- solve(H) %*% (t(Mzx) %*% Mzzinv %*% Mzy - alpha1 * Mxy)
  e  <- y - X %*% b
  
  ## V1
  VC1 <- (crossprod(e)/(n - k)) * solve(H)
  
  ## V2
  J  <- Mxzzzzx - alpha1 * crossprod(X, e) %*% crossprod(e, X) %*% solve(crossprod(e))
  S  <- (1 - alpha1) * J - alpha1 * H
  VC2 <- (crossprod(e)/(n - k)) * solve(H) %*% S %*% solve(H)
  
  ## V3: This is the asymptotic variance estimator we will use
  ## Note: Its value can be quite difference than V1.
  s2      <- crossprod(e)/(n - k)
  a_tilde <- crossprod(e, Z) %*% Mzzinv %*% crossprod(Z, e) %*% solve(crossprod(e))
  Ups     <- Z %*% Mzzinv %*% crossprod(Z, X)
  X_hat   <- X - e %*% crossprod(e, X) %*% solve(crossprod(e))
  V_hat   <- X_hat - Z %*% Mzzinv %*% crossprod(Z, X_hat)
  tau     <- k/n
  kappa   <- 0
  A1      <- 0
  A2      <- 0
  B1      <- 0
  for (ii in 1:n) {
    pii <- Z[ii, , drop=FALSE] %*% Mzzinv %*% t(Z[ii, , drop=FALSE])
    kappa  <- kappa + pii
    A1  <- A1 + (pii - tau) %*% Ups[ii, , drop=FALSE]
    A2  <- A2 + (e[ii]^2) * V_hat[ii, , drop=FALSE] / n
    B1  <- B1 + ((e[ii]^2-s2) * (t(V_hat[ii, , drop=FALSE]) %*% V_hat[ii, , drop=FALSE]))
  }
  kappa <- kappa/k
  B     <- k * (kappa - tau) * B1 / (n * (1 - 2 * tau + kappa * tau))
  A     <- crossprod(A1, A2)
  SB    <- s2 * (((1-a_tilde)^2) %*% (crossprod(X_hat, Z) %*% Mzzinv %*% crossprod(Z, X_hat)) +
                (a_tilde^2) %*% (crossprod(X_hat) - crossprod(X_hat, Z) %*% Mzzinv %*% crossprod(Z, X_hat)))
  VC3   <- solve(Mxzzzzx - a_tilde %*% Mxx) %*% (SB + A + t(A) + B) %*% solve(Mxzzzzx - a_tilde %*% Mxx)
  
  return(list(b=b,VC1=VC1, VC2=VC2, VC3=VC3))
}