# Heteroskedasticity-Consistent Standard Errors

hetero_se <- function(x,e, xxinv) {
  k <- dim(x)[2]
  n <- dim(x)[1]
  
  V <- crossprod((x * ((e^2) %*% t(rep(1,k)))), x)
  
  vhetero <- xxinv %*% V %*% t(xxinv)
  se <- sqrt(diag(vhetero))
  
  return(list(se=se, vhetero=vhetero, V=V))
} 