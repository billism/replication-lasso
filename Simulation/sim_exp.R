# Simulation : Exponential Design with pi=0.7
#########################################
## 2SLS
## Fuller
## Post-Lasso
## Post-Lasso-F
## Post_Lasso (Ridge)
## Post_Lasso-F (Ridge)
#########################################
options(warn=-1)
## Clear Workspace
rm(list = ls(all.names = TRUE))
gc()
## Library
library(Rcpp)
library(RcppArmadillo)
library(glmnet)
library(openxlsx)
source("./lasso_shooting.R")
source("./tsls.R")
source("./fuller.R")
sourceCpp("./lasso_shooting.cpp")

#########################################
#########      Parameters      ##########
#########################################

nSims <- 500          # Number of simulations
ns    <- c(100, 250)  # Sample sizes
nn    <- length(ns)
p     <- 100          # number of instruments
cps   <- c(30, 180)   # concentration parameters
ncp   <- length(cps)
s2z   <- 1            # /sigma_{z_h}
Sz    <- toeplitz(0.5^(0:(p-1)))  # /Sigma_{Z}
s2e   <- 1            # /sigma_{e}
cev   <- 0.6          # corr(e,v)
K     <- 15           # Number of Iterations to determine penalty loadings
## Exponential Design
#s        <- 5            # cutoff design with s=5
#PI_tild  <- c(rep(1, s), rep(0, p-s))            # (1,1,...,0,0,0,...)
pi       <- 0.7          # exponential design with pi=0.7
PI_tild  <- pi^(1:p - 1) # (1, 0.7, 0.7^2, 0.7^3, ... , 0.7^99)

#### Generate a load of matrices for storing results
source("./result_init.R")

###########################################
########      Simulation     ##############
###########################################

## sample size: 100, 250
## concentration parameter: 30, 180
## 4 combinations
for (ii in 1:nn) {
  for (jj in 1:ncp) {
    n   <-  ns[ii]
    cp  <- cps[jj]
    cat(paste0("=========================\n"))
    cat(paste0("Sample Size:        ", n, "\n"))
    cat(paste0("Concentration Para: ", cp, "\n"))
    cat(paste0("-------------------------\n"))
    ## Given n and cp, compute Variance-Covariance Matrix of e and v
    C        <- sqrt(cp / ((cp + n) * t(PI_tild) %*% Sz %*% PI_tild))   #Compute C given concentration
    #parameter and sample size
    PI       <- C * PI_tild                                             #First-stage coefficients
    s2v      <- 1 - C^2 * (t(PI_tild) %*% Sz %*% PI_tild)               # /sigma_{v} var(v)
    sev      <- cev * sqrt(s2e * s2v)                                   # /sigma_{e,v} cov(e,v)
    SU       <- matrix(cbind(s2e, sev, sev, s2v), 2)                    # /Sigma_{e,v}
    # variance-covariance matrix
    ## Some common components across simulations
    #### lambda for iteratively computing penalty loadings
    #### in Lasso and Lasso (Ridge)
    lambda0  <- 2 * 1.1 * sqrt(2 * log(2 * p * log(n)/.1))
    Rlambda0 <- 2 * 1.1 * sqrt(2 * log(2 * (p+1) * log(n/2)/.1))
    #### Alternatively: the following formula is given in Belloni et al (2012) Appendix A,
    #### the above employed version is just an approximation of the following one.
    # lambda0  <- 2 * 1.1 * sqrt(n) * qnorm(1 - (0.1/log(p))/(2 * p))   !!!TODO
    
    ## Start Simulation Loop
    runtime <- proc.time()
    for (kk in 1:nSims) {
      ## Print current simulation index for every 10 iterations
      if (kk %% 10 == 0) {cat(paste0(">> Current Repetition: ", kk, "\n"))}
      
      #########################################################
      ##  1. Generating Data according to the following DGP  ##
      #########################################################
      #### yi = beta*xi + ei
      #### xi = zi * PI + vi
      alpha  <- 1
      Z_Raw <- matrix(rnorm(n*p), n, p) %*% chol(Sz)   # generating z (nxp)
      U     <- matrix(rnorm(n*2), ncol=2) %*% chol(SU) # generating e, v (nx2)
      X_Raw <- Z_Raw %*% PI +  U[,2, drop=FALSE]       # generating x (nx1)
      Y_Raw <- alpha * X_Raw + U[,1, drop=FALSE]       # generating y (nx1)
      #### Demean
      Z     <- Z_Raw - colMeans(Z_Raw)
      X     <- X_Raw - mean(X_Raw)
      Y     <- Y_Raw - mean(Y_Raw)
      
      #######################################################
      ######          2. Estimation          ################
      #######################################################
      
      ##############################################
      # 2.1 2SLS and Fuller Estimation
      ##############################################
      if (p >= n) {
        #### If p >= n, randomly choose n-2 columns from Z
        p_use  <- sample(p, n-2)
        Z_use  <- Z[, p_use, drop=FALSE]
      } else {
        Z_use  <- Z
      }
      
      out_tsls <- tsls(Y, X, NULL, Z_use)
      out_full <- fuller(Y, X, NULL, Z_use)       # with default c=1
      
      b2sls[kk,ii,jj] = out_tsls$b                # Coef
      s2sls[kk,ii,jj] = sqrt(out_tsls$VC1)        # Standard Error
      bfull[kk,ii,jj] = out_full$b
      sfull[kk,ii,jj] = sqrt(out_full$VC3)
      
      ##############################################
      # 2.2 LASSO Estimation
      ##############################################
      #### Highest Correlation 2SLS
      Z1       <- Z[, which.max(cor(Z,X)), drop=FALSE]    # X is nx1
      out_tsls <- tsls(Y, X, NULL, Z1)
      
      b2sls1[kk,ii,jj] <- out_tsls$b
      s2sls1[kk,ii,jj] <- sqrt(out_tsls$VC1)
      
      #### Iteratively Choose Penalty Loadings
      v0    <- X
      pl    <- sqrt(crossprod(v0^2, Z^2))     # Initial Penalty Loadings
      b_tmp <- lasso_shooting_cpp(Z, X, lambda0 * pl)$b
      ind0  <- abs(b_tmp) > 0
      Z0    <- as.matrix(Z[, ind0])
      ind1  <- ind0
      #### It's possible the Lasso selects 0 instruments
      for (mm in 1:K) {
        if (sum(ind1) == 0) {break}
        v1    <- X - Z0 %*% lm(X ~ -1 + Z0)$coef
        pl    <- sqrt(crossprod(v1^2, Z^2))   # Update Penalty Loadings
        b_tmp <- lasso_shooting_cpp(Z, X, lambda0 * pl)$b
        ind1  <- abs(b_tmp) > 0
        Z0    <- as.matrix(Z[, ind1])
      }
      ### ind1 : A Logical Vector indicates the Instrument(s) Selected
      ### Z1   : Selected Instrument(s)
      if (sum(ind1 != 0)) { Z1 <- Z0}
      
      ####### Post-Lasso Estimation
      if ( sum(ind1) == 0 ) {
        blassoC[kk,ii,jj]  <- NaN
        slassoC[kk,ii,jj]  <- NaN
        blassoCF[kk,ii,jj] <- NaN
        slassoCF[kk,ii,jj] <- NaN
      } else {
        ### Run 2SLS with Selected Instruments
        bfs  <- lm(X ~ -1 + Z1)$coef
        efs  <- X - Z1 %*% as.matrix(bfs)
        ### Asymptotic Variance Estimator: Belloni et al (2012) 2.11
        Vfs  <- solve(crossprod(Z1)) %*% (t(Z1 * (efs^2 %*% matrix(rep(1,sum(ind1)), nrow=1))) %*% Z1) %*% solve(crossprod(Z1))
        ### F-Statistic
        FS[kk,ii,jj] = (t(bfs) %*% solve(Vfs) %*% bfs) / sum(ind1)
        
        out_2sls <- tsls(Y, X, NULL, Z1)
        out_full <- fuller(Y, X, NULL, Z1)
        blassoC[kk,ii,jj]  <- out_2sls$b
        slassoC[kk,ii,jj]  <- sqrt(out_2sls$VC1)
        blassoCF[kk,ii,jj] <- out_full$b
        slassoCF[kk,ii,jj] <- sqrt(out_full$VC3) 
      }
      
      #### Choose penalty by cross-validation
      #### alpha=1 for lasso; alpha=0 for ridge
      # cv <- cv.glmnet(x=Z,y=X, alpha=1, nfolds=10)
      # lambda <- rep(cv$lambda.min, p)      
      # b_tmp <- lasso_shooting_cpp(Z,X, lambda)$b
      # ind0 <- (abs(b_tmp) > 0)
      # Z0 <- as.matrix(Z[, ind0, drop=FALSE])
      # Z1 <- Z0
      # ind1 <- ind0
      # 
      # if ( sum(ind1) == 0) {
      #   blassoCV[kk,ii,jj]  <- NaN
      #   slassoCV[kk,ii,jj]  <- NaN
      #   blassoCFV[kk,ii,jj] <- NaN
      #   slassoCFV[kk,ii,jj] <- NaN
      #   FSV[kk,ii,jj] <- 0 
      # } else {
      #   bfs <- lm(X ~ -1 + Z1)$coef
      #   efs <- X - Z1 %*% as.matrix(bfs)
      #   Vfs  <- solve(crossprod(Z1)) %*% (t(Z1 * (efs^2 %*% matrix(rep(1,sum(ind1)), nrow=1))) %*% Z1) %*% solve(crossprod(Z1))
      #   FSV[kk,ii,jj] = t(bfs) %*% solve(Vfs) %*% bfs / sum(ind1)
      #   out_tsls <- tsls(Y,X,NULL,Z1)
      #   out_ful  <- fuller(Y,X,NULL,Z1)
      #   blassoCV[kk,ii,jj]  <- out_tsls$b
      #   slassoCV[kk,ii,jj]  <- sqrt(out_tsls$VC1)
      #   blassoCFV[kk,ii,jj] <- out_full$b
      #   slassoCFV[kk,ii,jj] <-  sqrt(out_full$VC3)      
      # }
      
      #### X-Dependent Penalty
      R <- 500 # number of simulations
      sim <- vector("numeric", length=R)
      for (l in 1:R) {
        g <- matrix(rep(rnorm(n), each=p), ncol=p, byrow=TRUE)
        sim[l] <- n * max(2 * colMeans(Z * g))
      }
      sigma0 <- sqrt(var(X))
      lambda <- rep(1.1*quantile(sim, probs=1-0.05)*sigma0, p)
      b_tmp <- lasso_shooting_cpp(Z,X, lambda)$b
      ind0 <- abs(b_tmp) > 0
      Z0 <- as.matrix(Z[,ind0])
      Z1 <- Z0
      ind1 <- ind0
      if ( sum(ind1) == 0) {
        blassoCn[kk,ii,jj]  <- NaN
        slassoCn[kk,ii,jj]  <- NaN
        blassoCFn[kk,ii,jj] <- NaN
        slassoCFn[kk,ii,jj] <- NaN
        FSn[kk,ii,jj] <- 0 
      } else {
        bfs <- lm(X ~ -1 + Z1)$coef
        efs <- X - Z1 %*% as.matrix(bfs)
        Vfs <- solve(crossprod(Z1)) %*% (t(Z1 * (efs^2 %*% matrix(rep(1,sum(ind1)), nrow=1))) %*% Z1) %*% solve(crossprod(Z1))
        FSn[kk,ii,jj] = t(bfs) %*% solve(Vfs) %*% bfs / sum(ind1)
        out_tsls <- tsls(Y,X,NULL,Z1)
        out_full<- fuller(Y,X,NULL,Z1)
        blassoCn[kk,ii,jj]  <-  out_tsls$b
        slassoCn[kk,ii,jj]  <- sqrt(out_tsls$VC1)
        blassoCFn[kk,ii,jj] <- out_full$b
        slassoCFn[kk,ii,jj] <-  sqrt(out_full$VC3)      
      }
      
      
      # X independent choice
      lambda <- rep(2*1.1*sqrt(n)*qnorm(1-0.1/(2*p))*sqrt(var(X)),p)   
      b_tmp  <- lasso_shooting_cpp(Z,X, lambda)$b
      ind0   <- abs(b_tmp) > 0
      Z0     <- as.matrix(Z[,ind0])
      ind1   <- ind0
      Z1     <- Z0
      if ( sum(ind1) == 0) {
        blassoCX[kk,ii,jj]  <- NaN
        slassoCX[kk,ii,jj]  <- NaN
        blassoCFX[kk,ii,jj] <- NaN
        slassoCFX[kk,ii,jj] <- NaN
        FSX[kk,ii,jj] <- 0 
      } else {
        bfs <- lm(X ~ -1 + Z1)$coef
        efs <- X - Z1%*%as.matrix(bfs)
        Vfs <- solve(crossprod(Z1)) %*% (t(Z1 * (efs^2 %*% matrix(rep(1,sum(ind1)), nrow=1))) %*% Z1) %*% solve(crossprod(Z1))
        FSX[kk,ii,jj] = t(bfs) %*% solve(Vfs) %*% bfs / sum(ind1)
        out_tsls <- tsls(Y,X, NULL,Z1)
        out_full <- fuller(Y,X, NULL,Z1)
        blassoCX[kk,ii,jj] <-  out_tsls$b
        slassoCX[kk,ii,jj] <- sqrt(out_tsls$VC1)
        blassoCFX[kk,ii,jj] <- out_full$b
        slassoCFX[kk,ii,jj] <-  sqrt(out_full$VC3)      
      }    
      
      ### Sup-Score Test
      for (mm in 1:na) {
        aEval <- aTest[mm]
        eTmp  <- Y-aEval*X
        ScoreVec    <- t(eTmp) %*% Z
        ScoreStd    <- sqrt(crossprod(eTmp^2, Z^2))
        ScaledScore <- ScoreVec / (1.1 * ScoreStd)
        supScore05[kk,mm,ii,jj] <- max(abs(ScaledScore)) < lambdaSS05
        supScore[kk,mm,ii,jj]   <- max(abs(ScaledScore))    # Sup-Score Test Statistic
      }
      
      #################################################################
      # Split Sample - Ridge - Lasso Estimation
      #################################################################
      ## Split Sample
      sss    <- n %/% 2         # Split Sample Size
      index  <- 1:n
      UseRidgeA <- sample(n, n/2)
      UseRidgeB <- index[-UseRidgeA]
      
      YA <- Y_Raw[UseRidgeA]
      XA <- X_Raw[UseRidgeA]
      ZA <- Z_Raw[UseRidgeA, , drop=FALSE]
      YB <- Y_Raw[UseRidgeB]
      XB <- X_Raw[UseRidgeB]
      ZB <- Z_Raw[UseRidgeB, , drop=FALSE]
      
      #### Leave-One-Out Cross Validation for Ridge Penalty
      #### Use glmnet Library for Better Performance
      #### Set alpha=0 for Ridge
      #### Note: Use Leave-one-out CV with Large Sample Size will be a lot slower,
      ####       Be careful for test purpose.
      cv.ridgeA <- cv.glmnet(x=ZA, y=XA, alpha=0, nfolds=10, lambda=10^seq(1,-4,-.05))
      LambdaRidgeA[kk,ii,jj] <- cv.ridgeA$lambda.min 
      
      cv.ridgeB <- cv.glmnet(x=ZB, y=XB, alpha=0, nfolds=10, lambda=10^seq(1,-4,-.05))
      LambdaRidgeA[kk,ii,jj] <- cv.ridgeB$lambda.min
      
      #### Ridge Fit
      RidgeFitB <- predict(cv.ridgeA, newx=ZB, s=cv.ridgeA$lambda.min)
      RidgeFitA <- predict(cv.ridgeB, newx=ZA, s=cv.ridgeB$lambda.min)
      
      #### WARNING: POSSIBLE MULTICOLLINEARITY
      #### This implies the Cross Validation Procesure chooses a very large lambda and 
      #### the Ridge procedure shrinks all coefficients to 0!!!
      if (sum(RidgeFitA==RidgeFitA[1])==length(RidgeFitA)) {
        cat("MULTICOLINEARITY FOUND!\n")
        assign(paste("ILLED_B_",kk), list(coef(cv.ridgeB, s = "lambda.min"), cv.ridgeB$lambda))
        RidgeFitA <- rnorm(length(RidgeFitA))}
      if (sum(RidgeFitB==RidgeFitB[1])==length(RidgeFitB)) {
        cat("MULTICOLINEARITY FOUND!\n")
        assign(paste("ILLED_A_",kk), list(coef(cv.ridgeA, s = "lambda.min"), cv.ridgeA$lambda))
        RidgeFitB <- rnorm(length(RidgeFitB))}
      
      #### Augmented Instrument
      ZB <- cbind(ZB, RidgeFitB)
      ZA <- cbind(ZA, RidgeFitA)
      
      ZA <- ZA - rep(1,sss) %*% t(as.matrix(colMeans(ZA)))
      XA <- XA - mean(XA)
      YA <- YA - mean(YA)
      
      ZB <- ZB - rep(1,sss) %*% t(colMeans(ZB))
      XB <- XB - mean(XB)
      YB <- YB - mean(YB)
      
      ## Highest Correlation 2SLS
      ZL1A <- ZA[, which.max(cor(ZA,XA)), drop=FALSE]
      out_2sls <- tsls(YA, XA, NULL, ZL1A)
      Rb2sls1A[kk,ii,jj] <- out_2sls$b
      Rs2sls1A[kk,ii,jj] <- sqrt(out_2sls$VC1)
      
      ZL1B <- ZB[, which.max(cor(ZB,XB)), drop=FALSE]
      out_2sls <- tsls(YB, XB, NULL, ZL1B)
      Rb2sls1B[kk,ii,jj] <- out_2sls$b
      Rs2sls1B[kk,ii,jj] <- sqrt(out_2sls$VC1)
      
      ## Sample A: Lasso
      #### Iteratively Choose Penalty Loadings
      v0    <- XA
      pl    <- sqrt(t(crossprod(v0^2, ZA^2)))     # Initial Penalty Loadings
      b_tmp <- lasso_shooting_cpp(ZA, XA, Rlambda0 * pl)$b
      ind0  <- (abs(b_tmp) > 0)
      ZL0   <- as.matrix(ZA[, ind0, drop=FALSE])
      ind1  <- ind0
      if ( sum(ind0) != 0) {
        for (mm in 1:K) {
          v1    <- XA - ZL0 %*% lm(XA ~ -1 + ZL0)$coef
          pl    <- sqrt(crossprod(v1^2, ZA^2))   # Update Penalty Loadings
          b_tmp <- lasso_shooting_cpp(ZA, XA, Rlambda0 * pl)$b
          ind1  <- (abs(b_tmp) > 0)
          ZL0    <- as.matrix(ZA[, ind1, drop=FALSE])
        }
      }
      ZL1A    <- ZL0
      
      #### Post Lasso Estimation
      if ( sum(ind1) == 0 ) {
        RblassoCA[kk,ii,jj]  <- NaN
        RslassoCA[kk,ii,jj]  <- NaN
        RblassoCFA[kk,ii,jj] <- NaN
        RslassoCFA[kk,ii,jj] <- NaN
      } else {
        ### Run 2SLS with Selected Instruments
        bfs  <- lm(XA ~ -1 + ZL1A)$coef
        efs  <- XA - ZL1A %*% as.matrix(bfs)
        ### Asymptotic Variance Estimator
        Vfs  <- solve(crossprod(ZL1A)) %*% (t(ZL1A * (efs^2 %*% matrix(rep(1,sum(ind1)), nrow=1))) %*% ZL1A) %*% solve(crossprod(ZL1A))
        ### F-Statistic
        RFSA[kk,ii,jj] = (t(bfs) %*% solve(Vfs) %*% bfs) / sum(ind1)
        
        out_2sls <- tsls(YA, XA, NULL, ZL1A)
        out_full <- fuller(YA, XA, NULL, ZL1A)
        RblassoCA[kk,ii,jj]  <- out_2sls$b
        RslassoCA[kk,ii,jj]  <- sqrt(out_2sls$VC1)
        RblassoCFA[kk,ii,jj] <- out_full$b
        RslassoCFA[kk,ii,jj] <- sqrt(out_full$VC3)
        WA <- RFSA[kk,ii,jj] * sum(ind1)          # Wald Statistic
      }
      ## Sample B: Lasso
      #### Iteratively Choose Penalty Loadings
      v0    <- XB
      pl    <- sqrt(crossprod(v0^2, ZB^2))     # Initial Penalty Loadings
      b_tmp <- lasso_shooting_cpp(ZB, XB, Rlambda0 * pl)$b
      ind0  <- (abs(b_tmp) > 0)
      ZL0   <- as.matrix(ZB[, ind0, drop=FALSE])
      ind1  <- ind0
      if ( sum(ind0) != 0) {
        for (mm in 1:K) {
          v1    <- XB - ZL0 %*% lm(XB ~ -1 + ZL0)$coef
          pl    <- sqrt(crossprod(v1^2, ZB^2))   # Update Penalty Loadings
          b_tmp <- lasso_shooting_cpp(ZB, XB, Rlambda0 * pl)$b
          ind1  <- (abs(b_tmp) > 0)
          ZL0    <- as.matrix(ZB[, ind1, drop=FALSE])
        }
      }
      ZL1B    <- ZL0
      
      #### Post-Lasso Estimation
      if ( sum(ind1) == 0 ) {
        RblassoCB[kk,ii,jj]  <- NaN
        RslassoCB[kk,ii,jj]  <- NaN
        RblassoCFB[kk,ii,jj] <- NaN
        RslassoCFB[kk,ii,jj] <- NaN
      } else {
        ### Run 2SLS with Selected Instruments
        bfs  <- lm(XB ~ -1 + ZL1B)$coef
        efs  <- XB - ZL1B %*% as.matrix(bfs)
        ### Asymptotic Variance Estimator
        Vfs  <- solve(crossprod(ZL1B)) %*% (t(ZL1B * (efs^2 %*% matrix(rep(1,sum(ind1)), nrow=1))) %*% ZL1B) %*% solve(crossprod(ZL1B))
        ### F-Statistic
        RFSA[kk,ii,jj] = (t(bfs) %*% solve(Vfs) %*% bfs) / sum(ind1)
        
        out_2sls <- tsls(YB, XB, NULL, ZL1B)
        out_full <- fuller(YB, XB, NULL, ZL1B)
        RblassoCB[kk,ii,jj]  <- out_2sls$b
        RslassoCB[kk,ii,jj]  <- sqrt(out_2sls$VC1)
        RblassoCFB[kk,ii,jj] <- out_full$b
        RslassoCFB[kk,ii,jj] <- sqrt(out_full$VC3)
        WB <- RFSB[kk,ii,jj] *  sum(ind1)          # Wald Statistic
      }
      
      ## Combine Estimators
      ncolA <- dim(ZL1A)[2]
      ncolB <- dim(ZL1B)[2]
      if (is.vector(ZL1A)) {ncolA <- 1}
      if (is.vector(ZL1B)) {ncolB <- 1}
      
      if ( ncolA == 0 && ncolB == 0) {
        RblassoCC[kk,ii,jj]  <- NaN
        RslassoCC[kk,ii,jj]  <- NaN
        RblassoCFC[kk,ii,jj] <- NaN
        RslassoCFC[kk,ii,jj] <- NaN
        RFSC[kk,ii,jj] <-  0
      } else if (ncolA == 0 && ncolB != 0) {
        RblassoCC[kk,ii,jj]  <- RblassoCB[kk,ii,jj]
        RslassoCC[kk,ii,jj]  <- RslassoCB[kk,ii,jj]
        RblassoCFC[kk,ii,jj] <- RblassoCFB[kk,ii,jj]
        RslassoCFC[kk,ii,jj] <- RslassoCFB[kk,ii,jj]
        RFSC[kk,ii,jj] <- RFSB[kk,ii,jj]
      } else if (ncolA != 0 && ncolB == 0) {
        RblassoCC[kk,ii,jj]  <- RblassoCA[kk,ii,jj]
        RslassoCC[kk,ii,jj]  <- RslassoCA[kk,ii,jj]
        RblassoCFC[kk,ii,jj] <- RblassoCFA[kk,ii,jj]
        RslassoCFC[kk,ii,jj] <- RslassoCFA[kk,ii,jj]
        RFSC[kk,ii,jj] <- RFSA[kk,ii,jj]
      } else {
        weightA <- RslassoCB[kk,ii,jj]^2 / (RslassoCA[kk,ii,jj]^2 + RslassoCB[kk,ii,jj]^2)
        weightB <- RslassoCA[kk,ii,jj]^2 / (RslassoCA[kk,ii,jj]^2 + RslassoCB[kk,ii,jj]^2)
        RblassoCC[kk,ii,jj]  <- weightA * RblassoCA[kk,ii,jj] + weightB*RblassoCB[kk,ii,jj]
        RslassoCC[kk,ii,jj]  <- sqrt((weightA * RslassoCA[kk,ii,jj])^2 + (weightB * RslassoCB[kk,ii,jj])^2)
        weightA <- RslassoCFB[kk,ii,jj]^2/(RslassoCFA[kk,ii,jj]^2 + RslassoCFB[kk,ii,jj]^2)
        weightB <- RslassoCFA[kk,ii,jj]^2/(RslassoCFA[kk,ii,jj]^2 + RslassoCFB[kk,ii,jj]^2)
        RblassoCFC[kk,ii,jj] <- weightA  *RblassoCFA[kk,ii,jj] + weightB*RblassoCFB[kk,ii,jj]
        RslassoCFC[kk,ii,jj] <- sqrt((weightA * RslassoCFA[kk,ii,jj])^2 + (weightB * RslassoCFB[kk,ii,jj])^2)
        RFSC[kk,ii,jj] <- (WA+WB) / (ncolA + ncolB)
      }
      
      ## Sup-Score Test
      for (mm in 1:RnaA) {
        aEval <- RaTestA[mm]
        eTmp <- YA-aEval*XA
        ScoreVec <- crossprod(eTmp, ZA)
        ScoreStd <- sqrt(crossprod(eTmp^2, ZA^2))
        ScaledScore <- ScoreVec/(1.1*ScoreStd)
        RsupScore05A[kk,mm,ii,jj] <- max(abs(ScaledScore)) < RlambdaSS05
        RsupScoreA[kk,mm,ii,jj]   <- max(abs(ScaledScore))
        eTmp <- YB-aEval*XB
        ScoreVec <- crossprod(eTmp, ZB)
        ScoreStd <- sqrt(crossprod(eTmp^2, ZB^2))
        ScaledScore <- ScoreVec/(1.1*ScoreStd)
        RsupScore05B[kk,mm,ii,jj] <- max(abs(ScaledScore)) < RlambdaSS05
        RsupScoreB[kk,mm,ii,jj]   <- max(abs(ScaledScore))
      }
    } ## End - Simulation Loop
    ## Calculate Time Per Repetition
    cat(paste0("Time Elapsed:    ", round((proc.time()[3] - runtime[3]), digits=2), "s\n"))
    cat(paste0("Time/Repetition: ", round((proc.time()[3] - runtime[3])/nSims, digits=2), "s\n"))
    runtime <- proc.time()
  }
}

#### Produce Report
source("./result_save.R")
#### Write to XLSX
options("openxlsx.borderColour" = "#4F80BD")
All_Results <- list("[100,30]" = Results_100_30,  "[100,180]" = Results_100_180,
                    "[250,30]" = Results_250_30,  "[250,180]" = Results_250_180)
write.xlsx(All_Results, file = "./Sim_Output_Exp.xlsx",
           colNames=TRUE, rowNames = TRUE, borders = "rows")
