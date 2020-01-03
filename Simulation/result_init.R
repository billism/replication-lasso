#########################################
####### Matrice Storing Results #########
#########################################
## =====================================
## Conventional Estimators: 2SLS and Fuller 
## =====================================
b2sls <- array(0, dim=c(nSims,nn,ncp))
bfull <- array(0, dim=c(nSims,nn,ncp))
s2sls <- array(0, dim=c(nSims,nn,ncp))
sfull <- array(0, dim=c(nSims,nn,ncp))
## =====================================
## IV-Lasso
## =====================================
## -------------------------------------
## 2SLS: Using Single Instrument with the
## Highest Corr to the Endogenous Variable
## (In case of no instrument is selected)
## -------------------------------------
b2sls1 <- array(0, dim=c(nSims,nn,ncp))
s2sls1 <- array(0, dim=c(nSims,nn,ncp))
## -------------------------------------
## Penalty Loading: Iterative Algorithm
## -------------------------------------
FS <- array(0, dim=c(nSims,nn,ncp))
blassoC  <- array(0, dim=c(nSims,nn,ncp))
blassoCF <- array(0, dim=c(nSims,nn,ncp))
slassoC  <- array(0, dim=c(nSims,nn,ncp))
slassoCF <- array(0, dim=c(nSims,nn,ncp))
## -------------------------------------
## Select Penalty with Cross-Validation
## -------------------------------------
FSV <- array(0, dim=c(nSims,nn,ncp))
blassoCV <- array(0, dim=c(nSims,nn,ncp))
blassoCFV <- array(0, dim=c(nSims,nn,ncp))
slassoCV <- array(0, dim=c(nSims,nn,ncp))
slassoCFV <- array(0, dim=c(nSims,nn,ncp))
## -------------------------------------
## X-Dependent Penalty
## -------------------------------------
FSX <- array(0, dim=c(nSims,nn,ncp))
blassoCX <- array(0, dim=c(nSims,nn,ncp))
blassoCFX <- array(0, dim=c(nSims,nn,ncp))
slassoCX <- array(0, dim=c(nSims,nn,ncp))
slassoCFX <- array(0, dim=c(nSims,nn,ncp))
## -------------------------------------
## X-Independent Penalty
## -------------------------------------
FSn <- array(0, dim=c(nSims,nn,ncp))
blassoCn <- array(0, dim=c(nSims,nn,ncp))
blassoCFn <- array(0, dim=c(nSims,nn,ncp))
slassoCn <- array(0, dim=c(nSims,nn,ncp))
slassoCFn <- array(0, dim=c(nSims,nn,ncp))
## -------------------------------------
## Sup-Score Test
## -------------------------------------
aTest <- seq(from= .5, to=1.5, by=0.01);        
na    <- length(aTest);                         # 101
supScore   <- array(0, dim=c(nSims,na,nn,ncp))  # Sup-Score Test Statistic
supScore05 <- array(0, dim=c(nSims,na,nn,ncp))
lambdaSS05 <- qnorm(1 - .05/(2*p))              # Critical Value for Sup-Score Test
## =========================================
## Sample Split, Ridge, Lasso
## =====================================
## -------------------------------------
## CV - Ridge Penalty
## -------------------------------------
LambdaRidgeA <- array(0, dim=c(nSims,nn,ncp))
LambdaRidgeB <- array(0, dim=c(nSims,nn,ncp))
## -------------------------------------
## 2SLS: Using Single Instrument with the
## Highest Corr to the Endogenous Variable
## (In case no instrument is selected)
## -------------------------------------
Rb2sls1A <- array(0, dim=c(nSims,nn,ncp))
Rs2sls1A <- array(0, dim=c(nSims,nn,ncp))
Rb2sls1B <- array(0, dim=c(nSims,nn,ncp))
Rs2sls1B <- array(0, dim=c(nSims,nn,ncp))
## -------------------------------------
## IV-LASSO (Ridge)
## -------------------------------------
RFSA <- array(0, dim=c(nSims,nn,ncp))
RblassoCA  <- array(0, dim=c(nSims,nn,ncp))
RblassoCFA <- array(0, dim=c(nSims,nn,ncp))
RslassoCA  <- array(0, dim=c(nSims,nn,ncp))
RslassoCFA <- array(0, dim=c(nSims,nn,ncp))
RFSB <-  array(0, dim=c(nSims,nn,ncp))
RblassoCB  <- array(0, dim=c(nSims,nn,ncp))
RblassoCFB <- array(0, dim=c(nSims,nn,ncp))
RslassoCB  <- array(0, dim=c(nSims,nn,ncp))
RslassoCFB <- array(0, dim=c(nSims,nn,ncp))
RFSC <-  array(0, dim=c(nSims,nn,ncp))
RblassoCC  <- array(0, dim=c(nSims,nn,ncp))
RblassoCFC <- array(0, dim=c(nSims,nn,ncp))
RslassoCC  <- array(0, dim=c(nSims,nn,ncp))
RslassoCFC <- array(0, dim=c(nSims,nn,ncp))
## -------------------------------------
#### Sup-Score Test (Ridge)
## -------------------------------------
RaTestA <- seq(.5, 1.5, 0.01)
RnaA    <- length(RaTestA)
RsupScoreA   <- array(0, dim=c(nSims,RnaA,nn,ncp))
RsupScore05A <- array(0, dim=c(nSims,RnaA,nn,ncp))
RaTestB <- seq(0.5, 1.5, 0.01)
RnaB    <- length(RaTestB)
RsupScoreB   <- array(0, dim=c(nSims,RnaA, nn,ncp))
RsupScore05B <- array(0, dim=c(nSims,RnaA, nn,ncp))
RlambdaSS05  <- qnorm(1 - .05/(2*(p+1)))
## ========================================