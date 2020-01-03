#####################################################################
####                        CSHomePrice                         #####
#####################################################################
options(warn=-1)
## Clear Workspace
rm(list = ls(all.names = TRUE))
gc()
## Library
library(dummies)
library(MASS)
library(openxlsx)
source("./lambda_search.R")
source("./lasso_shooting.R")
source("./hetero_se.R")
## Parameters
c <- 1.1       # LASSO Penalty Multiplier

#####################################################################
## Prepare Data
#####################################################################
#### Import Data
#### ================================================================
URL1  <- "./data/temp_CSHomePrice.mat.csv"
URL2  <- "./data/probs.csv"

CSHPI <- read.csv(URL1)                   # CSHomePriceIndex
probs <- read.csv(URL2, sep='\t')         # Probability Controls

#### Inspect Data
#### ================================================================
## dim(CSHPI)
## head(names(CHSPI))
## sum(sapply(CSHPI, is.integer))
## sum(sapply(CSHPI, is.double))

## dim(probs)
## head(names(probs))
## sum(sapply(probs, is.integer))
## sum(sapply(probs, is.double))

#### Match by Circuit and Year
#### ================================================================
#runtime <- proc.time()
probs.new <- data.frame(matrix(nrow=nrow(CSHPI), ncol=ncol(probs)-2))
colnames(probs.new) <- colnames(probs)[-(1:2)]

for (i in 1:nrow(CSHPI)) {
  year          <- CSHPI$year[i]
  circuit       <- CSHPI$circuit[i]
  index         <- which((probs$Syear == year) & (probs$Scircuit == circuit))
  probs.new[i,] <- probs[index, -(1:2)]
}

#cat(proc.time() - runtime)
#runtime <- proc.time()
#prr.new <- do.call(rbind, lapply(1:nrow(CSHPI), function(i) return(probs[which((probs$Syear == CSHPI$year[i]) & (probs$Scircuit == CSHPI$circuit[i])), -(1:2)])))
#cat(proc.time() - runtime)

#####################################################################
## Variables
#####################################################################
#### Create Dummies
#### ================================================================
###### Generate dummy variables for year and circuit
###### For a varialbe with n levels, we only need n-1 dummies
yearn <- as.factor(CSHPI$year)
iyearn <- dummy(yearn)[,-1]
circuitn <- as.factor(CSHPI$circuit)
icircuitn <- dummy(circuitn)[,-1]

###### Intercept
col1 <- rep(1, nrow(iyearn))

###### Join year and circuit
jcircuityear <- icircuitn * CSHPI$year

######
recode <- function(x) {
  x.f <- as.factor(x)
  levels(x.f) <- as.character(1:length(levels(x.f)))
  return(as.integer(x.f))
}

######## recode "circuit" to 1-11 and recode "year" to 1-18.
######## recover from new index to "cirtuit" and "year":
######## - circuit: index %/% 18
######## - year   : index %%  18
cy <- ((recode(CSHPI$circuit) - 1) * max(recode(CSHPI$year)) + recode(CSHPI$year))
Dt <- data.matrix(dummy(cy))

#### Exogenous Variables to be partialed out (x)
#### ================================================================
###### - whether there were relevant cases in that circuit-year
###### - the number of takings appellate decisions
###### - circuit-specific effects
###### - time-specific effects
###### - circuit-specific time trends
x <- data.matrix(data.frame("missing_cy_12" = CSHPI$missing_cy_12,
                "numcasecat_12" = CSHPI$numcasecat_12,
                "icircuitn" = icircuitn,
                "iyearn" = iyearn,
                "jcircuityear" = jcircuityear,
                "col1" = col1))

###### Endogenous Variable (d)
###### ================================================================
d <- data.matrix(CSHPI$numpro_casecat_12)

###### Dependent Variable (y)
###### ================================================================
y  <- data.matrix(CSHPI$logpriceindex)

###### Candidate Instruments (z)
###### ================================================================
source("./instruments_cs.R")
z <- data.matrix(z)

###### Probability Controls (xP)
###### ================================================================
source("./probs.R")
xP <- data.matrix(xP)

###### Daniel Chen's 2 Instruments
###### ================================================================
IND_dan <- t(rep(0, ncol(z)))   # row vector
IND_dan[c(1,2)] <- 1            # corresponding to the first 2 columns of z
IND_dan <- as.logical(IND_dan)  # make it boolean

#####################################################################
## Additional Pre-Processing
#####################################################################
###### Remove Candidate Instruments with mean < 0.05
###### ================================================================
cat(paste0("There are ", dim(z)[2], " potential instrumental varialbes in z.\n"))
I  <- colMeans(z)  >= .05
cat(paste0("There are ", sum(I), " columns in z with a mean >= .05.\n"))
z  <- z[, I, drop=FALSE]
IND_dan <- IND_dan[I]

###### Remove Probability Controls with mean < 0.05
###### ================================================================
I  <- colMeans(xP) >= .05
cat(paste0("There are ", sum(I), " columns in xP with a mean >= .05.\n"))
xP <- xP[, I, drop=FALSE]
x <- cbind(x, xP)

###### Collapse to circuit-year level
###### ================================================================
######## Since all the controls, instruments, and the endogenous variable
######## vary only at the circuit-year level
x <- solve(crossprod(Dt)) %*% crossprod(Dt, x)
y <- solve(crossprod(Dt)) %*% crossprod(Dt, y)
d <- solve(crossprod(Dt)) %*% crossprod(Dt, d)
z <- solve(crossprod(Dt)) %*% crossprod(Dt, z)
###### Partial out x
###### ================================================================
xxinv <- solve(crossprod(x))
My <- y - x %*% xxinv %*% crossprod(x, y)
Md <- d - x %*% xxinv %*% crossprod(x, d)
Mz <- z - x %*% xxinv %*% crossprod(x, z)

###### Remove Instruments with small standard deviation
I  <- apply(Mz, 2, sd) > 1e-6
cat(paste0("There are ", sum(I), " columns in Mz with a std >= 1e-6.\n"))
Mz <- Mz[, I, drop=FALSE]        # 147 Left
IND_dan <- IND_dan[I]

#####################################################################
## Select Instruments with Lasso
#####################################################################
n <- nrow(Mz)
k <- ncol(Mz)
#### Standardize For Lasso Selection
std_Md <- sd(Md)
std_Mz <- t(replicate(n, apply(Mz, 2, sd)))
Md <- Md / std_Md     # Md is one dimensional
Mz <- Mz / std_Mz     

#### - Use Lasso to select a desired number of instruments
#### - Select using data-dependent penalty "A"
#### - Select using data-dependent penalty "L"
#### ================================================================
#### Using Lasso to select 2 instruments
IND <- lambda_search(Mz, Md, 2, 4000)$Index
cat(paste0(sum(IND), " instruments selected!\n"))
if (sum(IND) != 0){
  for(rname in rownames(IND)[IND]) {
    cat(paste0(">>> ", rname, "\n"))
  }
}
zt  <- Mz[, IND, drop=FALSE]  # selected instruments
# regress Md on z1 and compute residual v
v   <- Md - zt %*% solve(crossprod(zt)) %*% crossprod(zt, Md)

#### ================================================================
#### Lasso with data-dependent penalty "A"
#### Refer to Belloni et al (2012) Appendix A
St     <- Mz * (v %*% t(rep(1,k)))
Ups    <- sqrt(colSums(St^2)/n)
lambda <- 2 * c * sqrt(2 * n * (log(2 * k)))
###### Note: apply penalty loadings to Mz or to lambda will lead to different results => scale it back
PIpA   <- lasso_shooting(Mz / (matrix(rep(1, n), nrow=n) %*% Ups), Md, lambda)$b
INDpA  <- abs(PIpA) > 1.0e-4
cat(paste0(sum(INDpA), " instruments selected!\n"))
if (sum(INDpA) != 0){
  for(rname in rownames(INDpA)[INDpA]) {
    cat(paste0(">>> ", rname, "\n"))
  }
}
zt     <- Mz[, INDpA, drop=FALSE]
v      <- Md - zt %*% solve(crossprod(zt)) %*% crossprod(zt, Md)

#### ================================================================
#### Lasso with data-dependent penalty "L"
St    <- Mz * (v %*% t(rep(1,k)))
UpsR  <- sqrt(colSums(St^2)/n)    # "R" for "Refined"
PIpL  <- lasso_shooting(Mz / (matrix(rep(1, n), nrow=n) %*% UpsR), Md, lambda)$b
INDpL <- abs(PIpL) > 1.0e-4
cat(paste0(sum(INDpL), " instruments selected!\n"))
if (sum(INDpL) != 0){
  for(rname in rownames(INDpL)[INDpL]) {
    cat(paste0(">>> ", rname, "\n"))
  }
}


#####################################################################
## OLS Estimation and 2SLS Estimation
#####################################################################
#### Scale back to original units
Md <- Md * std_Md
Mz <- Mz * std_Mz
#### Number of nonredundant partialed out x's
kx <- ncol(x)

#### OLS Estimate
#### ================================================================
OLS <- solve(crossprod(Md)) %*% crossprod(Md, My)
et  <- My - Md %*% OLS
sOLS <- sqrt((n-1)/(n-kx-1) * hetero_se(Md, et, solve(crossprod(Md)))$se)

#### 2SLS with Daniel's Instruments
#### ================================================================
zt  <- Mz[,IND_dan, drop=FALSE]
FS_dan  <- solve(crossprod(zt)) %*% crossprod(zt, Md)
ut  <- Md - zt %*% FS_dan
tmp <- hetero_se(zt, ut, solve(crossprod(zt)))
sFS_dan <- sqrt((n-ncol(zt)-kx)/(n-ncol(zt))) %*% tmp$se
VFS_dan <- ((n-ncol(zt)-kx)/(n-ncol(zt))) * tmp$vhetero
FS_F_dan <- t(FS_dan) %*% solve(VFS_dan) %*% FS_dan      # F-Statistic
Fs_p_dan <- 1 - pchisq(FS_F_dan, 2)                      # P-Value

TSLS_dan <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% (crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,My))
et <- My - Md %*% TSLS_dan
xxinv   <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% (crossprod(Md,zt) %*% solve(crossprod(zt)))
sTSLS_dan <- sqrt((n-1)/(n-kx-1)) * hetero_se(zt, et, xxinv)$se

#### 2SLS with 2 Instruments selected by Lasso
#### ================================================================
zt  <- Mz[,IND, drop=FALSE]
FS_2  <- solve(crossprod(zt)) %*% crossprod(zt, Md)
ut  <- Md - zt %*% FS_2
tmp <- hetero_se(zt, ut, solve(crossprod(zt)))
sFS_2 <- sqrt((n-ncol(zt)-kx)/(n-ncol(zt))) %*% tmp$se
VFS_2 <- ((n-ncol(zt)-kx)/(n-ncol(zt))) * tmp$vhetero
FS_F_2 <- t(FS_2) %*% solve(VFS_2) %*% FS_2      # F-Statistic
Fs_p_2 <- 1 - pchisq(FS_F_2, 2)                      # P-Value

TSLS_2 <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% (crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,My))
et <- My - Md %*% TSLS_2
xxinv   <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% (crossprod(Md,zt) %*% solve(crossprod(zt)))
sTSLS_2 <- sqrt((n-1)/(n-kx-1)) * hetero_se(zt, et, xxinv)$se

#### 2SLS with Instruments selected by Lasso (Penalty "A")
#### ================================================================
zt  <- Mz[,INDpA, drop=FALSE]
FS_A  <- solve(crossprod(zt)) %*% crossprod(zt, Md)
ut  <- Md - zt %*% FS_A
tmp <- hetero_se(zt, ut, solve(crossprod(zt)))
sFS_A <- sqrt((n-ncol(zt)-kx)/(n-ncol(zt))) %*% tmp$se
VFS_A <- ((n-ncol(zt)-kx)/(n-ncol(zt))) * tmp$vhetero
FS_F_A <- t(FS_A) %*% solve(VFS_A) %*% FS_A      # F-Statistic
Fs_p_A <- 1 - pchisq(FS_F_A, ncol(zt))           # P-Value

TSLS_A <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% (crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,My))
et <- My - Md %*% TSLS_A
xxinv   <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% (crossprod(Md,zt) %*% solve(crossprod(zt)))
sTSLS_A <- sqrt((n-1)/(n-kx-1)) * hetero_se(zt, et, xxinv)$se

#### 2SLS with Instruments selected by Lasso (Penalty "L")
#### ================================================================
zt  <- Mz[,INDpL, drop=FALSE]
FS_L  <- solve(crossprod(zt)) %*% crossprod(zt, Md)
ut  <- Md - zt %*% FS_L
tmp <- hetero_se(zt, ut, solve(crossprod(zt)))
sFS_L <- sqrt((n-ncol(zt)-kx)/(n-ncol(zt))) %*% tmp$se
VFS_L <- ((n-ncol(zt)-kx)/(n-ncol(zt))) * tmp$vhetero
FS_F_L <- t(FS_L) %*% solve(VFS_L) %*% FS_L      # F-Statistic
Fs_p_L <- 1 - pchisq(FS_F_L, ncol(zt))           # P-Value

TSLS_L <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% (crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,My))
et <- My - Md %*% TSLS_L
xxinv   <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% (crossprod(Md,zt) %*% solve(crossprod(zt)))
sTSLS_L <- sqrt((n-1)/(n-kx-1)) * hetero_se(zt, et, xxinv)$se

#### Specification Test
#### ================================================================
zt_dan <- Mz[,IND_dan, drop=FALSE]
et_dan <- My - Md %*% TSLS_dan
score_dan <- zt_dan * (et_dan %*% t(rep(1,2)))
score_L   <- zt * (et %*% t(rep(1, ncol(zt))))
se_dif    <- sqrt(sTSLS_dan^2 + sTSLS_L^2 - 2 * solve(crossprod(Md,zt_dan) %*%
                      solve(crossprod(zt_dan)) %*% crossprod(zt_dan,Md)) %*% 
                     (crossprod(Md,zt_dan) %*% solve(crossprod(zt_dan))) %*% 
                      crossprod(score_dan, score_L) %*% solve(crossprod(zt)) %*% 
                      crossprod(zt,Md) %*% solve(crossprod(Md,zt) %*% 
                      solve(crossprod(zt)) %*% crossprod(zt,Md)))
t_dif     <- (TSLS_dan - TSLS_L) / se_dif

#### 2SLS with Instruments selected by Lasso and Daneil (Union)
#### ================================================================
INDunion <- IND_dan | INDpL
zt  <- Mz[,INDunion, drop=FALSE]
FS_U  <- ginv(crossprod(zt)) %*% crossprod(zt, Md)
ut  <- Md - zt %*% FS_U
tmp <- hetero_se(zt, ut, solve(crossprod(zt)))
sFS_U <- sqrt((n-ncol(zt)-kx)/(n-ncol(zt))) %*% tmp$se
VFS_U <- ((n-ncol(zt)-kx)/(n-ncol(zt))) * tmp$vhetero
FS_F_U <- t(FS_U) %*% solve(VFS_U) %*% FS_U      # F-Statistic
Fs_p_U <- 1 - pchisq(FS_F_U, ncol(zt))           # P-Value

TSLS_U <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*% 
                (crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,My))
et <- My - Md %*% TSLS_U
xxinv   <- solve(crossprod(Md,zt) %*% solve(crossprod(zt)) %*% crossprod(zt,Md)) %*%
                (crossprod(Md,zt) %*% solve(crossprod(zt)))
sTSLS_U <- sqrt((n-1)/(n-kx-1)) * hetero_se(zt, et, xxinv)$se


#####################################################################
## Sup-Score Test
#####################################################################
aVec <- matrix(seq(-0.5, 1.5, 0.001), ncol=1)
SupScore <- matrix(rep(0, length(aVec)), nrow=length(aVec))
for (jj in 1:length(aVec)) {
  aT <- aVec[jj,1]
  eTmp <- My - Md %*% aT
  ScoreVect <- crossprod(eTmp, Mz)
  ScoreStd  <- sqrt(crossprod(eTmp^2, Mz^2))
  ScaledScore <- ScoreVect / (1.1 * ScoreStd)
  SupScore[jj, 1] <- max(abs(ScaledScore))
}
ind <- SupScore <= qnorm(1 - 0.05 / (2 * ncol(Mz)))
SupScoreInt <- cbind(min(aVec[ind]), max(aVec[ind]))

#####################################################################
## Save Results
#####################################################################
result_cs <- c(OLS,      sOLS,
            TSLS_dan, sTSLS_dan, FS_F_dan,
            TSLS_2,   sTSLS_2,   FS_F_2,
            TSLS_L,   sTSLS_L,   FS_F_L, sum(INDpL),
            TSLS_U,   sTSLS_U,    FS_F_U, sum(INDunion),
            t_dif)
result_cs <- matrix(result_cs, ncol=1)
rownames(result_cs) <- c("OLS", "s.e.",
                      "2SLS (Daniel)",  "s.e.", "FS_W",
                      "2SLS (2)",       "s.e.", "FS_W",
                      "Post-Lasso",     "s.e.", "FS_W", "S",
                      "Post-Lasso+",    "s.e.", "FS_W", "S",
                      "Spec. Test")
colnames(result_cs) <- c("Case-Shiller")

options("openxlsx.borderColour" = "#4F80BD")
write.xlsx(result_cs, file = "./Result_CSHPI.xlsx",
           colNames=TRUE, rowNames = TRUE, borders = "rows")




