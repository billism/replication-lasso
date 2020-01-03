#####################################################################
## Results
#####################################################################
Results <- matrix(NA, ncol=4, nrow=12)
colnames(Results) <- c("N(0)", "Bias", "MAD", "rp(0.05)")
rownames(Results) <- c("2SLS(100)", "FULL(100)",
                       "Post-LASSO (Iterative)", "Post-LASSO-F (Iterative)",
                       "Post-LASSO (CV)", "Post-LASSO-F (CV)",
                       "Post-LASSO (X-indep.)", "Post-LASSO-F (X-indep.)",
                       "Post-LASSO (X-dep.)", "Post-LASSO-F (X-dep.)",
                       "Post-Lasso (Ridge)",  "Post-Lasso-F (Ridge)")

for (ii in 1:nn) {
  for (jj in 1:ncp) {
    ## (1) 2SLS
    sfull[is.nan(sfull)] <- 10000
    Results[1,1] <- sum(is.nan(b2sls[,ii,jj]))
    Results[1,2] <- mean(b2sls[,ii,jj])-1
    Results[1,3] <- mean(abs(b2sls[,ii,jj]-1))
    Results[1,4] <- (nSims - sum(abs((b2sls[,ii,jj]-1)/s2sls[,ii,jj]) < qnorm(0.95)))/nSims
    ## (2)Fuller
    Results[2,1] <- sum(is.nan(bfull[,ii,jj]))
    Results[2,2] <- mean(bfull[,ii,jj])-1
    Results[2,3] <- mean(abs(bfull[,ii,jj]-1))
    Results[2,4] <- (nSims - sum(abs((bfull[,ii,jj]-1)/sfull[,ii,jj]) < qnorm(0.95)))/nSims
    ## (3) IV-Lasso (Iterative) 2SLS
    Results[3,1] <- sum(is.nan(blassoC[,ii,jj]))
    blassoC_aux  <- blassoC
    blassoC_aux[is.nan(blassoC_aux)] <- b2sls1[is.nan(blassoC_aux)]
    Results[3,2] <- mean(blassoC_aux[,ii,jj])-1
    Results[3,3] <- mean(abs(blassoC_aux[,ii,jj]-1))
    element_chosen <- which(!is.nan(blassoC[,ii,jj]))
    test1 <- sum(abs((blassoC[element_chosen,ii,jj]-1)/slassoC[element_chosen,ii,jj]) < qnorm(0.95))
    test2 <- sum(supScore05[-element_chosen,51,ii,jj])
    Results[3,4] <- (nSims - test1 - test2)/nSims
    ## (4) IV-Lasso (Iterative) Fuller
    Results[4,1] <- sum(is.nan(blassoCF[,ii,jj]))
    blassoCF_aux <- blassoCF
    blassoCF_aux[is.nan(blassoCF_aux)] <- b2sls1[is.nan(blassoCF_aux)]
    Results[4,2] <- mean(blassoCF_aux[,ii,jj])-1
    Results[4,3] <- mean(abs(blassoCF_aux[,ii,jj]-1))
    element_chosen <- which(!is.nan(blassoCF[,ii,jj]))
    test1 <- sum(abs((blassoCF[element_chosen,ii,jj]-1)/slassoCF[element_chosen,ii,jj]) < qnorm(0.95))
    test2 <- sum(supScore05[-element_chosen,51,ii,jj])
    Results[4,4] <- (nSims - test1 - test2)/nSims
    ## (5) IV-Lasso (CV) 2SLS
    # Results[5,1] <- sum(is.nan(blassoCV[,ii,jj]))
    # blassoCV_aux <- blassoCV
    # blassoCV_aux[is.nan(blassoCV_aux)] <- b2sls1[is.nan(blassoCV_aux)]
    # Results[5,2] <- mean(blassoCV_aux[,ii,jj])-1
    # Results[5,3] <- mean(abs(blassoCV_aux[,ii,jj]-1))
    # element_chosen <- which(!is.nan(blassoCV[,ii,jj]))
    # test1 <- sum(abs((blassoCV[element_chosen,ii,jj]-1)/slassoCV[element_chosen,ii,jj]) < qnorm(0.95))
    # test2 <- sum(supScore05[-element_chosen,51,ii,jj])
    # Results[5,4] <- (nSims - test1 - test2)/nSims
    # ## (6) IV-Lasso (CV) Fuller
    # Results[6,1] <- sum(is.nan(blassoCFV[,ii,jj]))
    # blassoCFV_aux <- blassoCFV
    # blassoCFV_aux[is.nan(blassoCFV_aux)] <- b2sls1[is.nan(blassoCFV_aux)]
    # Results[6,2] <- mean(blassoCFV_aux[,ii,jj])-1
    # Results[6,3] <- mean(abs(blassoCFV_aux[,ii,jj]-1))
    # element_chosen <- which(!is.nan(blassoCFV[,ii,jj]))
    # test1 <- sum(abs((blassoCFV[element_chosen,ii,jj]-1)/slassoCFV[element_chosen,ii,jj]) < qnorm(0.95))
    # test2 <- sum(supScore05[-element_chosen,51,ii,jj])
    # Results[6,4] <- (nSims - test1 - test2)/nSims
    ## (7) IV-Lasso (X-Independent) 2SLS
    Results[7,1] <- sum(is.nan(blassoCn[,ii,jj]))
    blassoCn_aux <- blassoCn
    blassoCn_aux[is.nan(blassoCn_aux)] <- b2sls1[is.nan(blassoCn_aux)]
    Results[7,2] <- mean(blassoCn_aux[,ii,jj])-1
    Results[7,3] <- mean(abs(blassoCn_aux[,ii,jj]-1))
    element_chosen <- which(!is.nan(blassoCn[,ii,jj]))
    test1 <- sum(abs((blassoCn[element_chosen,ii,jj]-1)/slassoCn[element_chosen,ii,jj]) < qnorm(0.95))
    test2 <- sum(supScore05[-element_chosen,51,ii,jj])
    Results[7,4] <- (nSims - test1 - test2)/nSims
    ## (8) IV-Lasso (X-Independent) Fuller
    Results[8,1]  <- sum(is.nan(blassoCFn[,ii,jj]))
    blassoCFn_aux <- blassoCFn
    blassoCFn_aux[is.nan(blassoCFn_aux)] <- b2sls1[is.nan(blassoCFn_aux)]
    Results[8,2] <- mean(blassoCFn_aux[,ii,jj])-1
    Results[8,3] <- mean(abs(blassoCFn_aux[,ii,jj]-1))
    element_chosen <- which(!is.nan(blassoCFn[,ii,jj]))
    test1 <- sum(abs((blassoCFn[element_chosen,ii,jj]-1)/slassoCFn[element_chosen,ii,jj]) < qnorm(0.95))
    test2 <- sum(supScore05[-element_chosen,51,ii,jj])
    Results[8,4] <- (nSims - test1 - test2)/nSims
    ## (9) IV-Lasso (X-Dependent) 2SLS
    Results[9,1] <- sum(is.nan(blassoCX[,ii,jj]))
    blassoCX_aux <- blassoCX
    blassoCX_aux[is.nan(blassoCX_aux)] <- b2sls1[is.nan(blassoCX_aux)]
    Results[9,2] <- mean(blassoCX_aux[,ii,jj])-1
    Results[9,3] <- mean(abs(blassoCX_aux[,ii,jj]-1))
    element_chosen <- which(!is.nan(blassoCX[,ii,jj]))
    test1 <- sum(abs((blassoCX[element_chosen,ii,jj]-1)/slassoCX[element_chosen,ii,jj]) < qnorm(0.95))
    test2 <- sum(supScore05[-element_chosen,51,ii,jj])
    Results[9,4] <- (nSims - test1 - test2)/nSims
    ## (10) IV-Lasso (X-Dependent) Fuller
    Results[10,1] <- sum(is.nan(blassoCFX[,ii,jj]))
    blassoCFX_aux <- blassoCFX
    blassoCFX_aux[is.nan(blassoCFX_aux)] <- b2sls1[is.nan(blassoCFX_aux)]
    Results[10,2] <- mean(blassoCFX_aux[,ii,jj])-1
    Results[10,3] <- mean(abs(blassoCFX_aux[,ii,jj]-1))
    element_chosen <- which(!is.nan(blassoCFX[,ii,jj]))
    test1 <- sum(abs((blassoCFX[element_chosen,ii,jj]-1)/slassoCFX[element_chosen,ii,jj]) < qnorm(0.95))
    test2 <- sum(supScore05[-element_chosen,51,ii,jj])
    Results[10,4] <- (nSims - test1 - test2)/nSims
    ## (11) IV-Lasso (Ridge) 2SLS
    Results[11,1] <- sum(is.nan(RblassoCC[,ii,jj]))
    RblassoCC_aux <- RblassoCC
    RblassoCC_aux[is.nan(RblassoCC_aux)] <- b2sls1[is.nan(RblassoCC_aux)]
    Results[11,2] <- mean(RblassoCC_aux[,ii,jj])-1
    Results[11,3] <- mean(abs(RblassoCC_aux[,ii,jj]-1))
    element_chosen <- which(!is.nan(RblassoCC[,ii,jj]))
    test1 <- sum(abs((RblassoCC[element_chosen,ii,jj]-1)/RslassoCC[element_chosen,ii,jj]) < qnorm(0.95))
    test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,ii,jj])+sum(RsupScore05B[-element_chosen,51,ii,jj]))
    Results[11,4] <- (nSims - test1 - test2)/nSims
    ## (12) IV-Lasso (Ridge) Fuller
    Results[12,1]  <- sum(is.nan(RblassoCFC[,ii,jj]))
    RblassoCFC_aux <- RblassoCFC
    RblassoCFC_aux[is.nan(RblassoCFC_aux)] <- b2sls1[is.nan(RblassoCFC_aux)]
    Results[12,2]  <- mean(RblassoCFC_aux[,ii,jj])-1
    Results[12,3]  <- mean(abs(RblassoCFC_aux[,ii,jj]-1))
    element_chosen <- which(!is.nan(RblassoCFC[,ii,jj]))
    test1 <- sum(abs((RblassoCFC[element_chosen,ii,jj]-1)/RslassoCFC[element_chosen,ii,jj]) < qnorm(0.95))
    test2 <- 1/2*(sum(RsupScore05A[-element_chosen,51,ii,jj])+sum(RsupScore05B[-element_chosen,51,ii,jj]))
    Results[12,4] <- (nSims - test1 - test2)/nSims

    
    assign(paste0("Results_",ns[ii],"_",cps[jj]), Results)
  }
}
