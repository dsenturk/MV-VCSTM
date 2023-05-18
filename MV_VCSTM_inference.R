MV_VCSTM_inference <- function(CARout # Output from function MV_VCSTM_CAR including the following: 
                                 # xiEst: estimates of region-specific PC scores (matrix of dimension nregion*(M1+M2))
                                 # alphaEst: estimates of spatial variance parameters (vector of length (M1+M2))
                                 # beta1Est: estimate of VCFs in the 1st dimension (matrix of dimension ngrid * 3)
                                 # beta2Est: estimate of VCFs in the 2nd dimension (matrix of dimension ngrid * 3)
                                 # beta1.mc: all posterior samples of VCFs in the 1st dimension (array of dimension (1000*ngrid*3))
                                 # beta2.mc: all posterior samples of VCFs in the 2nd dimension (array of dimension (1000*ngrid*3))
                                 # sigma1Est: estimate of measurement error variance in the 1st dimension (scalar)
                                 # sigma2Est: estimate of measurement error variance in the 2nd dimension (scalar)
                                 # nuEst: estimate of spatial smoothing parameter (scalar)
                                
){
  
#############################################################################
## Description: Function for obtaining prediction and inference for multivariate trajectories (Estimation Algorithm steps 5 in Section 2.2) 
## Definition:  nregion: number of regions, 
##              ngrid: number of follow-up time points.
## Args: see above
## Returns:     list()
##              Y1.hat: Region-specific trajectory predictions in the 1st dimension (matrix of dimension nregion*ngrid)
##              Y2.hat: Region-specific trajectory predictions in the 2nd dimension (matrix of dimension nregion*ngrid)
##              beta1.CB: Lower and upper bound of 95% simultaneous confidence band for all VCFs in 1st dimension (matrix of dimension ngrid*6)
##              beta2.CB: Lower and upper bound of 95% simultaneous confidence band for all VCFs in 2nd dimension (matrix of dimension ngrid*6)
############################################################################# 

# Install missing packages
list.of.packages <- c("refund", "fda", "mgcv", "MASS", "caTools", "locpol", 
                        "KernSmooth", "fANCOVA", "mgcv", "mvtnorm", "spdep", 
                        "fields", "R2WinBUGS", "pbugs")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 
  
# Load packages
library(refund)
library(fda)
library(mgcv)
library(MASS)
library(caTools)
library(locpol)
library(KernSmooth)
library(fANCOVA)
library(mvtnorm)
library(spdep)
library(fields)
library(Matrix)
library(nimble)
library(tidyverse)
library(fields)
library(R2WinBUGS)
library(pbugs)
  
## Define time points, number of regions, time points and eigen components in each dimension
nregion <- data.T$nregion
ngrid <- data.T$ngrid
gridPoints <- data.T$gridPoints


# Final number of eigencomponents
L <- dim(CARout$xiEst)[2]

# Estimated VCFs
# 1st dimension
beta01 <- CARout$beta1Est[,1]
beta11 <- CARout$beta1Est[,2]
beta21 <- CARout$beta1Est[,3]

# 2nd dimension
beta02 <- CARout$beta2Est[,1]
beta12 <- CARout$beta2Est[,2]
beta22 <- CARout$beta2Est[,3]

# All MCMC samples of estimated VCFs
beta01.mc <- CARout$beta1.mc[,,1]
beta11.mc <- CARout$beta1.mc[,,2]
beta21.mc <- CARout$beta1.mc[,,3]
beta02.mc <- CARout$beta2.mc[,,1]
beta12.mc <- CARout$beta2.mc[,,2]
beta22.mc <- CARout$beta2.mc[,,3]

# Get SD of VCFs from all MCMC samples
SD_beta01 <- apply(beta01.mc, 2, sd) 
SD_beta11 <- apply(beta11.mc, 2, sd)
SD_beta21 <- apply(beta21.mc, 2, sd)
SD_beta02 <- apply(beta02.mc, 2, sd)
SD_beta12 <- apply(beta12.mc, 2, sd)
SD_beta22 <- apply(beta22.mc, 2, sd)

## Reconstruct multivariate region-specific trajectories
Y1.hat <- Y2.hat <- matrix(rep(NA, nregion*ngrid), nrow = nregion)

for(i in 1:nregion){
  Y1.hat[i,] <- beta01 + Xmat[i, 1:2] %*% t(cbind(beta11, beta21)[,1:2]) + 
                 CARout$xiEst[i,(1:L)] %*% t(MCARout$psi11Est[,(1:L)])
  Y2.hat[i,] <- beta02 + Xmat[i, 1:2] %*% t(cbind(beta12, beta22)[,1:2]) + 
                 CARout$xiEst[i,(1:L)] %*% t(MCARout$psi22Est[,(1:L)])
}



# 95% simultanous confidence band
# 1st dimension
# beta01
cb <- quantile(apply(abs(beta01.mc - beta01)/(rep(1,1000) %*% t(SD_beta01)),1,max), .95)
beta01.CB.lower <- beta01 - cb*SD_beta01
beta01.CB.upper <- beta01 + cb*SD_beta01


# beta11
cb <- quantile(apply(abs(beta11.mc - beta11)/(rep(1,1000) %*% t(SD_beta11)),1,max), .95)
beta11.CB.lower <- beta11 - cb*SD_beta11
beta11.CB.upper <- beta11 + cb*SD_beta11


# beta21
cb <- quantile(apply(abs(beta21.mc - beta21)/(rep(1,1000) %*% t(SD_beta21)),1,max), .95)
beta21.CB.lower <- beta21 - cb*SD_beta21
beta21.CB.upper <- beta21 + cb*SD_beta21

#
beta1.CB <- cbind(beta01.CB.lower, beta01.CB.upper, beta11.CB.lower, beta11.CB.upper,
                 beta21.CB.lower, beta21.CB.upper)

# 2nd dimension
# beta02
cb <- quantile(apply(abs(beta02.mc - beta02)/(rep(1,1000) %*% t(SD_beta02)),1,max), .95)
beta02.CB.lower <- beta02 - cb*SD_beta02
beta02.CB.upper <- beta02 + cb*SD_beta02


# beta12
cb <- quantile(apply(abs(beta12.mc - beta12)/(rep(1,1000) %*% t(SD_beta12)),1,max), .95)
beta12.CB.lower <- beta12 - cb*SD_beta12
beta12.CB.upper <- beta12 + cb*SD_beta12


# beta22
cb <- quantile(apply(abs(beta22.mc - beta22)/(rep(1,1000) %*% t(SD_beta22)),1,max), .95)
beta22.CB.lower <- beta22 - cb*SD_beta22
beta22.CB.upper <- beta22 + cb*SD_beta22

#
beta2.CB <- cbind(beta02.CB.lower, beta02.CB.upper, beta12.CB.lower, beta12.CB.upper,
                  beta22.CB.lower, beta22.CB.upper)

## Construct output
out <- list(Y1.hat = Y1.hat, Y2.hat = Y2.hat, 
            beta1.CB = beta1.CB, beta2.CB = beta2.CB)
return(out)
}
  
