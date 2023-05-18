MV_VCSTM_univariate_decomposition <- function(data,  # data.frame with columns c("rid", "y1", "y2", "t", "dt", "x1", "x2"), for further analysis using MV-VCSTM
                                                         #         DATA.FRAME COLUMNS (data is stored in long format): 
                                                         #             rid: region IDs (vector of length ngrid*nregion) ngrid: number of follow-up time points; nregion: number of regions
                                                         #             y1: outcome is the 1st dimension (vector of length ngrid*nregion)
                                                         #             y2: outcome is the 2nd dimension (vector of length ngrid*nregion)
                                                         #             t: follow-up time (vector of length ngrid*nregion) 
                                                         #             dt: order of follow-up time in a region (vector of length ngrid*nregion)   
                                                         #             x1: 1st covariate (vector of length ngrid*nregion)
                                                         #             x2: 2nd covariate (vector of length ngrid*nregion)
                                              FVEvalue = 0.99   # A scalar between 0 to 1 to choose number of eigencomponents in each dimension. 
                                                                # To retain enough information at the initial step of the algorithm, it is better 
                                                                # to pick a large value for FVEvalue (close to 1). In our implementation, we use FVEvalue = 0.99. 
){
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

#############################################################################
## Description: Function for fitting multivariate varying coefficient function (Estimation Algorithm steps 1 in Section 2.2),
##              as well as univariate FPCA decomposition (Estimation Algorithm steps 2) of MV-VCSTM, 
##              including estimation of univariate eigenfunctions and eigenscores in each dimension. 
## Definition:  nregion: number of regions, 
##              ngrid: number of follow-up time points,
##              M1: estimated number of eigen components in the first dimension,
##              M2: estimated number of eigen components in the second dimension,
##              gridPoints: follow-up time points.
## Args:        see above
## Returns:     list()
##              dat.s1: data.frame in the first dimension with columns c("id", "gridPoints", "Y1", "dt"), for further analysis using MST-FPCA
#                 DATA.FRAME COLUMNS (data is stored in long format): 
#                  id: region IDs (vector of length ngrid*nregion) 
#                  gridPoints: follow-up time points (vector of length ngrid*nregion) 
#                  Y1: observed hospitalization rate data (vector of length ngrid*nregion) 
#                  dt: order of follow-up time points (vector of length ngrid*nregion) 
##              dat.s2: data.frame in the second dimension with columns c("id", "gridPoints", "Y2", "dt"), for further analysis using MST-FPCA
#                 DATA.FRAME COLUMNS (data is stored in long format) is the same as dat.s1. 
##              sig1Est: estimated error variance in the first-dimension (scalar)
##              sig2Est: estimated error variance in the second-dimension (scalar)
##              efun11: estimated univariate eigenfunctions in the first dimension (matrix of dimension ngrid*M1) 
##              efun22: estimated univaraite eigenfunctions in the second dimension (matrix of dimension ngrid*M2) 
##              e.scores11: estimated univariate PC eigenscores in the first dimension (matrix of dimension nregion*M1)
##              e.scores22: estimated univariate PC eigenscores in the second dimension (matrix of dimension nregion*M2)
##              M1: estimated number of eigencomponents chosen by FVE in the first diemension (scalar)
##              M2: estimated number of eigencomponents chosen by FVE in the second diemension (scalar)
#############################################################################


# Format data
df <- data
# Number of regions n
nregion <- length(unique(df$rid))

# Number of time points T
ngrid <- length(unique(df$t))
gridPoints <- unique(df$t)

############################################
# Step 1: Fit a multivariate varying coefficient model 
# to data under the working independence assumption
############################################
model <- gam(list(y1 ~ s(t, m = 2) + s(t, m = 2, by = x1) + s(t, m = 2, by = x2),  
                  y2 ~ s(t, m = 2) + s(t, m = 2, by = x1) + s(t, m = 2, by = x2)),
             data = df, method = "REML", family = mvn(d = 2))

# Obtain the residuals
fitted.v1 <- model$fitted.values[,1]
fitted.v2 <- model$fitted.values[,2]
X1.c <- df$y1 - fitted.v1
X2.c <- df$y2 - fitted.v2 



###########################################
## Step 2: Univariate FPCA in each dimension
###########################################

## FPCA in the 1st dimension 
dat.s1 <- df[,c("rid", "t", "dt", "y1")] 
colnames(dat.s1) <- c("id", "gridPoints", "dt", "Y1")


# Calculate the raw covariances
x0 <- rep(gridPoints, each = ngrid)
x1 <- rep(gridPoints, ngrid)
Xmat1 <- matrix(X1.c, nrow = ngrid)
Xcov1 <- tcrossprod(Xmat1) / nregion


# 2-d smoothing of the raw covariance
diag(Xcov1) <- NA
cov.mat1 <- gam(as.vector(Xcov1) ~ te(x0, x1, k = 10, bs = "ps"), method = "REML")
grids2d <- data.frame(x0 = x0, x1 = x1)
covmat.c1 <- matrix(predict(cov.mat1, newdata = grids2d), nrow = ngrid) 
cov1 <- (covmat.c1 + t(covmat.c1)) / 2 #  Symmetrize covariance function

# Estimate measurement errors
Xcov1 <- tcrossprod(Xmat1) / nregion 
diag <- diag(Xcov1)
loess_diag <- suppressWarnings(loess.as(gridPoints, diag, degree = 1,
                                        criterion = "gcv", user.span = NULL, plot = F))  
sig1Est <- mean(loess_diag$fitted - diag(cov1))

# Eigen decomposition
eigen_temp11 <- eigen(cov1, symmetric = TRUE)
eigen_temp11$values <- eigen_temp11$values[which(eigen_temp11$values > 0)]  # Obtain positive eigenvalues
eigen_temp11$vectors <- eigen_temp11$vectors[, 1:length(eigen_temp11$values)] 

# Normalizing
for(e in 1:length(eigen_temp11$values)){ 
  #e=1
  normal.factor <- trapz(gridPoints, eigen_temp11$vectors[, e]^2)
  eigen_temp11$vectors[, e] <- eigen_temp11$vectors[, e] / sqrt(2*normal.factor)
  eigen_temp11$values[e] <- eigen_temp11$values[e] * 2*normal.factor
}

# Number of eigencomponents in the 1st dimension
M1 <- length(which(cumsum(eigen_temp11$values) / sum(eigen_temp11$values) < FVEvalue)) + 1 

# Store the estimated univariate eigenfunctions in the first dimension
efun11 <- eigen_temp11$vectors[,1:M1]

# Estimate univariate PC scores in the first dimension
Yc11  <- t(Xmat1) # Centered data
e.scores11 <- matrix(0, nrow = nregion, ncol = M1)
for(i in 1:nregion){
  e.scores11[i,] <- apply(matrix(rep(Yc11[i,], M1), nrow = ngrid) * efun11[,1:M1],2, mean) 
}

## FPCA in the 2nd dimension 
dat.s2 <- df[,c("rid", "t", "dt", "y2")]
colnames(dat.s2) <- c("id", "gridPoints", "dt", "Y2")

# Calculate the raw covariances
Xmat2 <- matrix(X2.c, nrow = ngrid)
Xcov2 <- tcrossprod(Xmat2) / nregion


# 2-d smoothing of the raw covariance
diag(Xcov2) <- NA
cov.mat2 <- gam(as.vector(Xcov2) ~ te(x0, x1, k = 10, bs = "ps"), method = "REML")
grids2d <- data.frame(x0 = x0, x1 = x1)
covmat.c2 <- matrix(predict(cov.mat2, newdata = grids2d), nrow = ngrid) 
cov2 <- (covmat.c2 + t(covmat.c2)) / 2 #  Symmetrize covariance function

# Estimate measurement errors
Xcov2 <- tcrossprod(Xmat2) / nregion  
diag2 <- diag(Xcov2)
loess_diag <- suppressWarnings(loess.as(gridPoints, diag2, degree = 1,
                                        criterion = "gcv", user.span = NULL, plot = F))  
sig2Est <- mean(loess_diag$fitted - diag(cov2))

# Eigen decomposition
eigen_temp22 <- eigen(cov2, symmetric = TRUE)
eigen_temp22$values <- eigen_temp22$values[which(eigen_temp22$values > 0)]  # Obtain positive eigenvalues
eigen_temp22$vectors <- eigen_temp22$vectors[, 1:length(eigen_temp22$values)] 

# Normalizing
for(e in 1:length(eigen_temp22$values)){ 
  normal.factor <- trapz(gridPoints, eigen_temp22$vectors[, e]^2)
  eigen_temp22$vectors[, e] <- eigen_temp22$vectors[, e] / sqrt(2*normal.factor)
  eigen_temp22$values[e] <- eigen_temp22$values[e] * 2*normal.factor
}

# Number of eigencomponents in the 1st dimension
M2 <- length(which(cumsum(eigen_temp22$values) / sum(eigen_temp22$values) < FVEvalue)) + 1 

# Store the estimated univariate eigen functions in the 2nd dimension
efun22 <- eigen_temp22$vectors[,1:M2]

# Estimate univariate PC scores in the 2nd dimension
Yc22       <- t(Xmat2) # Centered data
e.scores22 <- matrix(0, nrow = nregion, ncol = M2)
for(i in 1:nregion){
  e.scores22[i,] <- apply(matrix(rep(Yc22[i,], M2), nrow = ngrid) * efun22[,1:M2],2, mean) 
}


## Construct output
out <- list(dat.s1 = dat.s1, dat.s2 = dat.s2, 
            sig1Est = sig1Est, sig2Est = sig2Est,
            efun11 = efun11, efun22 = efun22, 
            e.scores11 = e.scores11, e.scores22 = e.scores22,
            M1 = M1, M2 = M2)
return(out)
}





  


