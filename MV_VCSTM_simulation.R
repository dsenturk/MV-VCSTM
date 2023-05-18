MV_VCSTM_simulation <- function(numRegion=2, # Number of regions (scalar, input 1 if you want 367 regions, input 2 if you want 49 regions)
                                sigma=0.02 # Measurement error variance (scalar)
){
#############################################################################
## Description: Function for simulating one data set under the simulation design described
##              in Section 4.
## Args: see above
## Definition: nregion: number of regions; 
##             ngrid: number of follow-up time points; 
##             M: number of eigen components in each dimension.
## Returns: list()
#           data, data.frame with columns c("rid", "y1", "y2", "t", "dt", "x1", "x2"), for further analysis using MV-VCSTM
#           DATA.FRAME COLUMNS (data is stored in long format): 
#             rid: region IDs (vector of length ngrid*nregion) 
#             y1: observed outcome is the 1st dimension (vector of length ngrid*nregion)
#             y2: observed outcome is the 2nd dimension (vector of length ngrid*nregion)
#             t: follow-up time points (vector of length ngrid*nregion) 
#             dt: order of follow-up time points in each region (vector of length ngrid*nregion)   
#             x1: 1st covariate (vector of length ngrid*nregion)
#             x2: 2nd covariate (vector of length ngrid*nregion)
#  
#           Adj.Mat: adjacent matrix (matrix of dimension nregion*nregion)
#  
#           dat.True: list of all true values, which contains:
#             nregion: number of regions (scalar)
#             ngrid: number of follow-up time points (scalar)
#             gridPoints: follow up time points (vector of length ngrid)
#             M: number of PC components in each dimension (scalar)
#             psi1.True: multivariate eigenfunctions in the first dimension (matrix of dimension ngrid*M)
#             psi2.True: multivariate eigenfunctions in the second dimension (matrix of dimension ngrid*M)
#             nu.True: spatial smoothing parameter (scalar)
#             alpha1.True, alpha2.True, alpha3.True: spatial variance components (scalar)
#             sigma.True: measurement error variance (scalar)
#             region.PC: region-specific PC scores (matrix of dimension nregion*M)
#             Xmat: covariate matrix containing two covariates as two conlumns for each region (matrix of dimension nregion*2)
#             df.True, data.frame with columns c("rid", "y1.True", "y2.True", "t", "r.eff1.1","r.eff1.2","f.eff1.3", "r.eff2.1","r.eff2.2","f.eff2.3"),
#             for storing the true trajectories that will be used for evaluation of multivariate trajectory predictions and visualization
#             DATA.FRAME COLUMNS: 
#              rid: region IDs (vector of length ngrid*nregion)
#              y1.True: true hospitalization rate data (vector of length ngrid*nregion)
#              y2.True: true mortality rate data (vector of length ngrid*nregion)
#              t: follow-up time points (vector of length ngrid*nregion) 
#              r.eff1.1: first component of region-specific deviation in the first dimension U^(1)(t) (vector of length ngrid*nregion)
#              r.eff1.2: second component of region-specific deviation in the first dimension U^(1)(t) (vector of length ngrid*nregion)
#              r.eff1.3: third component of region-specific deviation in the first dimension U^(1)(t) (vector of length ngrid*nregion)
#              r.eff2.1: first component of region-specific deviation in the second dimension U^(2)(t) (vector of length ngrid*nregion)
#              r.eff2.2: second component of region-specific deviation in the second dimension U^(2)(t) (vector of length ngrid*nregion)
#              r.eff2.3: third component of region-specific deviation in the second dimension U^(2)(t) (vector of length ngrid*nregion)

  
################################################################################## 

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
  
# Define the grid points used for eigenfunctions
ngrid <- 24 # 2 year follow up
M <- 3 # number of eigen components in each dimension

# Generate multivariate eigenfunctions in each dimension
fbs <- create.fourier.basis(rangeval = c(0, 2), nbasis = 10)
gridPoints <- seq(0, 1, length.out = ngrid)

# Multivariate eigenfunction in the 1st dimension 
psi1.True <- eval.basis(gridPoints, fbs)[,c(4, 5, 8)]

# Multivariate eigenfunction in the 2nd dimension 
gridPoints2 <- seq(1, 2, length.out = ngrid)
psi2.True <- psi22.True <- eval.basis(gridPoints2, fbs)[,c(4, 9, 8)]
for(e in 1:M){ 
  set.seed(123+e)
  psi2.True[,e] <- (-1)^sample(1:2, 1, replace = T) * psi22.True[,e]
}


# Load adjacency matrix from maps
if(numRegion==1){
  # Merged HSA map
  load("AdjMat_367.RData")
  # Adjacency matrix from merged HSA map
  W <- AdjMat_367
} else {
  # US states map 
  load("49StateAdjMat.RData")
}

# Number of total regions
nregion <- nrow(W)

# Number of neighbors for each region
D <- rowSums(W) 


# Spatial correlaion parameter
nu.True <- 0.97

# Covariance matrix of region-specific multivariate PC scores
covmat <- solve(diag(D) - nu.True * W)

# Spatial variance component \alpha_l
alpha1.True <- 2.4  
alpha2.True <- 1.2
alpha3.True <- 0.6
alpha.True <- c(alpha1.True, alpha2.True, alpha3.True)

# Measurement error
sigma.True <- sigma

# Varying coefficient functions
beta0.1 <- function(x){
  return(0.4*exp(2*x-1))
}

beta1.1 <- function(x){
  return(x*(1-x))
}

beta2.1 <- function(x){
  return(0.25*(1-x)^2)
}

# 2nd dimension
beta0.2 <- function(x){
  return(0.5*sqrt(2*x)+0.2)
}

beta1.2 <- function(x){
  return((x-0.5)^2)
}

beta2.2 <- function(x){
  return(0.25*x^2)
}


#######################################################
# Construct Data set
#######################################################
# Region-specific deviation 
# Generate region-specific PC scores from multivariate normal distribution
set.seed(12)
region.PC1 <- mvrnorm(1, rep(0, nregion), covmat * alpha1.True)  
region.PC2 <- mvrnorm(1, rep(0, nregion), covmat * alpha2.True)  
region.PC3 <- mvrnorm(1, rep(0, nregion), covmat * alpha3.True)
# Store the generated PC scores into `region.PC` matrix
region.PC <- array(rep(0, nregion * M), c(nregion, M))
region.PC[,1] <- region.PC1
region.PC[,2] <- region.PC2
region.PC[,3] <- region.PC3

# Create data.frame to store dataset

df <- data.frame(matrix(ncol = 15, nrow = nregion * ngrid))
colnames(df) <- c("rid","y1", "y1.True", "y2", "y2.True", "t","dt","r.eff1.1","r.eff1.2","r.eff1.3",
                  "r.eff2.1","r.eff2.2","r.eff2.3", "x1", "x2")
df$rid <- rep(1:nregion, each = ngrid)
df$t <- rep(gridPoints,nregion)
df$dt <- rep(1:ngrid,nregion)

# Region-specific deviation 
# 1st dimension
df$r.eff1.1 <- region.PC[,1][df$rid] * psi1.True[,1][df$dt]
df$r.eff1.2 <- region.PC[,2][df$rid] * psi1.True[,2][df$dt]
df$r.eff1.3 <- region.PC[,3][df$rid] * psi1.True[,3][df$dt]
# 2nd dimension
df$r.eff2.1 <- region.PC[,1][df$rid] * psi2.True[,1][df$dt]
df$r.eff2.2 <- region.PC[,2][df$rid] * psi2.True[,2][df$dt]
df$r.eff2.3 <- region.PC[,3][df$rid] * psi2.True[,3][df$dt]


# Generate region-level covariates Xi
Cov <- diag((1-2^(-0.5)), 2) + 2^(-0.5)*c(1, 1) %*% t(c(1,1))
# Generate region-specific covariates from multivariate normal distribution
Xmat <- matrix(rep(NA, nregion*2), nrow = nregion)
set.seed(12)
for (i in 1:nregion){
  Xmat[i,] <- mvrnorm(1, rep(0, 2), Cov) 
}

df$x1 <- Xmat[,1][df$rid]
df$x2 <- Xmat[,2][df$rid]


# Generate outcome
# 1st dimension
set.seed(12)
measure.err1 <- rnorm(dim(df)[1], 0, sqrt(sigma.True))
df$y1 <- beta0.1(df$t) + df$x1 * beta1.1(df$t) + df$x2 * beta2.1(df$t) + 
  df$r.eff1.1 + df$r.eff1.2 + df$r.eff1.3 + measure.err1
df$y1.True <- beta0.1(df$t) + df$x1 * beta1.1(df$t) + df$x2 * beta2.1(df$t) + 
  df$r.eff1.1 + df$r.eff1.2 + df$r.eff1.3 


# 2nd dimension
set.seed(13)
measure.err2 <- rnorm(dim(df)[1], 0, sqrt(sigma.True))
df$y2 <- beta0.2(df$t) + df$x1 * beta1.2(df$t) + df$x2 * beta2.2(df$t) + 
  df$r.eff2.1 + df$r.eff2.2 + df$r.eff2.3 + measure.err2
df$y2.True <- beta0.2(df$t) + df$x1 * beta1.2(df$t) + df$x2 * beta2.2(df$t) + 
  df$r.eff2.1 + df$r.eff2.2 + df$r.eff2.3


# Store all the related true values into `data.True` list
data.True <- list(nregion = nregion, ngrid = ngrid, gridPoints = gridPoints, 
                  M = M, psi1.True = psi1.True, psi2.True = psi2.True, 
                  nu.True = nu.True, alpha1.True = alpha1.True, 
                  alpha2.True = alpha2.True, alpha3.True = alpha3.True, 
                  sigma.True = sigma.True, region.PC = region.PC, Xmat = Xmat,
                  df.True = df[,c(1, 3, 5:6, 8:13)])

# Generate output
out <- list(data = df[,c(1:2,4,6:7, 14:15)], Adj.Mat = W, data.True = data.True)

return(out)

}












