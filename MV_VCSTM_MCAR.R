MV_VCSTM_MCAR <- function(FPCAout, # Output from function MV_VCSTM_univariate_decomposition, including the follows:
                           # dat.s1: data.frame in first dimension in long format with four labeled columns (described below)
                            # DATA.FRAME COLUMNS: 
                             # id: region IDs (vector of length ngrid*nregion) 
                             # gridPoints: follow-up time points (vector of length ngrid*nregion) 
                             # Y1: observed hospitalization rate data (vector of length ngrid*nregion) 
                             # dt: order of the follow-up time points (vector of length ngrid*nregion) 
                           # dat.s2: data.frame in second dimension in long format with four labeled columns (as dat.s1)
                           # sig1Est: estimated measurement error variance in the first dimension (scalar)
                           # sig2Est: estimated measurement error variance in the second dimension (scalar)
                           # efun11: estimated univariate eigenfunctions in the first dimension (matrix of dimension ngrid*M1) 
                           # efun22: estimated univaraite eigenfunctions in the second dimension (matrix of dimension ngrid*M2) 
                           # e.scores11: estimated univariate PC eigenscores in the first dimension (matrix of dimension nregion*M1)
                           # e.scores22: estimated univariate PC eigenscores in the second dimension (matrix of dimension nregion*M2)
                           # M1: estimated number of eigencomponents chosen by FVE in the first diemension (scalar)
                           # M2: estimated number of eigencomponents chosen by FVE in the second diemension (scalar)
                          Adj.Mat # Adjacency matrix from the map (0-1 matrix of dimension nregion*nregion)
){
#############################################################################
## Description: Function for the first MCMC (Estimation Algorithm steps 3-4 in Section 2.2) 
##              including estimation the between eigencomponent matrix, as well as mutivariate eigenfucntions
## Definition:  nregion: number of regions, 
##              ngrid: number of follow-up time points,
##              M1: estimated number of eigen components in the first dimension,
##              M2: estimated number of eigen components in the second dimension.
## Args:        see above
## Returns:     list()
##              psi11Est: estimated multivariate eigenfunctions in the first dimension (matrix of dimension ngrid*(M1+M2))
##              psi22Est: estimated multivariate eigenfunctions in the second dimension (matrix of dimension ngrid*(M1+M2))
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
  
  
## Define time points, number of regions, and time points
nregion <- data.T$nregion  
ngrid <- data.T$ngrid
gridPoints <- data.T$gridPoints
  

# Eigencomponents chosen in each dimension
M1 <- FPCAout$M1
M2 <- FPCAout$M2

## WinBUGS model for the 1st MCMC (MCAR part)
MCARmodel = function(){
  # Likelihood
  for (i in 1:nregion){
    for (j in 1:noutcome){
      mat.Xi[i,j] ~ dnorm(Xi[i,j], invTao2)
    }
  }
  
  # Definition of the matrix Xi
  for (i in 1:nregion){
    for (j in 1:noutcome){
      Xi[i,j] <- inprod(tPsi[1:j,i], roDis[j, 1:j])
    }
  }
  
  # Psi matrix(spatially dependent random effects)
  # Impose each row of tPsi matrix a CAR prior
  for (j in 1:noutcome){
    tPsi[j,1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], 1, nu)
  }
  
  # Impose Uniform prior for nu
  for(i in 1:nregion){ceros[i] <- 0}
  nu ~ dunif(nu.inf, nu.sup)
  nu.inf <- min.bound(Cw[], adj[], num[], Mv[])
  nu.sup <- max.bound(Cw[], adj[], num[], Mv[])
  
  
  # Prior distributions of correlation between diseases
  #(cell not used the Cholesky decompostion are fixed to 0)
  for (i in 1:(noutcome-1)){
    for (j in (i+1):noutcome){
      roDis[i,j] <- 0
    }
  }
  for (i in 1:noutcome){
      roDis[i,i] ~ dlnorm(-3, 1)
  }
  for (i in 2:noutcome){
    for (j in 1:(i-1)){
      roDis[i,j] ~ dnorm(0, 5)
    }
  }
  
  # Other priors
  invTao2 ~ dgamma(2, sigma.est)
  tao2 <- 1/invTao2
} 


# Collect all the data needed for WinBUGS
# Define the score matrix
mat.Xi <- cbind(FPCAout$e.scores11, FPCAout$e.scores22)

# Number of total eigen components from two dimensions
noutcome <- M1 + M2 

# Adjacency matrix related elements needed for the built-in `proper.car` function in WinBUGS
W <- Adj.Mat
adj.weights.num <- as.carAdjacency(W)
CM <- as.carCM(adj.weights.num$adj, adj.weights.num$weights, adj.weights.num$num)
Cw <- CM$C
Mv <- CM$M


# Construct D matrix from number of neighbors for each region
D <- rowSums(W) 
D <- diag(D)

# Eestimated error variance
sigma.est <- mean(c(FPCAout$sig1Est, FPCAout$sig2Est))

# Data list for MCMC
MCARdata <- list(mat.Xi = mat.Xi, adj = adj.weights.num$adj, 
                 num = adj.weights.num$num, Cw = Cw, Mv = Mv, 
                 nregion = nregion, noutcome = noutcome, sigma.est = sigma.est)

# Define initial values 
## Define a function to generate strings
string.gen <- function(noutcome){
  string <- c()
  for (k in 1:noutcome){
    string <- c(string, rnorm(k-1), runif(1,0,2), rep(NA, (noutcome-k)))
  }
  return(string)
}

## Define a function to generate initial values
MCARinits = function(){
  list(tPsi = matrix(rnorm(nregion * noutcome), nrow = noutcome), 
       nu = runif(1, 0, 1), 
       roDis = matrix(string.gen(noutcome), ncol = noutcome, byrow = T), 
       invTao2 = rnorm(1,1/sigma.est, 0.0001))
}


# Paramaters to save
MCARparam = c("tao2","nu", "roDis")

# Model fit paralelly using WinBUGS
MCARmcmc = pbugs(data = MCARdata, inits = MCARinits, par = MCARparam, model = MCARmodel, 
              n.iter = 15000, n.burnin = 5000, 
              DIC = F, bugs.seed = 1, debug = FALSE, n.thin = 2, 
              bugs.directory="~/winbugs14_full_patched/WinBUGS14/") # Need to be filled with the directory where WinBUGS is installed!

# Collect results from MCMC
sum <- MCARmcmc$summary[,c('mean')]
roDis <- matrix(rep(NA, noutcome*noutcome), ncol = noutcome)
# roDis
for (i in 1:noutcome){
  for (j in 1:noutcome){
    if (i ==j){
      roDis[i,i] <- sum[paste('roDis[',i,',',i,']',sep = '')]
    } else{
      roDis[i,j] <- sum[paste('roDis[',i,',',j,']', sep = '')]
      roDis[j,i] <- 0
    }
  }
}

# Estimate between eigencomponents matrix Sigma_b
sigb <- roDis %*% t(roDis)

# Estimate spatial smoothing parameter nu
MCAR.nuEst <- sum['nu']

# Estimate the error variance
tao2Est <- sum['tao2']

## Using MFPCA algorithm to estimate multivariate eigenfunctions
# Eigen analysis of Sigma_b matrix
z.eg = svd(sigb)

# Estimate multivariate eigenfunctions from their univariate counterparts
psi11 <- matrix(0, nrow = ngrid, ncol = (M1+M2))
psi22 <- matrix(0, nrow = ngrid, ncol = (M1+M2))
for(j in 1:(M1+M2)){
  # j = 1
  psi11[,j] = FPCAout$efun11 %*% z.eg$u[1:M1,j] 
  psi22[,j] = FPCAout$efun22 %*% z.eg$u[(M1+1):(M1+M2),j]
}

# Normalizing
for(e in 1:(M1+M2)){ 
  normal.factor <- trapz(gridPoints, psi11[, e]^2)
  psi11[, e] <- psi11[, e] / sqrt(2*normal.factor)
}

for(e in 1:(M1+M2)){ 
  normal.factor <- trapz(gridPoints, psi22[, e]^2)
  psi22[, e] <- psi22[, e] / sqrt(2*normal.factor)
}

psi11Est <- psi11
psi22Est <- psi22


# Construct output
out <- list(psi11Est = psi11Est, psi22Est = psi22Est)

return(out)
}


