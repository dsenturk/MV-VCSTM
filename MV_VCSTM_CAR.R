MV_VCSTM_CAR <- function(FPCAout, # Output from function MST_FPCA_univariate_decomposition, including the follows:
                           # dat.s1: data.frame in first dimension in long format with four labeled columns (described below)
                            # DATA.FRAME COLUMNS: 
                            # id: region IDs (vector of length ngrid*nregion) 
                            # gridPoints: follow-up time points (vector of length ngrid*nregion) 
                            # Y1: observed hospitalization rate data (vector of length ngrid*nregion) 
                            # dt: follow-up time points in order (vector of length ngrid*nregion) 
                           # dat.s2: data.frame in second dimension in long format with four labeled columns (as dat.s1)
                           # sig1Est: estimated error variance in the first dimension (scalar)
                           # sig2Est: estimated error variance in the second dimension (scalar)
                           # efun11: estimated univariate eigenfunctions in the first dimension (matrix of dimension ngrid*M1) 
                           # efun22: estimated univaraite eigenfunctions in the second dimension (matrix of dimension ngrid*M2) 
                           # e.scores11: estimated univariate PC eigenscores in the first dimension (matrix of dimension nregion*M1)
                           # e.scores22: estimated univariate PC eigenscores in the second dimension (matrix of dimension nregion*M2)
                           # M1: estimated number of eigencomponents chosen by FVE in the first diemension (scalar)
                           # M2: estimated number of eigencomponents chosen by FVE in the second diemension (scalar)
                         MCARout, # Output from function MV_VCSTM_MCAR, including the following:
                           # psi11Est: estimated multivariate eigenfunctions in the first dimension (matrix of dimension ngrid*(M1+M2))
                           # psi22Est: estimated multivariate eigenfunctions in the second dimension (matrix of dimension ngrid*(M1+M2))
                         L, # An integer between 1 and (M1+M2) defining the number of eigencomponents included in the final CAR fitting.
                            # The estimated spatial variance parameters `alphaEst` are used as a guidance in choosing L in applications.
                            # In our implementation, we choose L = M1 + M2.
                         Adj.Mat # Adjacency matrix from the map (0-1 matrix of dimension nregion*nregion)
){
  
#############################################################################
## Description: Function for the second MCMC (Estimation Algorithm steps 4 in Section 2.2) 
##              including estimation of multivariate PC scores.
## Definition:  nregion: number of regions, 
##              ngrid: number of follow-up time points,
##              M: number of eigen components in each dimension.
## Args:        see above
## Returns:     list()
##              xiEst: estimates of region-specific PC score matrix (matrix of dimension nregion*(M1+M2))
##              alphaEst: estimates of spatial variance parameters (vector of length (M1+M2))
##              beta1Est: estimate of VCFs in the 1st dimension (matrix of dimension ngrid * 3)
##              beta2Est: estimate of VCFs in the 2nd dimension (matrix of dimension ngrid * 3)
##              beta1.mc: all posterior samples of VCFs in the 1st dimension (array of dimension (1000*ngrid*3))
##              beta2.mc: all posterior samples of VCFs in the 2nd dimension (array of dimension (1000*ngrid*3))
##              sigma1Est: estimate of measurement error variance in the 1st dimension (scalar)
##              sigma2Est: estimate of measurement error variance in the 2nd dimension (scalar)
##              nuEst: estimate of spatial smoothing parameter (scalar)
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

# Input the generated two dimensional data
# 1st dimension
Y1 <- matrix(data$y1, nrow = nregion, byrow = TRUE)
# 2nd dimension
Y2 <- matrix(data$y2, nrow = nregion, byrow = TRUE)

# Covariates
x0 <- rep(1, nregion)
x1 <- data.T$Xmat[,1]
x2 <- data.T$Xmat[,2]

# Eigencomponents chosen in each dimension
M1 <- FPCAout$M1
M2 <- FPCAout$M2

###########################################################
## Using MCMC (2nd MCMC) to estimate varying coefficient functions, 
## as well as multivariate PC scores by MCMC in WinBUGS
#########################################################

# WinBUGS model for the 2nd MCMC (CAR part)

Spline.CARmodel <- function()
{                           
  #Model structure
  for (i in 1:nregion){
    for (j in 1:nmonth){
      Y[i, j] ~ dnorm(Mean1[i,j], tauepsilon1)
      Mean1[i,j] <- x0[i]*(inprod2(a01[1:2], V1[j,1:2]) + inprod2(b01[1:nknots], Z[j,1:nknots])) + 
        x1[i]*(inprod2(a11[1:2], V1[j,1:2]) + inprod2(b11[1:nknots], Z[j,1:nknots])) +
        x2[i]*(inprod2(a21[1:2], V1[j,1:2]) + inprod2(b21[1:nknots], Z[j,1:nknots])) +  
        inprod2(xi[1:noutcome,i], psi1.mat[j,1:noutcome])
      
      
      Y[i+nregion, j] ~ dnorm(Mean2[i,j], tauepsilon2)
      Mean2[i,j] <- x0[i]*(inprod2(a02[1:2], V1[j,1:2]) + inprod2(b02[1:nknots], Z[j,1:nknots])) + 
        x1[i]*(inprod2(a12[1:2], V1[j,1:2]) + inprod2(b12[1:nknots], Z[j,1:nknots])) +
        x2[i]*(inprod2(a22[1:2], V1[j,1:2]) + inprod2(b22[1:nknots], Z[j,1:nknots])) + 
        inprod2(xi[1:noutcome,i], psi2.mat[j,1:noutcome])
      
    }
  }
  
  
  # Impose CAR prior on mutivariate PC scores xi's
  for (k in 1:noutcome){
    xi[k,1:nregion] ~ car.proper(ceros[], Cw[], adj[], num[], Mv[], invAlpha[k], nu)
    invAlpha[k] ~ dgamma(1, 1) 
    alpha[k] <- 1/invAlpha[k]
  }
  
  
  # Impose Uniform prior on the spatial smoothing parameter
  nu ~ dunif(nu.inf, nu.sup)
  
  # Define constants
  for(i in 1:nregion){ceros[i] <- 0}
  nu.inf <- min.bound(Cw[], adj[], num[], Mv[])
  nu.sup <- max.bound(Cw[], adj[], num[], Mv[])  
  
  
  #Distribution of the truncated polynomial coefficients
  for (i in 1:nknots){
    b01[i]~dnorm(0,taub)
    b11[i]~dnorm(0,taub)
    b21[i]~dnorm(0,taub)
    b02[i]~dnorm(0,taub)
    b12[i]~dnorm(0,taub)
    b22[i]~dnorm(0,taub)
  }
  
  #Prior distributions for monomial coefficients
  for (l in 1:2){
    a01[l]~dnorm(0,1.0E-6)
    a11[l]~dnorm(0,1.0E-6)
    a21[l]~dnorm(0,1.0E-6)
    a02[l]~dnorm(0,1.0E-6)
    a12[l]~dnorm(0,1.0E-6)
    a22[l]~dnorm(0,1.0E-6)
  }
  
  
  #Prior distributions for precision parameters
  taub ~ dgamma(1.0E-6,1.0E-6)
  tauepsilon1 ~ dgamma(1.0E-6,1.0E-6)
  tauepsilon2 ~ dgamma(1.0E-6,1.0E-6)
  
  
  
  #Deterministic transformations
  sigmaepsilon1 <- 1/sqrt(tauepsilon1)
  sigmaepsilon2 <- 1/sqrt(tauepsilon2)
  sigmab <- 1/sqrt(taub)
  
  for (j in 1:nmonth){
    beta01[j] <- a01[1]*V1[j,1] + a01[2]*V1[j,2] + b01[1]*Z[j,1]+b01[2]*Z[j,2]+b01[3]*Z[j,3]
    beta11[j] <- a11[1]*V1[j,1] + a11[2]*V1[j,2] + b11[1]*Z[j,1]+b11[2]*Z[j,2]+b11[3]*Z[j,3]
    beta21[j] <- a21[1]*V1[j,1] + a21[2]*V1[j,2] + b21[1]*Z[j,1]+b21[2]*Z[j,2]+b21[3]*Z[j,3]
    
    beta02[j] <- a02[1]*V1[j,1] + a02[2]*V1[j,2] + b02[1]*Z[j,1]+b02[2]*Z[j,2]+b02[3]*Z[j,3]
    beta12[j] <- a12[1]*V1[j,1] + a12[2]*V1[j,2] + b12[1]*Z[j,1]+b12[2]*Z[j,2]+b12[3]*Z[j,3]
    beta22[j] <- a22[1]*V1[j,1] + a22[2]*V1[j,2] + b22[1]*Z[j,1]+b22[2]*Z[j,2]+b22[3]*Z[j,3]
    
  } 
}

# Collect data needed for WinBUGS
# Number of total follow-up months
nmonth <- ngrid

# Total number of eigencomponents from two dimensions
noutcome <- L

# Adjacency matrix related elements needed for the built-in `proper.car` function in WinBUGS
W <- Adj.Mat
adj.weights.num <- as.carAdjacency(W)
CM <- as.carCM(adj.weights.num$adj, adj.weights.num$weights, adj.weights.num$num)
Cw <- CM$C
Mv <- CM$M

# Construct D matrix from number of neighbors for each region
D <- rowSums(W) 
D <- diag(D)

# Constructed outcome matrix from elements in two dimensions
Y <- rbind(Y1, Y2)


# Estimated multivariate eigenfunctions
# 1st dimension
psi1.mat <- MCARout$psi11Est
# 2nd dimension
psi2.mat <- MCARout$psi22Est


# Model VCF via thin-plate splines
# Define the matrix of fixed effects
V1 <- cbind(rep(1,ngrid),gridPoints)

# Define knots used for the p-spline at equal spaced quantiles of the covariate
nknots <- 3 # Input the number of knots
knots <- quantile(unique(gridPoints), seq(0,1,length=(nknots+2))[-c(1,(nknots+2))])

# Obtain the design matrix for random coefficients using a radial basis
Z_K <- (abs(outer(gridPoints,knots,"-")))^3
OMEGA_all <- (abs(outer(knots,knots,"-")))^3
svd.OMEGA_all <- svd(OMEGA_all)
sqrt.OMEGA_all <- t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
Z <- t(solve(sqrt.OMEGA_all,t(Z_K)))


# WinBUGS data list
bug_data <- list(nmonth = nmonth, nregion = nregion, noutcome = noutcome,
                 V1 = V1, Z = Z, nknots = nknots, Y = Y, 
                 x0=x0, x1=x1, x2=x2, 
                 adj = adj.weights.num$adj, num = adj.weights.num$num, 
                 Cw = Cw, Mv = Mv, 
                 psi1.mat = psi1.mat, psi2.mat = psi2.mat)

# Initital values
bug_init <- function(){
  list(a01=c(0,0),a11=c(0,0),a21=c(0,0),
       a02=c(0,0),a12=c(0,0),a22=c(0,0),
       b01=rep(0,nknots),b11=rep(0,nknots),b21=rep(0,nknots),
       b02=rep(0,nknots),b12=rep(0,nknots),b22=rep(0,nknots),
       invAlpha = rep(1, noutcome), nu = rnorm(1, .9, .0001), 
       tauepsilon1=50, tauepsilon2=50,
       xi = matrix(rnorm(nregion * noutcome, 0, 0.0001), nrow = noutcome),
       taub = 0.01)
}



# Paramaters to save
bug_param <- c("beta01", "beta11", "beta21", "beta02", "beta12", "beta22",
               "alpha", "sigmaepsilon1","sigmaepsilon2", "nu", "xi")


# Model fit paralelly using WinBUGS
n.iter = 15000; n.burnin = 5000; n.thin = 10
Carfit = pbugs(data = bug_data, inits = bug_init, par = bug_param, 
            model = Spline.CARmodel, n.iter = n.iter, n.burnin = n.burnin, DIC = F, 
            bugs.seed = 1, debug = FALSE, n.thin = n.thin, n.chains = 1,
            bugs.directory="~/winbugs14_full_patched/WinBUGS14/") # Need to be filled with the directory where WinBUGS is installed!


# Obtain MCMC sample means
outsum <- Carfit$summary[,c('mean')]

# Estimate VCFs
beta1Est <- beta2Est <- matrix(rep(NA, ngrid*3), nrow = ngrid)
for (i in 1:ngrid){
  beta1Est[i,1] <- outsum[paste('beta01[',i,']', sep = '')]
  beta1Est[i,2] <- outsum[paste('beta11[',i,']', sep = '')]
  beta1Est[i,3] <- outsum[paste('beta21[',i,']', sep = '')]
  beta2Est[i,1] <- outsum[paste('beta02[',i,']', sep = '')]
  beta2Est[i,2] <- outsum[paste('beta12[',i,']', sep = '')]
  beta2Est[i,3] <- outsum[paste('beta22[',i,']', sep = '')]
}


# Estimate region-specific multivariate PC scores
xiEst <- matrix(rep(NA, nregion*noutcome), nrow = noutcome)
for (i in 1:noutcome){
  for (j in 1:nregion){
    xiEst[i,j] <- outsum[paste('xi[',i,',',j,']', sep = '')]
  }
}
xiEst <- t(xiEst)

# Estimate spatial variance papameter \alpha_\ell's
alphaEst <- rep(0, L)
for (i in 1:L){
  alphaEst[i] <- outsum[paste('alpha[',i,']', sep = '')]
}

# Estimate measurement error variance 
sigma1Est <- outsum["sigmaepsilon1"] # 1st dimension
sigma2Est <- outsum["sigmaepsilon2"] # 2nd dimension

# Estimat spatial smoothing parameter \nu
nuEst <- outsum["nu"]

# Collect estimated VCFs from all MCMC samples (used for inference)
beta1.mc <- array(rep(NA, ((n.iter - n.burnin)/n.thin)*ngrid*3), 
                       c((n.iter - n.burnin)/n.thin, ngrid, 3))
beta2.mc <- array(rep(NA, ((n.iter - n.burnin)/n.thin)*ngrid*3), 
                       c((n.iter - n.burnin)/n.thin, ngrid, 3))

beta1.mc[,,1] <- Carfit$sims.list$beta01
beta1.mc[,,2] <- Carfit$sims.list$beta11
beta1.mc[,,3] <- Carfit$sims.list$beta21
beta2.mc[,,1] <- Carfit$sims.list$beta02
beta2.mc[,,2] <- Carfit$sims.list$beta12
beta2.mc[,,3] <- Carfit$sims.list$beta22



## Construct output
out <- list(xiEst = xiEst, alphaEst = alphaEst, 
            beta1Est = beta1Est, beta2Est = beta2Est,
            beta1.mc = beta1.mc, beta2.mc = beta2.mc,
            sigma1Est = sigma1Est, sigma2Est = sigma2Est, 
            nuEst = nuEst)
return(out)
}



