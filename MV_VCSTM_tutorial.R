## MV_VCSTM_tutorial.R
#############################################################################
## Description: A step-by-step implementation of MV-VCSTM and the associated  
## procedures described in "Multivariate Varying Coefficient Spatiotemporal 
## Model".  
#############################################################################
## Functions implemented: 
## MV_VCSTM_simulation.R (Simulate a two-dimensional dataset), 
## MV_VCSTM_univariate_decomposition.R (Step 1-2 of the estimation algorithm),
## MV_VCSTM_MCAR.R (Step 3-4), 
## MV_VCSTM_CAR.R (Step 5), 
## MV_VCSTM_inference.R (Construct simultaneous credible intervals for varying coefficient functions (VCFs), as well as predictions).
#############################################################################
## Tutorial Outline:
## 1. Simulate two-dimensional outcome data (MV_VCSTM_simulation.R)
## 2. Perform MV_VCSTM estimation (MV_VCSTM_univariate_decomposition.R, MV_VCSTM_MCAR.R, MV_VCSTM_CAR.R)
## 3. Provide inference on varying coefficient functions (VCFs) (MV_VCSTM_inference.R)
## 4. Visualization of MV_VCSTM results
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

#############################################################################
# 1. Simulate two-dimensional data
#############################################################################

# Simulate one dataset from the simulation design described 
data.G <- MV_VCSTM_simulation(numRegion = 2, sigma = .02)  # MST_FPCA_simulation.R (numRegion = 1 if you want 367 regions, 
                                                                                  # numRegion = 2 if you want 49 regions,
                                                                                  # sigma = 0.02 where 0.02 can be any scalar representing the setup error variance)

# Data frame used for estimation and inference
data <- data.G[[1]]

# Adjacency matrix from the map
Adj.Mat <- data.G[[2]]

# Save the underlying true values when generating the data
data.T <- data.G[[3]]

#############################################################################
# 2. Perform MV-VCSTM estimation
#############################################################################

# NOTE: Performing MV-VCSTM estimation steps 1-2 (univariate FPCA decomposition) with 49 or 367 regions will take about one to two minutes.

FPCAout <- MV_VCSTM_univariate_decomposition(data, FVEvalue = 0.99)  # MV_VCSTM_univariate_decompostion.R
# Note that FVEvalue should be a scalar between 0 to 1, which is to choose the number of eigencomponents in each dimension. 
# To retain enough information at the initial step of the algorithm, it is better to use a large FVEvalue (close to 1). 
# In our implementation, we use FVEvalue = 0.99.

# Number of eigencomponents chosen in each dimension
M1 <- FPCAout$M1  # 1st dimension
M2 <- FPCAout$M2  # 2nd dimension


# NOTE: Performing MV-VCSTM estimation steps 3-4 (1st MCMC and estimate multivariate eigenfunctions) with 49 regions will take about one minute (367 regions will take about 6 minutes).

MCARout <- MV_VCSTM_MCAR(FPCAout, Adj.Mat) # MV-VCSTM_MCAR.R  

# NOTE: Performing MV-VCSTM estimation steps 5 (2nd MCMC to estimate multivariate PC scores) with 49 regions will take about 12 minutes (367 regions will take about 90 minutes).

CARout <- MV_VCSTM_CAR(FPCAout, MCARout, L = M1 + M2, Adj.Mat) # MV_VCSTM_CAR.R

# Note that L needs to be an integer between 1 and (M1 + M2), which can be defined by the user. 
# L is defining the number of eigencomponents included in the final CAR fitting.
# The estimated spatial variance parameters alphaEst are used as a guidance in choosing L in applications as mentioned in Section 2 of the paper.
# In our implementation, we choose L = M1 + M2.


#############################################################################
# 3. Inferences on varying coefficient functions (VCFs) and make predictions
#############################################################################
# NOTE: Performing inferences on VCFs and making predictions with 49 or 367 regions will take less than one minute.

Inferenceout <- MV_VCSTM_inference(CARout) # MV_VCSTM_inference.R

#############################################################################
# 4. Visualization of MST-FPCA results
#############################################################################  
# Define the grid points used for the multivariate eigenfunctions
ngrid <- data.T$ngrid # 2 year follow up
gridPoints <- data.T$gridPoints

# Number of components chosen from the 1st and 2nd dimension
M1 <- FPCAout$M1
M2 <- FPCAout$M2

# True varying coefficient functions (VCFs)
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

# True multivariate eigenfunctions
psi1.True <- data.T$psi1.True # 1st dimension
psi2.True <- data.T$psi2.True # 2nd dimension

###########################################################################
# Plot estimates of varying coefficient functions 
# and 95% simultaneous confidence bands
###########################################################################

par(mfrow=c(2,3))

# Plot beta0.1(t)
beta01CB <- Inferenceout$beta1.CB[,c(1,2)]
plot(gridPoints, CARout$beta1Est[,1], 'l', col = "black", lwd = 2, 
     main = expression(paste("(a)"," ", widehat(beta)[0]^{(1)}, "(t)")),
     ylim = range(beta01CB), xlim = c(0,1), xaxs = "i", cex.main = 1.8, 
     xlab = "", ylab = "", cex.axis = 1.4) 
title(xlab = "t", ylab = "", 
      line = 2, cex.lab = 1.5)
lines(gridPoints, beta0.1(gridPoints), lwd = 2, lty = 1, col = "grey52")
lines(gridPoints, beta01CB[,1], lwd = 2, lty = 2)
lines(gridPoints, beta01CB[,2], lwd = 2, lty = 2)

# Plot beta1.1(t)
beta11CB <- Inferenceout$beta1.CB[,c(3,4)]
plot(gridPoints, CARout$beta1Est[,2], 'l', col = "black", lwd = 2, 
     main = expression(paste("(b)"," ", widehat(beta)[1]^{(1)}, "(t)")),
     ylim = range(beta11CB), xlim = c(0,1), xaxs = "i", cex.main = 1.8,
     xlab = "", ylab = "", cex.axis = 1.4) 
title(xlab = "t", ylab = "", 
      line = 2, cex.lab = 1.5)
lines(gridPoints,beta1.1(gridPoints), lwd = 2, lty = 1, col = "grey52")
lines(gridPoints, beta11CB[,1], lwd = 2, lty = 2)
lines(gridPoints, beta11CB[,2], lwd = 2, lty = 2)

# Plot beta2.1(t)
beta21CB <- Inferenceout$beta1.CB[,c(5,6)]
plot(gridPoints, CARout$beta1Est[,3], 'l', col = "black", lwd = 2 ,
     main = expression(paste("(c)"," ", widehat(beta)[2]^{(1)}, "(t)")),
     ylim = range(beta21CB), xlim = c(0,1), xaxs = "i", cex.main = 1.8,
     xlab = "", ylab = "", cex.axis = 1.4) 
title(xlab = "t", ylab = "", 
      line = 2, cex.lab = 1.5)
lines(gridPoints, beta2.1(gridPoints), lwd = 2, lty = 1, col = "grey52")
lines(gridPoints, beta21CB[,1], lwd = 2, lty = 2)
lines(gridPoints, beta21CB[,2], lwd = 2, lty = 2)

# Plot beta0.2(t)
beta02CB <- Inferenceout$beta2.CB[,c(1,2)]
plot(gridPoints, CARout$beta2Est[,1], 'l', col = "black", lwd = 2, 
     main = expression(paste("(d)"," ", widehat(beta)[0]^{(2)}, "(t)")),
     ylim = range(beta02CB), xlim = c(0,1), xaxs = "i", cex.main = 1.8,
     xlab = "", ylab = "", cex.axis = 1.4) 
title(xlab = "t", ylab = "", 
      line = 2, cex.lab = 1.5)
lines(gridPoints, beta0.2(gridPoints), lwd = 2, lty = 1, col = "grey52")
lines(gridPoints, beta02CB[,1], lwd = 2, lty = 2)
lines(gridPoints, beta02CB[,2], lwd = 2, lty = 2)

# Plot beta1.2(t)
beta12CB <- Inferenceout$beta2.CB[,c(3,4)]
plot(gridPoints, CARout$beta2Est[,2], 'l', col = "black", lwd = 2, 
     main = expression(paste("(e)"," ", widehat(beta)[1]^{(2)}, "(t)")),
     ylim = range(beta12CB), xlim = c(0,1), xaxs = "i", cex.main = 1.8,
     xlab = "", ylab = "", cex.axis = 1.4) 
title(xlab = "t", ylab = "", 
      line = 2, cex.lab = 1.5)
lines(gridPoints, beta1.2(gridPoints), lwd = 2, lty = 1, col = "grey52")
lines(gridPoints, beta12CB[,1], lwd = 2, lty = 2)
lines(gridPoints, beta12CB[,2], lwd = 2, lty = 2)

# Plot beta2.2(t)
beta22CB <- Inferenceout$beta2.CB[,c(5,6)]
plot(gridPoints, CARout$beta2Est[,3], 'l', col = "black", lwd = 2, 
     main = expression(paste("(f)"," ", widehat(beta)[2]^{(2)}, "(t)")),
     ylim = range(beta22CB), xlim = c(0,1), xaxs = "i", cex.main = 1.8,
     xlab = "", ylab = "", cex.axis = 1.4) 
title(xlab = "t", ylab = "", line = 2, cex.lab = 1.5)
lines(gridPoints, beta2.2(gridPoints), lwd = 2, lty = 1, col = "grey52")
lines(gridPoints, beta22CB[,1], lwd = 2, lty = 2)
lines(gridPoints, beta22CB[,2], lwd = 2, lty = 2)


###########################################################
# Plot estimates of multivariate eigenfunctions
# in two dimensions
###########################################################

# Estimated multivariate eigenfunctions
psi1Est <- MCARout$psi11Est
psi2Est <- MCARout$psi22Est


## Align eigenfunctions
# Note: the estimated eigenfunctions may flip signs. Here we match the sign of the
#      estimated eigenfunctions to the true eigenfunctions for plotting  

# Construct a function to calculate the mean squared deviation errors of the eigenfunction
# with the 1st output as the error value and the 2nd output as the indicator 
# with 1 indicating the sign of the eigenfunction has not been flipped, and 2 indicating
# the sign of the eigenfunction has been flipped.

MSDE.eifun <- function(x,y){
  err1 <- trapz(gridPoints, (x-y)^2)
  err2 <- trapz(gridPoints, (x+y)^2)
  return(c(min(err1, err2), which.min(c(err1, err2))))
}

# Calculate the mean squared deviation errors using the constructed `MSDE.eifun` function
# 1st dimension
MSDE.psi11 <- MSDE.eifun(psi1Est[,1], psi1.True[,1])
MSDE.psi12 <- MSDE.eifun(psi1Est[,2], psi1.True[,2]) 
MSDE.psi13 <- MSDE.eifun(psi1Est[,3], psi1.True[,3])
# 2nd dimension
MSDE.psi21 <- MSDE.eifun(psi2Est[,1], psi2.True[,1])
MSDE.psi22 <- MSDE.eifun(psi2Est[,2], psi2.True[,2]) 
MSDE.psi23 <- MSDE.eifun(psi2Est[,3], psi2.True[,3])

# Align 
# 1st dimension
fun1.sign.1 <- MSDE.psi11[2] * (-2) + 3 
# Note: If the sign of the eigenfunction has been flipped (indicator = 2), then *(-2)+3 would produce (-1) to give us back the 
# original sign; If the sign of the eigenfunction has not been flipped (indicator = 1), then *(-2)+3 would produce (1) 
# which is the original sign
psi1Est[,1] <- psi1Est[,1] * fun1.sign.1
fun1.sign.2 <- MSDE.psi12[2] * (-2) + 3
psi1Est[,2] <- psi1Est[,2] * fun1.sign.2
fun1.sign.3 <- MSDE.psi13[2] * (-2) + 3
psi1Est[,3] <- psi1Est[,3] * fun1.sign.3

# 2nd dimension
fun2.sign.1 <- MSDE.psi21[2] * (-2) + 3
psi2Est[,1] <- psi2Est[,1] * fun2.sign.1
fun2.sign.2 <- MSDE.psi22[2] * (-2) + 3
psi2Est[,2] <- psi2Est[,2] * fun2.sign.2
fun2.sign.3 <- MSDE.psi23[2] * (-2) + 3
psi2Est[,3] <- psi2Est[,3] * fun2.sign.3

# Plot estimates of multivariate eigenfunctions
# 1st dimension eigenfunctions
par(mfrow = c(2,2))
matplot(gridPoints, psi1.True, "l", ylim = c(-2,2), xaxs = "i", 
        main = "(a) True: 1st dimension", cex.main = 2, 
        xlab = "", ylab = "", lwd = 2, col = c("black", "red", "blue"))
legend("bottomleft", legend = c("1", "2", "3"), col = c("black", "red", "blue"), 
       lty = 1:3, cex = 1)
title(xlab = "time", ylab = expression(paste(psi^(1),(t))), 
      line = 2, cex.lab = 1.6)

matplot(gridPoints, psi1Est[,1:M], "l", ylim = c(-2,2), xaxs = "i", 
        main = "(b) Estimated: 1st dimension", cex.main = 2,
        xlab = "", ylab = "", lwd = 2, col = c("black", "red", "blue"))
title(xlab = "time",ylab = expression(paste(widehat(psi)^(1),(t))), 
      line = 2, cex.lab = 1.6)

# second dimension eigenfunctions
matplot(gridPoints, psi2.True, "l", ylim = c(-2,2), xaxs = "i", 
        main = "(c) True: 2nd dimension", cex.main = 2,
        xlab = "", ylab = "", lwd = 2, col = c("black", "red", "blue"))
legend("bottomleft", legend = c("1", "2", "3"), col = c("black", "red", "blue"), 
       lty = 1:3, cex = 1)
title(xlab = "time",ylab = expression(paste(psi^(2),(t))), 
      line = 2, cex.lab = 1.6)

matplot(gridPoints, psi2Est[,1:M], "l", ylim = c(-2,2), xaxs = "i", 
        main = "(d) Estimated: 2nd dimension", cex.main = 2,
        xlab = "", ylab = "", lwd = 2, col = c("black", "red", "blue"))
title(xlab = "time", ylab = expression(paste(widehat(psi)^(2),(t))), 
      line = 2, cex.lab = 1.6)


##############################################################
# Plot estimates of multivariate trajectories 
# in two dimensions
##############################################################

# Region-specific predicted trajectory for the first region
# Region id for plots
R.ID <- 1
# True trajectory for the `R.ID`'s region
R.traj1.T <- matrix(data.T$df.True$y1.True, nrow = nregion, byrow = TRUE)[R.ID,] #1st dimension
R.traj2.T <- matrix(data.T$df.True$y2.True, nrow = nregion, byrow = TRUE)[R.ID,] #2nd dimension
# Predicted trajectory for the `R.ID`'s region
R.traj1.Est <- Inferenceout$Y1.hat[R.ID,] #1st dimension
R.traj2.Est <- Inferenceout$Y2.hat[R.ID,] #2nd dimension

par(mfrow = c(2,1))
ylim1 <- min(R.traj1.T, R.traj1.Est)
ylim2 <- max(R.traj1.T, R.traj2.Est)
plot(gridPoints, R.traj1.T, "l", lwd = 2, ylim = c(ylim1,ylim2), xaxs = "i", 
     main = "1st dimension", cex.main = 2, xlab = "", ylab = "", col = "black")
title(xlab = "Time", ylab = "Trajectory", line = 2, cex.lab = 1.6)
lines(gridPoints, R.traj1.Est, col = "grey52", lwd = 2, lty = 1)
legend("topright", legend = c("True", "Predicted"), col = c("black", "grey52"), 
       lty = c(1,1), cex = 1)

ylim1 <- min(R.traj2.T, R.traj2.Est)
ylim2 <- max(R.traj2.T, R.traj2.Est)
plot(gridPoints, R.traj2.T, "l", lwd = 2, ylim = c(ylim1,ylim2), xaxs = "i",
     main = "2nd dimension", cex.main = 2, xlab = "", ylab = "", col = "black")
title(xlab = "Time",ylab = "Trajectory", line = 2, cex.lab = 1.6)
lines(gridPoints, R.traj2.Est, col = "grey52", lwd = 2, lty = 1)



