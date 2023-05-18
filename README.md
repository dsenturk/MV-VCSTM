CONTENTS OF THIS FOLDER ——————————————

MV_VCSTM_tutorial.R : A step-by-step implementation of MV-VCSTM and the associated procedures described in "Multivariate Varying Coefficient Spatiotemporal Model".

MV_VCSTM_simulation.R : Function for simulating one data set under the simulation design described in Section 4 of "Multivariate Varying Coefficient Spatiotemporal Model".

MV_VCSTM_univariate_decomposition.R : Function for fitting multivariate varying coefficient function model and employ univariate FPCA decomposition on the residuals (estimation steps 1-2 in the estimation algorithm) described in "Multivariate Varying Coefficient Spatiotemporal Model", including estimation of univariate eigenfunctions and PC scores.

MV_VCSTM_MCAR.R : Function for the 1st MCMC estimation (estimation step 3-4 in the estimation algorithm) of MV-VCSTM model described in "Multivariate Varying Coefficient Spatiotemporal Model", including estimation of within eigencomponent variation, as well as the multivariate eigenfunctions.

MV_VCSTM_CAR.R : Function for the 2nd MCMC estimation (estimation step 5 in the estimation algorithm) of MV-VCSTM model described in "Multivariate Varying Coefficient Spatiotemporal Model", including estimation of varying coefficient functions (VCFs), spatial variance parameters, measurement error variance and multivariate PC scores.

MV_VCSTM_inference.R : Function for constructing simultaneous credible bands for VCFs, as well as obtaining prediction for multivariate trajectories described in Section 2.2 of "Multivariate Varying Coefficient Spatiotemporal Model".

INTRODUCTION ——————————————

The contents of this folder allow for implementation of the MV-VCSTM estimation and inference described in "Multivariate Varying Coefficient Spatiotemporal Model". Users can simulate a sample data frame (MV_VCSTM_simulation.R) and apply the proposed estimation algorithm (MV_VCSTM_univariate_decomposition.R, MV_VCSTM_MCAR.R, MV_VCSTM_CAR.R). Also, we include tools to construct simultaneous credible bands for VCFs and perform predictions for multvariate trajectories (MV_VCSTM_inference.R). Detailed instructions on how to perform the aforementioned procedures, make inference on VCFs and predictions of multivariate trajectories, and visualize results are included in MV_VCSTM_tutorial.R.

REQUIREMENTS ——————————————

The included R programs require R 4.0.5 (R Core Team, 2021) and the packages listed in MV_VCSTM_tutorial.R, as well as pre-installed WinBUGS software. In order to properly call WinBUGS from R, the computer system should be Windows.

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using commands in MV_VCSTM_tutorial.R, as well as pre-install WinBUGS software.
