# MSE calculation for predicted signals vs true signals between MCMC and SPM

library(R.matlab)

# Function for calculating MSE
calculate_MSE <- function(SPM, MCMC){
  # Calculate MSE for SPM
  mse_spm <- mean((truth - SPM)^2)
  
  # Loop through MCMC results and calculate MSE for each seed
  mse_mcmc <- numeric()
  for (j in 1:length(MCMC)){
    mse_mcmc[j] <- mean((truth - MCMC[[j]])^2)
  }
  
  # Calculate mean MSE and SE
  mean_mcmc <- mean(mse_mcmc)
  se_mcmc <- sd(mse_mcmc)/sqrt(length(MCMC))
  
  result <- round(c(mean_mcmc,se_mcmc,mse_spm), digits = 4)
  names(result) <- c("MCMC Mean","MCMC SE","SPM")
  
  return(result)
}

########## Balloon DGM ############
####### Diag A, zero delays ###########
# Load in true signal
truth <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagA_zero.mat")
truth <- truth$y.signal

# Load in SPM predicted signal
SPM <- readMat("Balloon_Simulation/Output/balloon_diagA_zero_out.mat")
SPM <- SPM$DCM.y

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/balloon_diagA_zero_pred.RData")

# Calculate MSE
calculate_MSE(SPM,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

####### Diag AB, zero delays ###########
# Load in true signal
truth <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagAB_zero.mat")
truth <- truth$y.signal

# Load in SPM predicted signal
SPM <- readMat("Balloon_Simulation/Output/balloon_diagAB_zero_out.mat")
SPM <- SPM$DCM.y

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/balloon_diagAB_zero_pred.RData")

# Calculate MSE
calculate_MSE(SPM,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

####### Diag A, 1s delays ###########
# Load in true signal
truth <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagA_nonzero.mat")
truth <- truth$y.signal

# Load in SPM predicted signal
SPM <- readMat("Balloon_Simulation/Output/balloon_diagA_nonzero_out.mat")
SPM <- SPM$DCM.y

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/balloon_diagA_nonzero_pred.RData")

# Calculate MSE
calculate_MSE(SPM,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

####### Diag AB, 1s delays ###########
# Load in true signal
truth <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagAB_nonzero.mat")
truth <- truth$y.signal

# Load in SPM predicted signal
SPM <- readMat("Balloon_Simulation/Output/balloon_diagAB_nonzero_out.mat")
SPM <- SPM$DCM.y

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/balloon_diagAB_nonzero_pred.RData")

# Calculate MSE
calculate_MSE(SPM,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

########## Canonical DGM ############
####### Diag A ###########
# Load in true signal
load("Balloon_Simulation/Data/canonical_sim_signal_diagA.RData")
truth <- true_signal

# Load in SPM predicted signal
SPM <- readMat("Balloon_Simulation/Output/canonical_diagA_out.mat")
SPM <- SPM$DCM.y

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/canonical_diagA_pred.RData")

# Calculate MSE
calculate_MSE(SPM,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

####### Diag AB ###########
# Load in true signal
load("Balloon_Simulation/Data/canonical_sim_signal_diagAB.RData")
truth <- true_signal

# Load in SPM predicted signal
SPM <- readMat("Balloon_Simulation/Output/canonical_diagAB_out.mat")
SPM <- SPM$DCM.y

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/canonical_diagAB_pred.RData")

# Calculate MSE
calculate_MSE(SPM,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################


