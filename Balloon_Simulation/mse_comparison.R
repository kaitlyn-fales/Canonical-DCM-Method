# MSE calculation for predicted signals vs true signals between MCMC and SPM

library(R.matlab)

# Function for calculating MSE
calculate_MSE <- function(SPM, MCMC){
  # Calculate MSE for SPM
  mse_spm <- numeric()
  for (i in 1:length(SPM)){
    mse_spm[i] <- mean((truth - SPM[[i]])^2)
  }
  
  # Loop through MCMC results and calculate MSE 
  mse_mcmc <- numeric()
  for (j in 1:length(MCMC)){
    mse_mcmc[j] <- mean((truth - MCMC[[j]])^2)
  }
  
  # Calculate mean MSE and SE
  mean_mcmc <- mean(mse_mcmc)
  mean_spm <- mean(mse_spm)
  
  result <- round(c(mean_mcmc,mean_spm), digits = 4)
  names(result) <- c("CDCM Mean","SPM Mean")
  
  return(result)
}

########## Balloon DGM ############
####### Diag A, zero delays ###########
# Load in true signal
truth <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagA_zero.mat")
truth <- truth$y.signal

# Load in SPM predicted signals into a list
SPM_pred_signal <- list()
for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/balloon_diagA_zero_rep%03d_out.mat", k))
  SPM_pred_signal[[k]] <- result$DCM.y
}

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/balloon_diagA_zero_pred.RData")

# Calculate MSE
calculate_MSE(SPM_pred_signal,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

####### Diag AB, zero delays ###########
# Load in true signal
truth <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagAB_zero.mat")
truth <- truth$y.signal

# Load in SPM predicted signals into a list
SPM_pred_signal <- list()
for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/balloon_diagAB_zero_rep%03d_out.mat", k))
  SPM_pred_signal[[k]] <- result$DCM.y
}

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/balloon_diagAB_zero_pred.RData")

# Calculate MSE
calculate_MSE(SPM_pred_signal,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

####### Diag A, 1s delays ###########
# Load in true signal
truth <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagA_nonzero.mat")
truth <- truth$y.signal

# Load in SPM predicted signals into a list
SPM_pred_signal <- list()
for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/balloon_diagA_nonzero_rep%03d_out.mat", k))
  SPM_pred_signal[[k]] <- result$DCM.y
}

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/balloon_diagA_nonzero_pred.RData")

# Calculate MSE
calculate_MSE(SPM_pred_signal,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

####### Diag AB, 1s delays ###########
# Load in true signal
truth <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagAB_nonzero.mat")
truth <- truth$y.signal

# Load in SPM predicted signals into a list
SPM_pred_signal <- list()
for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/balloon_diagAB_nonzero_rep%03d_out.mat", k))
  SPM_pred_signal[[k]] <- result$DCM.y
}

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/balloon_diagAB_nonzero_pred.RData")

# Calculate MSE
calculate_MSE(SPM_pred_signal,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

########## Canonical DGM ############
####### Diag A ###########
# Load in true signal
load("Balloon_Simulation/Data/canonical_sim_signal_diagA.RData")
truth <- true_signal

# Load in SPM predicted signals into a list
SPM_pred_signal <- list()
for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/canonical_diagA_rep%03d_out.mat", k))
  SPM_pred_signal[[k]] <- result$DCM.y
}

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/canonical_diagA_pred.RData")

# Calculate MSE
calculate_MSE(SPM_pred_signal,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################

####### Diag AB ###########
# Load in true signal
load("Balloon_Simulation/Data/canonical_sim_signal_diagAB.RData")
truth <- true_signal

# Load in SPM predicted signals into a list
SPM_pred_signal <- list()
for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/canonical_diagAB_rep%03d_out.mat", k))
  SPM_pred_signal[[k]] <- result$DCM.y
}

# Load in MCMC predicted signal
load("Balloon_Simulation/Output/canonical_diagAB_pred.RData")

# Calculate MSE
calculate_MSE(SPM_pred_signal,pred_signal)

rm(list = setdiff(ls(), "calculate_MSE"))
######################################


