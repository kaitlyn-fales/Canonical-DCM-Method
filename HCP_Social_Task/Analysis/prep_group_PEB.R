# PEB for combining subjects into group DCM
library(R.matlab)

rm(list = ls())

# Load in diagnostics - only include converged chains
load("HCP_Social_Task/Analysis/MCMC_diagnostics.RData")

# Hypothesis
A_idxs <- matrix(c(2,1,
                   1,2,
                   1,1,
                   2,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,2,1,
                   2,1,2,
                   2,1,1,
                   2,2,2), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1,
                   2,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

# Number of regions and inputs
m <- 2; n_u  <- 2

# Function for structuring raw parameters (diagonal not reparameterized)
struct_raw_paramMats = function(m, n_u, idxs, nu){
  
  A = matrix(data = 0,nrow = m, ncol = m)
  for(i in 1:nrow(idxs$A_idxs)){
    A[idxs$A_idxs[i,1], idxs$A_idxs[i,2]] = with(nu,nu_A[i])
  }
  
  B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
  for(i in 1:nrow(idxs$B_idxs)){
    B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
  }
  
  C = matrix(data=0, nrow = m, ncol = n_u)
  for(i in 1:nrow(idxs$C_idxs)){
    C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
  }
  
  return(list(A=A,B=B,C=C))
  
} # diagonal not modified for SPM - need raw values

# Compile stan program from parent directory
canonical_dcm = cmdstanr::cmdstan_model("canonical_dcm.stan")

# Source functions from parent directory
source("canonical_dcm_functions.R")

# Set up MCMC posterior outputs for SPM PEB
###############################################

for (i in 1:nrow(diagnostics)){
  
    # Load your R posterior draws
    load(paste0("HCP_Social_Task/Output/sub-", diagnostics$subject[i],
                "_", diagnostics$condition[i], "_draws.RData"))
    
    # Load your data
    load(paste0("HCP_Social_Task/Data/sub-", diagnostics$subject[i],
                "_", diagnostics$condition[i], ".RData"))
    
    # Summarise posterior draws
    summary <- posterior::summarise_draws(draws_df)
    
    ############## Do posterior means for nu #########################
    # Do posterior means for nu
    nu = list(nu_A = summary$mean[4:7],
              nu_B = summary$mean[8:11], 
              nu_C = summary$mean[12:13])
    
    paramMats <- struct_raw_paramMats(m = m, n_u = n_u, idxs = idxs, nu = nu)
    ##################################################################
    
    ############## Do posterior covariance ###########################
    draws_mat <- as_draws_matrix(draws_df)[,4:13]
    
    # Change order of columns of draws matrix to match the order of SPM
    SPM_indices <- c(3,1,2,4,7,5,6,8,9,10)
    draws_mat <- draws_mat[,SPM_indices]
    
    # Add zero columns to pad matrix for the params not est in hypothesis to match SPM
    zero_pad4 <- matrix(0, ncol = 4, nrow = nrow(draws_mat))
    zero_pad2 <- matrix(0, ncol = 2, nrow = nrow(draws_mat))
    
    # Adding 4 zero col for B[1] since there is no B param for input1
    # Adding 2 zero col for final two param in C since only two estimated
    draws_mat <- cbind(draws_mat[,1:4],zero_pad4,draws_mat[,5:10],zero_pad2)
    
    Cp <- round(cov(draws_mat), digits = 8)
    rownames(Cp) <- NULL
    colnames(Cp) <- NULL
    #################################################################
    
    ############## Use VL to get approx free energy (ELBO) ##########
    # Put data into form for sampler
    stan_dat <- get_stan_dat(dat, idxs)
    
    # Initial values are posterior means
    inits <- list(lp__ = c(summary$mean[1]),
                  sigma = summary$mean[2:3],
                  nu_A = summary$mean[4:7],
                  nu_B = summary$mean[8:11], 
                  nu_C = summary$mean[12:13],
                  z0 = summary$mean[14:15],
                  beta = summary$mean[16:17])
    
    pf_out <- capture.output({canonical_dcm$pathfinder(data = stan_dat, 
                                                       init = list(inits),
                                                       seed = 1234,
                                                       num_paths = 1)
                             })
    
    # Find the line with 'Best Iter'
    best_line <- pf_out[grep("Best Iter:", pf_out)]
    
    # Extract the number after 'ELBO ('
    elbo_val <- str_extract(best_line, "(?<=ELBO \\().*?(?=\\))")
    
    # Convert to numeric
    elbo <- as.numeric(elbo_val)
    
    #################################################################
    
    # Save as .mat file
    writeMat(paste0("HCP_Social_Task/Output/sub_", diagnostics$subject[i],
                    "_", diagnostics$condition[i], ".mat"),
             A = paramMats$A,
             B1 = paramMats$B[[1]],
             B2 = paramMats$B[[2]],
             C = paramMats$C,
             Cp = Cp,
             DCM_F = elbo)
}

