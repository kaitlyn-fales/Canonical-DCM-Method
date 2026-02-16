remove(list=ls())

library(R.matlab)
library(deSolve)
library(dplyr)

load("HCP_Social_Task/SPM/Results/SPM_results.RData")

########### Set up ############################################################
subjects <- c(100307,100408,101107,101309,101915,103111,103414,103818,105014,
              105115,106016,108828,110411,111312,111716,113619,113922,114419,
              115320,116524,117122,118528,118730,118932,120111,122317,122620,
              123117,123925,124422,125525,126325,127630,127933,128127,128632,
              129028,130013,130316,131217,131722,133019,133928,135225,135932,
              136833,138534,139637,140925,144832,146432,147737,148335,148840,
              149337,149539,149741,151223,151526,151627,153025,154734,156637,
              159340,160123,161731,162733,163129,176542,178950,188347,189450,
              190031,192540,196750,198451,199655,201111,208226,211417,211720,
              212318,214423,221319,239944,245333,280739,298051,366446,397760,
              414229,499566,654754,672756,751348,756055,792564,856766,857263,
              899885)

conditions <- c('phaseLR_maskL','phaseLR_maskR','phaseRL_maskL','phaseRL_maskR')

# 2 nodes, 2 experimental inputs
m = 2; n_u = 2

# Indices of parameters
A_idxs <- matrix(c(1,1,
                   2,1,
                   1,2,
                   2,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,1,1,
                   2,2,1,
                   2,1,2,
                   2,2,2), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1,
                   2,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

struct_paramMats = function(m, n_u, idxs, nu){
  
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
  
}

# ODE for neural activation that needs to be solved
linear <- function(t, z, params, input) {
  # t is the current time point in the integration, 
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]]
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+",lapply(1:nrow(u), function(i) u[i,1]*B[[i]]))
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}

# Putting together to form U matrix of experimental inputs
input_u = function(t){
  c(approxfun(x = times,y = u[,1],rule = 2)(t),
    approxfun(x = times,y = u[,2],rule = 2)(t))
}

################################################################################

# Empty list to store all subject achieved SNR
SNR_list <- list()

# Seed for reproducibility
set.seed(1234)

# Loop through and generate synthetic dataset
for (i in 1:length(subjects)){
  
  sub <- subjects[i]
  
  # Empty dataframe to store achieved SNR
  subject_SNR <- matrix(NA, nrow = length(conditions), ncol = 2)
  
  # Empty vectors to store final scaling (if needed)
  final_scale_flag <- numeric()
  final_scale_vec <- numeric()
  
  for (j in 1:length(conditions)){
    
    # Unique seed to determine noise
    set.seed(1234 + i*100 + j)
    
    # Load data for u matrix
    load(paste0("HCP_Social_Task/Data/sub-",sub,"_",conditions[j],".RData"))
    u <- rbind(0,dat$u)
    times <- c(0,dat$times)
    
    # Posterior means
    SPM <- SPM_results[[conditions[j]]]
    
    # Posterior means - MCMC
    nu_A = SPM$mean_alpha[1:4]
    nu_B = SPM$mean_alpha[5:8]
    nu_C = SPM$mean_alpha[9:10]
    
    nu <- list(nu_A, nu_B, nu_C)
    
    # Model params
    paramMats = struct_paramMats(m = m,n_u = n_u,idxs = idxs,nu = nu)
    
    # Solve the ODE for z
    out_z <- ode(y = rep(0.1,m),
                 times = times,
                 func = linear,
                 parms= with(paramMats, list(A=A, B=B, C=C)),
                 input = input_u,
                 atol = 1e-6, 
                 rtol = 1e-6
    )
    
    # Get rid of first column of z (same as times)
    out_z <- out_z[,-1]
    
    # hrf function (using the difference of two gammas - canonical)
    HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
    HRF_mu = function(mu,tp){
      convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
    }
    
    BOLD = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
    
    SNR <- 1.679376 # target (real SNR will be larger than this to avoid triggering SPM internal scaling)
    y_obs <- matrix(NA, nrow = length(times) - 1, ncol = m)
    variance <- numeric()
    
    # Add Gaussian noise but ensure range <= 4
    R_signal <- max(BOLD) - min(BOLD)
    max_allowed_noise_span <- 4 - R_signal
    if (max_allowed_noise_span <= 0)
      stop("Clean signal already has range >= 4: cannot add noise without triggering scaling.")
    
    for (k in 1:m) {
      # target variance for given SNR
      variance[k] <- (var(BOLD[, k]) + (mean(BOLD[, k]))^2) / SNR
      
      # draw Gaussian noise
      noise <- rnorm(length(times) - 1, mean = 0, sd = sqrt(variance[k]))
      
      # compute span of noise and scale it if it would exceed allowed range
      noise_span <- max(noise) - min(noise)
      if (noise_span > max_allowed_noise_span) {
        scale_factor <- max_allowed_noise_span / noise_span
        noise <- noise * scale_factor
      }
      
      # add noise to clean signal
      y_obs[, k] <- BOLD[, k] + noise
    }
    
    # Compute achieved SNR (mean across columns) 
    # SNR = var(signal) / var(noise)
    achieved_SNR <- numeric(m)
    for (k in 1:m) {
      achieved_SNR[k] <- var(BOLD[, k]) / var(y_obs[, k] - BOLD[, k])
    }
    
    # Final safety scaling like SPM, with flag to indicate when needed
    scale <- max(y_obs) - min(y_obs)
    final_scale_factor <- 4 / max(scale, 4)
    y_obs <- y_obs * final_scale_factor
    
    final_scale_flag[j] <- ifelse(final_scale_factor < 1, 1, 0)
    final_scale_vec[j] <- final_scale_factor
    
    # Save SNR results
    subject_SNR[j,] <- achieved_SNR
    
    # Save/export data
    dat <- list(times = times[-1], u = u[-1,], y_obs = y_obs)
    save(dat, file = paste0("HCP_Social_Task/Analysis/Sensitivity_Sim/SPMpar_MCMCdgm/Data/sub-",sub,"_",conditions[j],".RData"))
    
    # Write observed signal to a MATLAB file for SPM 
    writeMat(paste0("HCP_Social_Task/Analysis/Sensitivity_Sim/SPMpar_MCMCdgm/Data/sub_",sub,"_",conditions[j],".mat"),
             U = dat$u,
             Y = dat$y_obs)
    
  }
  
  # Add subject ID and conditions names to SNR df
  subject_SNR <- data.frame(subject_SNR,final_scale_flag,final_scale_vec)
  subject_SNR$subject <- sub
  subject_SNR$condition <- conditions
  
  SNR_list[[i]] <- subject_SNR
  
}

# Save SNR_list
SNR_df <- bind_rows(SNR_list)
colnames(SNR_df)[1:2] <- c("V5_SNR","pSTS_SNR")
save(SNR_df, file = "HCP_Social_Task/Analysis/Sensitivity_Sim/SPMpar_MCMCdgm/Data/sim_dat_SNR.RData")

