remove(list=ls())

library(R.matlab)
library(dplyr)

# Scanning specifications
nscan <- 274
TR <- 0.72

# Time
max_time = TR*nscan
times <- seq(0,max_time,TR)

################ List of subject IDs to loop through #######################
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
############################################################################

# Empty list to store all subject achieved SNR
SNR_list <- list()

# Seed for reproducibility
set.seed(1234)

# Loop through all subject and conditions
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
    
    # Load in SPM-generated signal
    SPM <- readMat(paste0("HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/Data/sub_",sub,"_",conditions[j],".mat"))
    BOLD <- SPM$y.signal
    u <- SPM$U
    m <- ncol(BOLD)
    
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
    dat <- list(times = times[-1], u = u, y_obs = y_obs)
    save(dat, file = paste0("HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/Data/sub-",sub,"_",conditions[j],".RData"))
    
    # Write observed signal to a MATLAB file for SPM (overwrite the signal file)
    writeMat(paste0("HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/Data/sub_",sub,"_",conditions[j],".mat"),
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
save(SNR_df, file = "HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/Data/sim_dat_SNR.RData")

