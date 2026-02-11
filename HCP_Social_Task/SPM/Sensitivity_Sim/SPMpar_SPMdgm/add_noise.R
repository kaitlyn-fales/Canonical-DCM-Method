remove(list=ls())

library(R.matlab)

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

for (i in 1:length(subjects)){
  
  sub <- subjects[i]
  
  for (j in 1:length(conditions)){
    
    # Load in SPM-generated signal
    SPM <- readMat(paste0("HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/Data/sub_",sub,"_",conditions[j],".mat"))
    BOLD <- SPM$y.signal
    u <- SPM$U
    m <- ncol(BOLD)
    
    # Add some white Gaussian noise according to SNR
    SNR = 1.679376
    variance <-  numeric()
    y_obs <- matrix(NA, nrow = length(times)-1, ncol = m)
    for (k in 1:m){
      variance[k] <- (var(BOLD[,k])+(mean(BOLD[,k]))^2)/SNR
      y_obs[,k] <- BOLD[,k] + rnorm(length(times)-1, mean = 0, sd = sqrt(variance[k]))
    }
    
    # Enforce scaling like SPM
    scale <- max(y_obs) - min(y_obs)
    scale <- 4 / max(scale, 4)
    y_obs <- y_obs * scale
    
    # Save/export data
    dat <- list(times = times[-1], u = u, y_obs = y_obs)
    save(dat, file = paste0("HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/Data/sub-",sub,"_",conditions[j],".RData"))
    
    # Write observed signal to a MATLAB file for SPM (overwrite the signal file)
    writeMat(paste0("HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/Data/sub_",sub,"_",conditions[j],".mat"),
             U = dat$u,
             Y = dat$y_obs)
    
  }
  
}



