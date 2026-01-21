# Data preparation for SPM - convert .RData files in Data folder to .mat files

setwd("/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM")

library(R.matlab)

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
############################################################################

# Phase-encoding (run/session)
phase <- c("LR","RL")

# Type of mask (right or left)
mask_type <- c("L","R")

# Combinations of all phases and mask types
combos <- expand.grid(phase = phase,mask = mask_type)

# For loop to load .RData data files and write as .mat files for SPM
for (i in 1:length(subjects)){
  
  for (j in 1:nrow(combos)){
    
    # Load .RData file
    load(paste0("../Data/sub-",subjects[i],"_phase",combos$phase[j],"_mask",combos$mask[j],".RData"))
    
    # Write observed signal to a MATLAB file for SPM
    writeMat(paste0("Data/sub_",subjects[i],"_phase",combos$phase[j],"_mask",combos$mask[j],".mat"),
             U = dat$u,
             Y = dat$y_obs)
  }
  
}







