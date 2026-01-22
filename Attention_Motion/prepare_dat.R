######### Preparation of u matrix for SPM attention to motion data ###############

# Packages
library(R.matlab)

# Load in stimulus onset scans
stimuli <- readMat("Attention_Motion/SPM_VBA_Comparison/factors.mat")

# Getting stimuli according to model design (in scans)
photic.os <- sort(c(stimuli$att,stimuli$natt,stimuli$stat))
motion.os <- sort(c(stimuli$att,stimuli$natt)) 
attention.os <- c(stimuli$att)

# Setting duration (in scans)
dur <- 10

# Times
TR <- 3.22; nscan <- 360; max_time = TR*nscan
times <- seq(0,max_time,by = TR)

# Stimulus indicator function
stim.indicator <- function(os,dur){
  stim <- rep(0,length(times)-1)
  for (i in 1:length(os)){
    stim[os[i]:(os[i] + dur - 1)] <- 1
  }
  return(stim)
}

# Get stimuli
Photic <- stim.indicator(photic.os,dur)
Motion <- stim.indicator(motion.os,dur)
Attention <- stim.indicator(attention.os,dur)

# U matrix
u <- cbind(Photic,Motion,Attention)

##########################################
# Load in data
dat <- readMat("Attention_Motion/SPM_VBA_Comparison/y_obs_SPM.mat")
colnames(dat$y.obs) <- c("V1","V5","SPC")

motion_dat <- list(times = times[-1], u = u, y_obs = dat$y.obs)
save(motion_dat, file = "Attention_Motion/motion_dat.RData")






