########### Preprocess data for social task DCM ##################

setwd("/storage/work/krf5429/HCP_Social_Task")

# Packages
library(RNifti)
library(nlme)

# Directory for results
data_dir <- "Data"

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
phase <- "LR"

# Type of mask (right or left)
mask_type <- "L"
mask <- readNifti(paste0("mask_",mask_type,".nii"))

# Base directory for raw data
base_dir <- "/scratch/krf5429/Data"

# Scanning specifications
nscan <- 274
TR <- 0.72

# Time
max_time = TR*nscan
times <- seq(0,max_time,TR)

# Creation of stimulus indicator function
stim.indicator <- function(os,dur){
  stim <- rep(0,length(times)-1)
  for (i in 1:length(os)){
    stim[os[i]:(os[i] + dur[i] - 1)] <- 1
  }
  return(stim)
}

########## Loop to process subject data ################################

for (i in 1:length(subjects)){
  # Subject to loop through
  sub <- subjects[i]
  
  ########### Generate stimulus indicator for social task ########### 
  # Load file for task
  os_mental <- read.delim(paste0(base_dir,"/",sub,"/tfMRI_SOCIAL_",phase,"/EVs/mental.txt"),
                          header = F)
  os_rnd <- read.delim(paste0(base_dir,"/",sub,"/tfMRI_SOCIAL_",phase,"/EVs/rnd.txt"),
                       header = F)
  
  # Stimulus onset times (in scans)
  motion <- ceiling(sort(c(os_mental$V1,os_rnd$V1))/TR)
  animate <- ceiling(sort(c(os_mental$V1))/TR)
  
  # Duration times (in scans)
  dur1 <- ceiling(c(os_mental$V2,os_rnd$V2)/TR)
  dur2 <- ceiling(os_mental$V2/TR)
  
  # Running stim indicator function for each stimulus
  stim1 <- stim.indicator(motion,dur1)
  stim2 <- stim.indicator(animate,dur2)
  
  # Combine together
  u <- cbind(stim1,stim2)
  colnames(u) <- c("All Motion","Animate")
  ########################################################################
  
  ########### Prep data by pulling first PC from each VOI ################
  dat <- readNifti(paste0(base_dir,"/",sub,"/tfMRI_SOCIAL_",phase,"/tfMRI_SOCIAL_",
                          phase,"_hp0_clean_rclean_tclean.nii.gz"))
  y_obs <- matrix(NA, nrow = dim(dat)[4], ncol = 2)
  colnames(y_obs) <- c(paste0(mask_type,"_V5"),paste0(mask_type,"_pSTS"))
  
  for (i in 1:2){
    Vxyz <- t(which(mask == i, arr.ind=TRUE))
    BOLD <- matrix(data = NA, nrow = dim(dat)[4], ncol = ncol(Vxyz))
    for (j in 1:ncol(BOLD)){
      BOLD[,j] <- dat[Vxyz[1,j],Vxyz[2,j],Vxyz[3,j],]
    }
    sds <- apply(BOLD, 2, sd)
    BOLD <- BOLD[, sds > 0]
    VOI <- prcomp(BOLD, center = T, scale. = T)$x[,1]
    y_obs[,i] <- VOI
  }
  
  # Enforce scaling like SPM
  scale <- max(y_obs) - min(y_obs)
  scale <- 4 / max(scale, 4)
  y_obs <- y_obs * scale
  
  # Plot data and export
  png(paste0(data_dir,"/Figures/sub-",sub,"_phase",phase,"_mask",mask_type,".png"), 
      width = 800, height = 600, type = "cairo")
  par(mfrow = c(1,2))
  for (i in 1:2) plot(y_obs[,i], type = "l", 
                      xlab = "Scans", 
                      ylab = colnames(y_obs)[i],
                      ylim = c(min(y_obs),max(y_obs)))
  mtext(paste0("Subject ",sub," - ",mask_type," Hemisphere Data for ",phase," Phase Encoding"),
        line = -2, outer = T)
  dev.off()
  ########################################################################
  
  # Export data file
  dat <- list(times = times[-1], u = u, y_obs = y_obs, scale = scale)
  save(dat, file = paste0(data_dir,"/sub-",sub,"_phase",phase,"_mask",mask_type,".RData"))
}

########################################################################





