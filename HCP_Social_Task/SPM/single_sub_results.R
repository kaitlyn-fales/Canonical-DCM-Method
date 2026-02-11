# Examine single subject differences in output parameters
library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(Cairo)

# Load diagnostics df
load("HCP_Social_Task/Analysis/diagnostics_compilation.RData")

# Summary by subject, then sample proportionally from the subjects by group
set.seed(1234)
subject_sample <- diagnostics_df %>%
  distinct(subject) %>%
  slice_sample(n = 10)

# Phase-encoding (run/session)
phase <- c("LR","RL")

# Type of mask (right or left)
mask_type <- c("L","R")

# Combinations of all phases and mask types
combos <- expand.grid(phase = phase,mask = mask_type)

############ Functions ##############
# Unstructure params function
unstruct_paramMats = function(params, idxs){
  nu_A = params$A[idxs$A_idxs]
  nu_C = params$C[idxs$C_idxs]
  nu_B = c()
  for(i in 1:nrow(idxs$B_idxs)){
    nu_B[i] = params$B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]]
  }
  c(nu_A = nu_A, nu_B = nu_B, nu_C = nu_C)
}

# Function for summary
get_summary <- function(mu,sigma2){
  
  # Set alpha level
  alpha = 0.1
  
  intervals <- matrix(NA, nrow=length(mu), ncol=2)
  for (i in 1:nrow(intervals)){
    q_lower <- qnorm(alpha/2, mean = mu[i], sd = sqrt(sigma2[i]))
    q_upper <- qnorm(1-alpha/2, mean = mu[i], sd = sqrt(sigma2[i]))
    intervals[i,] <- c(q_lower, q_upper)
  }
  
  df <- data.frame(mu,intervals)
  colnames(df) <- c("mean","q5","q95")
  
  df$CI_length <- abs(df$q95-df$q5)
  
  return(df)
}

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
#####################################

# For loop to run through subjects
for (i in 1:nrow(subject_sample)){
  # Choose subject - run all code below for each subject in dataframe
  sub <- subject_sample$subject[i]
  
  summaries <- list()
  for (j in 1:nrow(combos)){
    DCM <- readMat(paste0("HCP_Social_Task/SPM/Output/DCM_sub_", sub,
                          "_phase", combos$phase[j], "_mask", combos$mask[j], "_out.mat"))
  
    # Posterior means
    means <- list(A = DCM$Ep.A, B = DCM$Ep.B, C = DCM$Ep.C)
    
    # Grab posterior means
    means$B <- list(means$B[,,1],means$B[,,2])
    
    # Unstructure means
    mu <- unstruct_paramMats(means,idxs)
    
    # All posterior variances coming from A, B, C (hemodynamic params are last 4)
    Cp <- diag(DCM$Cp)
    Cp <- Cp[1:(length(Cp)-4)]
    
    # Get rid of nonzero entries - correspond to placeholders not est in A, B, C
    Cp <- Cp[Cp != 0]
    
    # Get results
    summary <- get_summary(mu,Cp)
    
    # Add names in 
    variable <- c("nu_A[1]","nu_A[2]","nu_A[3]","nu_A[4]","nu_B[1]",
                   "nu_B[2]","nu_B[3]","nu_B[4]","nu_C[1]","nu_C[2]")
    
    # Make data frame
    summary <- data.frame(variable,summary)
    rownames(summary) <- NULL
    
    summary$phase_mask <- paste(combos$phase[j],combos$mask[j],sep = "_")
    summaries[[j]] <- summary
  }
  
  df <- bind_rows(summaries)
  
  df <- df %>%
    mutate(
      param_group = str_extract(variable, "nu_[ABC]"),
      param_index = str_extract(variable, "(?<=\\[)\\d+(?=\\])")
    )
  
  # Make sure the order is preserved nicely
  df$param_index <- factor(df$param_index, levels = sort(unique(df$param_index)))
  df$phase_mask <- factor(df$phase_mask, levels = c("LR_L", "LR_R", "RL_L", "RL_R"))
  
  # Now plot for all param_group 
  ggplot(df,
         aes(x = param_index, y = mean, color = phase_mask)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_vline(xintercept = as.numeric(factor(unique(df$param_index))) + 0.5,
               linetype = "dashed", color = "gray40") +
    geom_point(position = position_dodge(width = 0.8), size = 2.5) +
    geom_errorbar(aes(ymin = q5, ymax = q95),
                  position = position_dodge(width = 0.8), width = 0.8) +
    facet_wrap(~param_group, scales = "free_y", labeller = label_parsed) +
    labs(x = "Parameter index", y = "Estimate", color = "Phase mask") +
    ggtitle(paste0("Subject ",sub)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  dims <- dev.size("in") 
  ggsave(paste0("HCP_Social_Task/SPM/Results/sub_",sub,"_post_means.png"))
  
  ##########################################################################
  
  CairoPNG(paste0("HCP_Social_Task/SPM/Results/sub_", sub, "_est_vs_obs.png"),
           width = dims[1], height = dims[2], units = "in", dpi = 96)
  par(mfrow = c(2,2), mar = c(4, 4, 4, 1), oma = c(2, 1, 2, 1))
  
  # Loop through combos, plot est vs obs and export
  for (j in 1:nrow(combos)){
    load(paste0("HCP_Social_Task/Data/sub-",sub,"_phase",combos$phase[j],"_mask",combos$mask[j],".RData"))
    times <- c(0,dat$times)
    y_obs <- dat$y_obs
    
    # Load predicted signal
    DCM <- readMat(paste0("HCP_Social_Task/SPM/Output/DCM_sub_", sub,
                          "_phase", combos$phase[j], "_mask", combos$mask[j], "_out.mat"))
    y_pred <- DCM$DCM.y
    
    # Compute SPM predicted shift (beta) - not available from program directly
    beta_offset <- colMeans(y_obs - y_pred)
    y_obs <- sweep(y_obs, 2, beta_offset, "-")
    
    # Choose colors for each region
    cols_all <- brewer.pal(8, "Dark2")
    cols <- cols_all[c(1, 3)]
    
    # ylim range
    y_range <- range(c(y_pred, y_obs))
    ylim_exp <- c(y_range[1], y_range[2] * 1.15)
    
    # Set up plot
    plot(times[-1], y_obs[,1], type = "n", 
         xlab = "Time (s)", ylab = "BOLD Response",
         ylim = ylim_exp,
         main = paste0("Phase ", combos$phase[j],
                       ", Mask ", combos$mask[j]))
    
    # Add observed (dashed) and predicted (solid) lines by region
    for (i in 1:ncol(y_pred)) {
      lines(times[-1], y_obs[,i], col = cols[i], lty = 3, lwd = 1.5)   # dashed
      lines(times[-1], y_pred[,i], col = cols[i], lty = 1, lwd = 1.5)  # solid
    }
    
    # Add a legend
    legend("top", legend = c("V5","p-STS"), horiz = T,
           col = cols, lty = 1, lwd = 2, bty = "n")
    
    # Assign based on phase and mask
    output <- list(y_obs = y_obs, y_pred = y_pred)
    
    assign(paste0("phase",combos$phase[j],"_mask",combos$mask[j]), output)
  }
  
  # Add Subject title to plot
  mtext(paste0("Predicted (Solid) vs. Observed (Dotted): Subject ",sub), line = 0, outer = T)
  
  dev.off()
  
  # Save estimated signal to compare with SPM later
  SPM_result  <- list(phaseLR_maskL = phaseLR_maskL,
                      phaseLR_maskR = phaseLR_maskR,
                      phaseRL_maskL = phaseRL_maskL,
                      phaseRL_maskR = phaseRL_maskR)
  
  save(SPM_result, file = paste0("HCP_Social_Task/SPM/Results/sub_",sub,"_pred_vs_obs.RData"))
  
}

