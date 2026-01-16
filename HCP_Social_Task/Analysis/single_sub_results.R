# Examine single subject differences in output parameters
setwd("/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/Analysis")

library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(Cairo)

# Load diagnostics df after running diagnostics_EDA.R
load("diagnostics_compilation_final.RData")

# Summary by subject, then sample proportionally from the subjects by group
subject_summary <- diagnostics_df %>%
  group_by(subject) %>%
  summarise(prop_converged = mean(converged)) %>%
  arrange(desc(prop_converged))

set.seed(1234)
n_total <- 5

# Compute how many samples per group
group_counts <- subject_summary %>%
  count(prop_converged, name = "group_size") %>%
  mutate(
    n_to_sample = pmax(1, round(n_total * group_size / sum(group_size)))  # ensure at least 1
  )

# Apply group-wise sampling
sampled_subjects <- subject_summary %>%
  inner_join(group_counts, by = "prop_converged") %>%
  group_by(prop_converged) %>%
  group_modify(~ dplyr::slice_sample(.x, n = .x$n_to_sample[1])) %>%
  ungroup()

sampled_subjects

# Choose subject - run all code below for each subject in sampled_subjects dataframe
sub <- sampled_subjects$subject[1]

# Phase-encoding (run/session)
phase <- c("LR","RL")

# Type of mask (right or left)
mask_type <- c("L","R")

# Combinations of all phases and mask types
combos <- expand.grid(phase = phase,mask = mask_type)

summaries <- list()
for (j in 1:nrow(combos)){
  load(paste0("../Output/sub-",sub,"_phase",combos$phase[j],"_mask",combos$mask[j],"_draws.RData"))
  summary <- posterior::summarise_draws(draws_df)[4:15,c(1:2,6:7)]
  summary <- summary %>%
    mutate(across(where(is.numeric), ~round(., digits = 3)))
  summary$phase_mask <- paste(combos$phase[j],combos$mask[j],sep = "_")
  summaries[[j]] <- summary
}

df <- bind_rows(summaries)

df <- df %>%
  mutate(
    param_group = str_extract(variable, "nu_[ABC]|z0"),
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
ggsave(paste0("sub_",sub,"_post_means.png"))

########### Generate predicted signal from posterior means ###############

########### Set up experimental details ######################################
# 2 nodes, 2 experimental inputs
m = 2; n_u = 2

# Indices of parameters
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

struct_paramMats = function(m, n_u, idxs, nu){
  
  A = matrix(data = 0,nrow = m, ncol = m)
  for(i in 1:nrow(idxs$A_idxs)){
    A[idxs$A_idxs[i,1], idxs$A_idxs[i,2]] = with(nu,nu_A[i])
  }
  diag(A) = -0.5*exp(diag(A))
  
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

##########################################################################

CairoPNG(paste0("sub_", sub, "_est_vs_obs.png"),
         width = dims[1], height = dims[2], units = "in", dpi = 96)
par(mfrow = c(2,2), mar = c(3, 4, 4, 1), oma = c(1, 1, 2, 1))

# Loop through combos, plot est vs obs and export to .mat file for SPM
for (j in 1:nrow(combos)){
  load(paste0("../Data/sub-",sub,"_phase",combos$phase[j],"_mask",combos$mask[j],".RData"))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  y_obs <- dat$y_obs
  
  ######### Solve neural signal z #############################################
  # Get estimated signal based on posterior mean parameters (MCMC)
  
  # Posterior means - MCMC
  MCMC_means <- filter(df, phase_mask == paste(combos$phase[j],combos$mask[j],sep = "_"))
  nu_A = MCMC_means$mean[1:4]
  nu_B = MCMC_means$mean[5:8]
  nu_C = MCMC_means$mean[9:10]
  
  nu <- list(nu_A, nu_B, nu_C)
  
  # Initial value - MCMC
  z0 <- MCMC_means$mean[11:12] 
  
  # Model params
  paramMats = struct_paramMats(m = m,n_u = n_u,idxs = idxs,nu = nu)
  
  # Solve the ODE for z
  out_z <- ode(y = z0,
               times = times,
               func = linear,
               parms= with(paramMats, list(A=A, B=B, C=C)),
               input = input_u,
               atol = 1e-6, 
               rtol = 1e-6
  )
  
  # Get rid of first column of z (same as times)
  out_z <- out_z[,-1]
  
  ##########################################################################
  
  ######### Hemodynamic model #####################################
  # hrf function (using the difference of two gammas - canonical)
  HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
  HRF_mu = function(mu,tp){
    convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
  }
  
  y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  ##########################################################################
  
  # Choose colors for each region
  cols_all <- brewer.pal(8, "Dark2")
  cols <- cols_all[c(1, 3)]
  
  # ylim range
  y_range <- range(c(y_pred, y_obs))
  ylim_exp <- c(y_range[1], y_range[2] * 1.15)
  
  # Set up plot
  plot(times[-1], y_obs[,1], type = "n",  # "n" means no points/lines yet
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
  
  # Save estimated signal to compare with SPM later
  save(y_pred, file = paste0("sub_",sub,"_phase",combos$phase[j],"_mask",combos$mask[j],"_pred.RData"))
  
  # Write observed signal to a MATLAB file for SPM
  writeMat(paste0("sub_",sub,"_phase",combos$phase[j],"_mask",combos$mask[j],".mat"),
           U = dat$u,
           Y = dat$y_obs)
}

# Add Subject title to plot
mtext(paste0("Predicted (Solid) vs. Observed (Dotted): Subject ",sub), line = 0, outer = T)

dev.off()
