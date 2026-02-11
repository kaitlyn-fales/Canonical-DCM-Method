# Results comparison between MCMC and SPM

library(tidyverse)
library(deSolve)
library(RColorBrewer)
library(grid)
library(gridExtra)

# MCMC
load("HCP_Social_Task/Analysis/Results/MCMC_results.RData")

# SPM
load("HCP_Social_Task/SPM/Results/SPM_results.RData")

################ Comparing posterior means and HPD intervals ##################

plot_titles <- c("Phase LR Encoding and Mask L",
                 "Phase LR Encoding and Mask R",
                 "Phase RL Encoding and Mask L",
                 "Phase RL Encoding and Mask R")

# For loop to make plots
for (i in 1:length(MCMC_results)){
  MCMC <- MCMC_results[[i]]
  SPM <- SPM_results[[i]]
  
  MCMC$Method <- "Canonical DCM"
  SPM$Method <- "SPM"
  combined <- rbind(MCMC, SPM)
  
  p <- ggplot(combined, aes(x = factor(parameter), y = mean_alpha, group = Method, colour = Method)) +
    geom_point(aes(shape = Method), size = 2.5, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = q2.5_alpha, ymax = q97.5_alpha, group = Method),
                  width = 0.6, position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = c(4.5, 8.5), linetype = "dashed") +
    ggtitle(plot_titles[i]) +
    labs(x = "", y = "") + ylim(c(-2.2,2.2)) +
    scale_color_brewer(palette = "Dark2") +
    annotate("text", x = 2.5, y = -2, label = "hat(alpha)[A]", parse = TRUE, size = 7) +
    annotate("text", x = 6.5, y = -2, label = "hat(alpha)[B]", parse = TRUE, size = 7) +
    annotate("text", x = 9.5, y = -2, label = "hat(alpha)[C]", parse = TRUE, size = 7) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16),
      legend.position = c(0, 1),       
      legend.justification = c(0, 1), 
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 14),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )
  
  # Assign a name for gridExtra
  assign(paste0("p",i), p)
  
}

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2,
             top = textGrob("Group-Level Posterior Mean Estimates and 95% HPD Intervals",
                            gp = gpar(fontface = "bold", fontsize = 16)),
             bottom = textGrob("Parameter",
                               gp = gpar(fontface = "bold", fontsize = 14)),
             left = textGrob("Group-Level Posterior Mean Estimate",
                             gp = gpar(fontface = "bold", fontsize = 14), rot = 90))


###################################################################

################ Generate predicted signal (MCMC) #################
rm(list = setdiff(ls(), c("MCMC_results","SPM_results")))

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

par(mfrow = c(2,2))
for (j in 1:length(MCMC_results)){
  # Load data for u matrix
  load(paste0("HCP_Social_Task/Data/sub-100307_",names(MCMC_results)[j],".RData"))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  
  # Posterior means
  MCMC <- MCMC_results[[j]]
  
  # Posterior means - MCMC
  nu_A = MCMC$mean_alpha[1:4]
  nu_B = MCMC$mean_alpha[5:8]
  nu_C = MCMC$mean_alpha[9:10]
  
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
  
  y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  
  # Choose colors for each region
  cols_all <- brewer.pal(8, "Dark2")
  cols <- cols_all[c(1, 3)]
  
  # ylim range
  y_range <- range(c(y_pred))
  ylim_exp <- c(y_range[1], y_range[2] * 1.15)
  
  name <- strsplit(names(MCMC_results)[j], split = "_")
  phase <- substr(name[[1]][1], nchar(name[[1]][1]) - 1, nchar(name[[1]][1]))
  mask <- substr(name[[1]][2], nchar(name[[1]][2]), nchar(name[[1]][2]))
  
  # Set up plot
  plot(times[-1], y_pred[,1], type = "n", 
       xlab = "Time (s)", ylab = "Predicted BOLD",
       ylim = ylim_exp, 
       main = paste0("Phase ",phase," Encoding, Mask ",mask))
  
  # Add observed (dashed) and predicted (solid) lines by region
  for (i in 1:ncol(y_pred)) {
    lines(times[-1], y_pred[,i], col = cols[i], lty = 1, lwd = 1.5)  # solid
  }
  
}

mtext("MCMC Group Level Predicted Signal", side = 3, line = -1.5, outer = T)
###################################################################

################ Generate predicted signal (SPM) ##################
rm(list = setdiff(ls(), c("MCMC_results","SPM_results")))

par(mfrow = c(2,2))
for (j in 1:length(SPM_results)){
  # Load data for u matrix
  load(paste0("HCP_Social_Task/Data/sub-100307_",names(SPM_results)[j],".RData"))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  
  # Choose colors for each region
  cols_all <- brewer.pal(8, "Dark2")
  cols <- cols_all[c(1, 3)]
  
  # ylim range
  y_range <- range(c(BOLD))
  ylim_exp <- c(y_range[1], y_range[2] * 1.15)
  
  name <- strsplit(names(SPM_results)[j], split = "_")
  phase <- substr(name[[1]][1], nchar(name[[1]][1]) - 1, nchar(name[[1]][1]))
  mask <- substr(name[[1]][2], nchar(name[[1]][2]), nchar(name[[1]][2]))
  
  # Set up plot
  plot(times[-1], BOLD[,1], type = "n", 
       xlab = "Time (s)", ylab = "Predicted BOLD",
       ylim = ylim_exp,
       main = paste0("Phase ",phase," Encoding, Mask ",mask))
  
  # Add observed (dashed) and predicted (solid) lines by region
  for (i in 1:ncol(BOLD)) {
    lines(times[-1], BOLD[,i], col = cols[i], lty = 1, lwd = 1.5)  # solid
  }
}

mtext("SPM Group Level Predicted Signal", side = 3, line = -1.5, outer = T)

###################################################################

