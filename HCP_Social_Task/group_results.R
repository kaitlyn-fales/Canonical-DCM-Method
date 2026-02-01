# Results comparison between MCMC and SPM

library(tidyverse)
library(RColorBrewer)
library(grid)
library(gridExtra)

# MCMC
load("HCP_Social_Task/Analysis/Results/MCMC_results.RData")

# SPM
load("HCP_Social_Task/SPM/SPM_results.RData")

plot_titles <- c("Phase LR Encoding and Mask L",
                 "Phase LR Encoding and Mask R",
                 "Phase RL Encoding and Mask L",
                 "Phase RL Encoding and Mask R")

# For loop to make plots
for (i in 1:length(MCMC_results)){
  MCMC <- MCMC_results[[i]]
  SPM <- SPM_results[[i]]
  
  # Only plot off-diagonal params in A and B - not directly comparable for three reasons
  # 1. different parameterizations and 2. hemodynamic model affect the scale of A diagonal params
  # 3. Regulate the dynamics, don't have same interpreted units as all other params which are in Hz
  MCMC <- MCMC[-c(3:4,7:8),]
  SPM <- SPM[-c(3:4,7:8),]
  
  MCMC$Method <- "Canonical DCM"
  SPM$Method <- "SPM"
  combined <- rbind(MCMC, SPM)
  
  p <- ggplot(combined, aes(x = factor(parameter), y = mean_alpha, group = Method, colour = Method)) +
    geom_point(aes(shape = Method), size = 5, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = q5_alpha, ymax = q95_alpha, group = Method),
                  width = 0.5, position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = c(2.5, 4.5), linetype = "dashed") +
    ggtitle(plot_titles[i]) +
    labs(x = "", y = "") + ylim(c(-1.3,1.5)) +
    scale_color_brewer(palette = "Dark2") +
    annotate("text", x = 1.5, y = -1.2, label = "hat(nu)[A]", parse = TRUE, size = 8) +
    annotate("text", x = 3.5, y = -1.2, label = "hat(nu)[B]", parse = TRUE, size = 8) +
    annotate("text", x = 5.5, y = -1.2, label = "hat(nu)[C]", parse = TRUE, size = 8) +
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
             top = textGrob("Group-Level Posterior Mean Estimates and 90% HPD Intervals",
                            gp = gpar(fontface = "bold", fontsize = 16)),
             bottom = textGrob("Parameter",
                               gp = gpar(fontface = "bold", fontsize = 14)),
             left = textGrob("Group-Level Posterior Mean Estimate",
                             gp = gpar(fontface = "bold", fontsize = 14), rot = 90))


