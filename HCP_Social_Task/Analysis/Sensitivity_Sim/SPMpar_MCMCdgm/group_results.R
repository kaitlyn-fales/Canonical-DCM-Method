# Results comparison between MCMC and SPM

library(tidyverse)
library(deSolve)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(R.matlab)

# Truth (MCMC real data post means)
load("HCP_Social_Task/SPM/Results/SPM_results.RData")
truth <- SPM_results

# MCMC
load("HCP_Social_Task/Analysis/Sensitivity_Sim/SPMpar_MCMCdgm/MCMC_results.RData")

# SPM
load("HCP_Social_Task/Analysis/Sensitivity_Sim/SPMpar_MCMCdgm/SPM_results.RData")


################ Comparing posterior means and HPD intervals ##################

plot_titles <- c("Phase LR Encoding and Mask L",
                 "Phase LR Encoding and Mask R",
                 "Phase RL Encoding and Mask L",
                 "Phase RL Encoding and Mask R")

# For loop to make plots
for (i in 1:length(MCMC_results)){
  true_vals <- truth[[i]]$mean_alpha
  
  MCMC <- cbind(MCMC_results[[i]],true_vals)
  SPM <- cbind(SPM_results[[i]],true_vals)
  
  MCMC$Method <- "Canonical DCM"
  SPM$Method <- "SPM"
  combined <- rbind(MCMC, SPM)
  
  p <- ggplot(combined, aes(x = factor(parameter), y = mean_alpha, group = Method, colour = Method)) +
    geom_point(aes(shape = Method), size = 2.5, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = q2.5_alpha, ymax = q97.5_alpha, group = Method),
                  width = 0.6, position = position_dodge(width = 0.6)) +
    geom_point(aes(x = factor(parameter), y = true_vals), size = 3, shape = 8, colour = "black") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = c(4.5, 8.5), linetype = "dashed") +
    ggtitle(plot_titles[i]) +
    labs(x = "", y = "") + ylim(c(-2.2,2.2)) +
    scale_color_brewer(palette = "Dark2") +
    annotate("text", x = 2.5, y = -2, label = "hat(alpha)[A]", parse = TRUE, size = 7) +
    annotate("text", x = 6.5, y = -2, label = "hat(alpha)[B]", parse = TRUE, size = 7) +
    annotate("text", x = 9.5, y = -2, label = "hat(alpha)[C]", parse = TRUE, size = 7) +
    annotate("text", x = 1, y = 1.2, label = "True value", hjust = 0, size = 4) +
    annotate("point", x = 0.7, y = 1.2, shape = 8, size = 3, color = "black") +
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

