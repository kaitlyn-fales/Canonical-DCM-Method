# Results comparison between MCMC and SPM

library(tidyverse)
library(deSolve)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(R.matlab)
library(gtable)

# Truth (SPM real data post means)
load("HCP_Social_Task/SPM/Results/SPM_results.RData")
truth <- SPM_results

# MCMC
load("HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/MCMC_results.RData")

# SPM
load("HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/SPM_results.RData")


################ Comparing posterior means and HPD intervals ##################

plot_titles <- c("Phase LR Encoding and Mask L",
                 "Phase LR Encoding and Mask R",
                 "Phase RL Encoding and Mask L",
                 "Phase RL Encoding and Mask R")

plot_list <- vector("list", length(MCMC_results))

# For loop to make plots
for (i in seq_along(MCMC_results)) {
  true_vals <- truth[[i]]$mean_alpha
  
  MCMC <- MCMC_results[[i]]
  SPM  <- SPM_results[[i]]
  
  MCMC$Method <- "CDCM"
  SPM$Method  <- "SPM"
  
  combined <- rbind(MCMC, SPM)
  
  true_df <- data.frame(
    parameter = MCMC$parameter,
    true_vals = true_vals,
    Method = "True value"
  )
  
  p <- ggplot(combined, aes(x = factor(parameter), y = mean_alpha,
                            group = Method, colour = Method, shape = Method)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = q2.5_alpha, ymax = q97.5_alpha),
                  width = 0.6, position = position_dodge(width = 0.6)) +
    geom_point(data = true_df,
               aes(x = factor(parameter), y = true_vals,
                   colour = Method, shape = Method),
               inherit.aes = FALSE,
               size = 3) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = c(4.5, 8.5), linetype = "dashed") +
    ggtitle(plot_titles[i]) +
    labs(x = "", y = "", colour = "", shape = "") +
    ylim(c(-2.2, 2.3)) +
    scale_color_manual(
      values = c("CDCM" = "#1B9E77",
                 "SPM" = "#D95F02",
                 "True value" = "black"),
      breaks = c("CDCM", "SPM", "True value")
    ) +
    scale_shape_manual(
      values = c("CDCM" = 16,
                 "SPM" = 17,
                 "True value" = 8),
      breaks = c("CDCM", "SPM", "True value")
    ) +
    annotate("text", x = 2.5, y = -2, label = "hat(alpha)[A]", parse = TRUE, size = 7) +
    annotate("text", x = 6.5, y = -2, label = "hat(alpha)[B]", parse = TRUE, size = 7) +
    annotate("text", x = 9.5, y = -2, label = "hat(alpha)[C]", parse = TRUE, size = 7) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    ) +
    guides(
      colour = guide_legend(override.aes = list(
        linetype = c(0, 0, 0)
      )),
      shape = guide_legend(override.aes = list(
        linetype = c(0, 0, 0)
      ))
    )
  
  plot_list[[i]] <- p
}

# Function to extract legend from a ggplot object
get_legend <- function(myplot) {
  plot_grob <- ggplotGrob(myplot)
  legend_index <- which(sapply(plot_grob$grobs, function(x) x$name) == "guide-box")
  plot_grob$grobs[[legend_index]]
}

# Extract legend from first plot
shared_legend <- get_legend(plot_list[[1]])

# Remove legends from all plots
plot_list_nolegend <- lapply(plot_list, function(p) {
  p + theme(legend.position = "none")
})

# Arrange the 4 plots
plots_grid <- arrangeGrob(
  grobs = plot_list_nolegend,
  ncol = 2,
  nrow = 2,
  top = textGrob("SPM Group-Level Parameters and Data Generating Model",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Parameter",
                    gp = gpar(fontface = "bold", fontsize = 15)),
  left = textGrob("Group-Level Posterior Mean Estimate",
                  gp = gpar(fontface = "bold", fontsize = 15), rot = 90)
)

# Stack plots + shared legend
final_plot <- arrangeGrob(
  plots_grid,
  shared_legend,
  ncol = 1,
  heights = c(10, 1.2)
)

grid.newpage()
grid.draw(final_plot)


###################################################################

################ Compiling information for table ##################

result <- list()

for (i in 1:length(truth)){
  true_vals <- truth[[i]]$mean_alpha
  
  MCMC <- cbind(MCMC_results[[i]],true_vals)
  SPM <- cbind(SPM_results[[i]],true_vals)
  
  MCMC$Method <- "CDCM"
  SPM$Method <- "SPM"
  combined <- rbind(MCMC, SPM)
  
  combined$coverage <- ifelse(combined$true_vals >= combined$q2.5_alpha & combined$true_vals <= combined$q97.5_alpha, 1, 0)
  combined$correct_sign <- ifelse(sign(combined$true_vals) == sign(combined$mean_alpha), 1, 0)
  combined$dist_to_truth <- abs(combined$mean_alpha - combined$true_vals)
  
  result_df <- combined %>% group_by(Method) %>% summarize(par_coverage = mean(coverage),
                                                           correct_sign = mean(correct_sign),
                                                           avg_dist_truth = round(mean(dist_to_truth), digits = 3),
                                                           med_dist_truth = round(median(dist_to_truth), digits = 3))
  result_df$condition <- names(truth)[i]
  
  result[[i]] <- result_df
}

bind_rows(result)

###################################################################
