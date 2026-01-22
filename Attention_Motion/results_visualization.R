# Examining results from SPM attention to motion data across methods

# Packages
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(posterior)

# Summarize results from MCMC
load("Attention_Motion/Output/attention_motion_draws.RData")

results <- summarize_draws(draws_df)
MCMC <- cbind(results[c(5:8,12:14),c(1:2,6:7)],CI_length = abs(results$q95-results$q5)[c(5:8,12:14)])
MCMC <- cbind(MCMC[,1],round(MCMC[,2:5], digits = 4))
write.csv(MCMC, "Attention_Motion/Output/MCMC_summary.csv")

# To use for est vs obs - initial value means
z0_mean <- results[15:17,1:2]
save(z0_mean, file = "Attention_Motion/Output/MCMC_z0_mean.RData")

# To use for est vs obs - self loop param on A
diag_A <- results[9:11,1:2]
save(diag_A, file = "Attention_Motion/Output/MCMC_diag_A.RData")

# Trace plots
bayesplot::mcmc_trace(draws_df[,5:14])

# Load in summaries
SPM <- read.csv("Attention_Motion/SPM_VBA_Comparison/SPM_summary.csv")
VBA <- read.csv("Attention_Motion/SPM_VBA_Comparison/VBA_summary.csv")

# Rename so all columns match
names(SPM)[names(SPM) == 'X'] <- 'par'
SPM$Method <- rep("SPM",nrow(SPM))
VBA <- data.frame(par = SPM$par,VBA[,2:5])
VBA$Method <- rep("VBA",nrow(VBA))
MCMC <- data.frame(par = SPM$par,MCMC[,2:5])
MCMC$Method <- rep("Canonical DCM",nrow(MCMC))

# One df with all results
df <- rbind(SPM,VBA,MCMC)

par.labs <- c("V5 -> V1", 
              "V1 -> V5", 
              "SPC -> V5",
              "V5 -> SPC",
              "Motion on V1 -> V5",
              "Attention on V1 -> V5",
              "Photic -> V1")
names(par.labs) <- SPM$par

# Plot results
df %>% 
  ggplot(mapping = aes(x = par, y = mean, color = Method)) +
  geom_point(shape = 19, size = 2, position=position_dodge(width=0.8)) +
  geom_errorbar(aes(ymax = q95, ymin = q5), 
                position=position_dodge(width=0.8),
                width = 0.4, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = NULL, y = "Posterior Mean Estimate") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.87,0.25),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size = 12)) +
  facet_wrap(~par, scales = "free_x", nrow = 2,
             labeller = labeller(par = par.labs))


