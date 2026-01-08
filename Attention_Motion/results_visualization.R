# Examining results from SPM attention to motion data across methods

setwd("/storage/work/krf5429/attention_motion")

# Packages
library(gridExtra, lib.loc = "/storage/work/krf5429/R_packages")
library(grid, lib.loc = "/storage/work/krf5429/R_packages")
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# Load in summaries
SPM <- read.csv("SPM_summary.csv")
VBA <- read.csv("VBA_summary.csv")
MCMC <- read.csv("MCMC_summary.csv")

# Rename so all columns match
names(SPM)[names(SPM) == 'X'] <- 'par'
SPM$Method <- rep("SPM",nrow(SPM))
VBA <- data.frame(par = SPM$par,VBA[,2:5])
VBA$Method <- rep("VBA",nrow(VBA))
MCMC <- data.frame(par = SPM$par,round(MCMC[,3:6], digits = 4))
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


