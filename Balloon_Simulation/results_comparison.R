library(posterior)

# Simple model with diag A est
truth <- c(0.4,0.3,-0.1,0.15,-0.2,0.7)
nu = truth

# Diagonal of A estimated
# Load in file names
path = "/storage/work/krf5429/Canonical-DCM-Method/Balloon_Simulation/Output"
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("draws", temp)][6:10]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

coverage <- numeric()
length <- numeric()
for (i in 1:length(files)){
  load(files[i])
  summary <- summarise_draws(draws_df)
  summary <- summary[4:9,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  coverage[i] <- mean(summary$coverage)
  length[i] <- mean(summary$length)
}

# Plot results from one chain to examine
set.seed(1234)
chain <- sample(1:5, size = 1)

x = c(1:6)
load(files[chain])
results <- summarize_draws(draws_df, c(default_summary_measures(),default_mcse_measures()))
df <- data.frame(x,nu,results[4:9,c(2,6:8)])

# Export resulting post. means for est vs. obs
post_means <- results[-(1:3),1:2]
save(post_means, file = "balloon_post_means_diagA.RData")

ggplot(df, aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("MCMC Posterior Mean Fit from Balloon Hemodynamic Data Generating Model (Diag A)") +
  labs(x = "Parameter", y = "Estimate") + ylim(c(-0.3,1.8)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 16)) 

############################################################################
rm(list = ls())

# Simple model with diag A and diag B est
truth <- c(0.4,0.3,-0.1,0.15,-0.2,0.05,0.7)
nu = truth

# Diagonal of A and diag B estimated
# Load in file names
path = "/storage/work/krf5429/Canonical-DCM-Method/Balloon_Simulation/Output"
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("draws", temp)][1:5]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

coverage <- numeric()
length <- numeric()
for (i in 1:length(files)){
  load(files[i])
  summary <- summarise_draws(draws_df)
  summary <- summary[4:10,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  coverage[i] <- mean(summary$coverage)
  length[i] <- mean(summary$length)
}

# Plot results from one chain to examine
set.seed(12345)
chain <- sample(1:5, size = 1)

x = c(1:7)
load(files[chain])
results <- summarize_draws(draws_df, c(default_summary_measures(),default_mcse_measures()))
df <- data.frame(x,nu,results[4:10,c(2,6:8)])

# Export resulting post. means for est vs. obs
post_means <- results[-(1:3),1:2]
save(post_means, file = "balloon_post_means_diagAB.RData")

ggplot(df, aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("MCMC Posterior Mean Fit from Balloon Hemodynamic Data Generating Model (Diag A and B)") +
  labs(x = "Parameter", y = "Estimate") + ylim(c(-0.3,1.8)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 16)) 













