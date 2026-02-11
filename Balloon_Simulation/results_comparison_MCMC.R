library(posterior)
library(grid)
library(gridExtra)


############### Balloon data generating model #################################

# Simple model with diag A est, zero hemodynamics
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7,0,0,0,0)
nu = truth

# Diagonal of A estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("draws", temp)][16:20]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
coverage <- numeric()
length_nu <- numeric()
length <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- summarise_draws(draws_df)
  summary <- summary[4:13,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:6])
  coverage[i] <- mean(summary$coverage)
  
  # Length
  length_nu[i] <- mean(summary$length[1:6])
  length[i] <- mean(summary$length)
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_coverage = mean(coverage),
                                     mean_length_nu = mean(length_nu),
                                     mean_length = mean(length))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q5   = mean(q5),
    q95  = mean(q95),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

# Export resulting post. means for est vs. obs
post_means <- df[,3:4]
save(post_means, file = paste(path,"balloon_post_means_diagA_zero.RData",sep = "/"))

# Plot resulting average estimated nu vs true nu
p1 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A, Zero Hemodynamics") +
  labs(x = "", y = "") + ylim(c(-0.8,1.8)) + 
  annotate("text", x = 2.5, y = -0.8, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -0.8, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -0.8, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 1.7, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.7, y = 1.7, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 1.55, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.7, y = 1.55, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14)) 

############################################################################
rm(list = setdiff(ls(), "p1"))

# Simple model with diag A and diag B est and zero hemodynamics
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,-0.5*-0.5*exp(0.15)*exp(0.05-1),0.7,0,0,0,0)
nu = truth

# Diagonal of A and diag B estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("draws", temp)][6:10]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
coverage <- numeric()
length_nu <- numeric()
length <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- summarise_draws(draws_df)
  summary <- summary[4:14,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:7])
  coverage[i] <- mean(summary$coverage)
  
  # Length
  length_nu[i] <- mean(summary$length[1:7])
  length[i] <- mean(summary$length)
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_coverage = mean(coverage),
                                     mean_length_nu = mean(length_nu),
                                     mean_length = mean(length))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q5   = mean(q5),
    q95  = mean(q95),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

# Export resulting post. means for est vs. obs
post_means <- df[,3:4]
save(post_means, file = paste(path,"balloon_post_means_diagAB_zero.RData",sep = "/"))

# Plot resulting average estimated nu vs true nu
p2 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 6.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A and B, Zero Hemodynamics") +
  labs(x = "", y = "") + ylim(c(-0.8,1.8)) + 
  annotate("text", x = 2.5, y = -0.8, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -0.8, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 7, y = -0.8, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 1.7, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.7, y = 1.7, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 1.55, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.7, y = 1.55, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14)) 

######################################################################
rm(list = setdiff(ls(), c("p1","p2")))

# Simple model with diag A est, nonzero hemodynamics
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7,0,0,0,0)
nu = truth

# Diagonal of A estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("draws", temp)][11:15]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
coverage <- numeric()
length_nu <- numeric()
length <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- summarise_draws(draws_df)
  summary <- summary[4:13,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:6])
  coverage[i] <- mean(summary$coverage)
  
  # Length
  length_nu[i] <- mean(summary$length[1:6])
  length[i] <- mean(summary$length)
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_coverage = mean(coverage),
                                     mean_length_nu = mean(length_nu),
                                     mean_length = mean(length))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q5   = mean(q5),
    q95  = mean(q95),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

# Export resulting post. means for est vs. obs
post_means <- df[,3:4]
save(post_means, file = paste(path,"balloon_post_means_diagA_nonzero.RData",sep = "/"))

# Plot resulting average estimated nu vs true nu
p3 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A, Non-Zero Hemodynamics") +
  labs(x = "", y = "") + ylim(c(-0.8,1.8)) + 
  annotate("text", x = 2.5, y = -0.8, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -0.8, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -0.8, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 1.7, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.7, y = 1.7, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 1.55, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.7, y = 1.55, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14)) 

############################################################################
rm(list = setdiff(ls(), c("p1","p2","p3")))

# Simple model with diag A and diag B est and nonzero hemodynamics
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,-0.5*-0.5*exp(0.15)*exp(0.05-1),0.7,0,0,0,0)
nu = truth

# Diagonal of A and diag B estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("draws", temp)][1:5]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
coverage <- numeric()
length_nu <- numeric()
length <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- summarise_draws(draws_df)
  summary <- summary[4:14,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:7])
  coverage[i] <- mean(summary$coverage)
  
  # Length
  length_nu[i] <- mean(summary$length[1:7])
  length[i] <- mean(summary$length)
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_coverage = mean(coverage),
                                     mean_length_nu = mean(length_nu),
                                     mean_length = mean(length))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q5   = mean(q5),
    q95  = mean(q95),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

# Export resulting post. means for est vs. obs
post_means <- df[,3:4]
save(post_means, file = paste(path,"balloon_post_means_diagAB_nonzero.RData",sep = "/"))

# Plot resulting average estimated nu vs true nu
p4 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 6.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A and B, Non-Zero Hemodynamics") +
  labs(x = "", y = "") + ylim(c(-0.8,1.8)) + 
  annotate("text", x = 2.5, y = -0.8, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -0.8, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 7, y = -0.8, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 1.7, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.7, y = 1.7, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 1.55, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.7, y = 1.55, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14)) 

# Plot both
grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2,
             top = textGrob("MCMC Parameter Estimates from Balloon Hemodynamic Data Generating Model",
                            gp = gpar(fontface = "bold", fontsize = 16)),
             bottom = textGrob("Parameter",
                               gp = gpar(fontface = "bold", fontsize = 14)),
             left = textGrob("Average Posterior Mean Estimate",
                             gp = gpar(fontface = "bold", fontsize = 14), rot = 90))


################################################################################

rm(list = ls())

############### Canonical data generating model #################################

# Simple model with diag A est
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7,0.1,0.1,0,0)
nu = truth

# Diagonal of A estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("draws", temp)][11:15]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
coverage <- numeric()
length_nu <- numeric()
length <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- summarise_draws(draws_df)
  summary <- summary[4:13,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:6])
  coverage[i] <- mean(summary$coverage)
  
  # Length
  length_nu[i] <- mean(summary$length[1:6])
  length[i] <- mean(summary$length)
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_coverage = mean(coverage),
                                     mean_length_nu = mean(length_nu),
                                     mean_length = mean(length))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q5   = mean(q5),
    q95  = mean(q95),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

# Export resulting post. means for est vs. obs
post_means <- df[,3:4]
save(post_means, file = paste(path,"canonical_post_means_diagA.RData",sep = "/"))

# Plot resulting average estimated nu vs true nu
p1 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A") +
  labs(x = "", y = "") + ylim(c(-0.8,1.8)) + 
  annotate("text", x = 2.5, y = -0.8, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -0.8, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -0.8, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 1.7, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.7, y = 1.7, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 1.55, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.7, y = 1.55, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14)) 

############################################################################
rm(list = setdiff(ls(), "p1"))

# Simple model with diag A and diag B est
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.05,0.7,0.1,0.1,0,0)
nu = truth

# Diagonal of A and diag B estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("draws", temp)][16:20]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
coverage <- numeric()
length_nu <- numeric()
length <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- summarise_draws(draws_df)
  summary <- summary[4:14,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:7])
  coverage[i] <- mean(summary$coverage)
  
  # Length
  length_nu[i] <- mean(summary$length[1:7])
  length[i] <- mean(summary$length)
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_coverage = mean(coverage),
                                     mean_length_nu = mean(length_nu),
                                     mean_length = mean(length))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q5   = mean(q5),
    q95  = mean(q95),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

# Export resulting post. means for est vs. obs
post_means <- df[,3:4]
save(post_means, file = paste(path,"canonical_post_means_diagAB.RData",sep = "/"))

# Plot resulting average estimated nu vs true nu
p2 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 6.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A and B") +
  labs(x = "", y = "") + ylim(c(-0.8,1.8)) + 
  annotate("text", x = 2.5, y = -0.8, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -0.8, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 7, y = -0.8, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 1.7, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.7, y = 1.7, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 1.55, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.7, y = 1.55, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14)) 

# Plot both
grid.arrange(p1, p2, ncol=2, nrow=1,
             top = textGrob("MCMC Parameter Estimates from Canonical Hemodynamic Data Generating Model",
                            gp = gpar(fontface = "bold", fontsize = 16)),
             bottom = textGrob("Parameter",
                               gp = gpar(fontface = "bold", fontsize = 14)),
             left = textGrob("Average Posterior Mean Estimate",
                             gp = gpar(fontface = "bold", fontsize = 14), rot = 90))


################################################################################


