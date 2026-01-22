library(posterior)
library(cmdstanr)
library(grid)
library(gridExtra)

# Simple model
truth <- c(0.4,0.3,0.1,0.1,-0.2,0.3,0.1,0.1,0,0)
nu = truth

# Analytic
##############################################################
# Load in file names
path = paste0(getwd(),"/Computational_Efficiency_Sim/Output")
temp = list.files(path = path, 
                  pattern="\\.csv$")

# Simple analytic diag fix
temp <- temp[grepl("^simple_mdl_analytic_diag_var", temp)]

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
  result <- read_cmdstan_csv(files[i])
  summary <- posterior::summarise_draws(result$post_warmup_draws)
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

# Time
time <- read_cmdstan_csv(files)
round(colMeans(time$time$chains[,2:4])/60,digits=2)

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

p1 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("Piecewise Analytic Solution") +
  labs(x = "", y = "") + ylim(c(-0.8,1)) + 
  annotate("text", x = 2.5, y = -0.6, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -0.6, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -0.6, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 0.9, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.7, y = 0.9, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 0.75, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.7, y = 0.75, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14))

rm(list = setdiff(ls(), c("path","nu","truth","p1")))
####################################################################

# CKRK
################################################################
# Load in file names
temp = list.files(path = path, 
                  pattern="\\.csv$")

# Simple ckrk diag fix
temp <- temp[grepl("^simple_mdl_ckrk_diag_var", temp)]

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
  result <- read_cmdstan_csv(files[i])
  summary <- posterior::summarise_draws(result$post_warmup_draws)
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

# Time
time <- read_cmdstan_csv(files)
round(colMeans(time$time$chains[,2:4])/60,digits=2)

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

p2 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("CKRK Numeric Solver") +
  labs(x = "", y = "") + ylim(c(-0.8,1)) + 
  annotate("text", x = 2.5, y = -0.6, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -0.6, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -0.6, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 0.9, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.7, y = 0.9, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 0.75, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.7, y = 0.75, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14))

# Plot both
grid.arrange(p1, p2, ncol=2, nrow=1,
             top = textGrob("Estimated Diagonal of A Matrix",
                            gp = gpar(fontface = "bold", fontsize = 16)),
             bottom = textGrob("Parameter",
                               gp = gpar(fontface = "bold", fontsize = 14)),
             left = textGrob("Average Posterior Mean Estimate",
                             gp = gpar(fontface = "bold", fontsize = 14), rot = 90))

rm(list = setdiff(ls(), "path"))

#########################################################################


# Complex model
truth <- c(c(rep(c(0.4,0.3,0.5,0.1),2),0.4,0.3),rep(0.1,6),rep(-0.2,3),rep(0.3,3),
           rep(0.1,6),rep(0,6))
nu = truth

# Analytic
#####################################################################
# Load in file names
temp = list.files(path = path, 
                  pattern="\\.csv$")

# Complex analytic diag fix
temp <- temp[grepl("^complex_mdl_analytic_diag_var", temp)]

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
  result <- read_cmdstan_csv(files[i])
  summary <- posterior::summarise_draws(result$post_warmup_draws)
  summary <- summary[8:41,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:22])
  coverage[i] <- mean(summary$coverage)
  
  # Length
  length_nu[i] <- mean(summary$length[1:22])
  length[i] <- mean(summary$length)
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_coverage = mean(coverage),
                                     mean_length_nu = mean(length_nu),
                                     mean_length = mean(length))
round(coverage_length_result,digits = 3)

# Time
time <- read_cmdstan_csv(files)
round(colMeans(time$time$chains[,2:4])/60,digits=2)

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

p1 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 16.5, linetype = "dashed") +
  geom_vline(xintercept = 19.5, linetype = "dashed") +
  ggtitle("Piecewise Analytic Solution") +
  labs(x = "", y = "") + ylim(c(-2,2)) + 
  annotate("text", x = 8.5, y = -1.6, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 18, y = -1.6, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 21, y = -1.6, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 1.9, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.8, y = 1.9, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 1.6, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.8, y = 1.6, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14))


rm(list = setdiff(ls(), c("path","nu","truth","p1")))
#########################################################################

# CKRK
#####################################################################
# Load in file names
temp = list.files(path = path, 
                  pattern="\\.csv$")

# Complex analytic diag fix
temp <- temp[grepl("^complex_mdl_ckrk_diag_var", temp)]

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
  result <- read_cmdstan_csv(files[i])
  summary <- posterior::summarise_draws(result$post_warmup_draws)
  summary <- summary[8:41,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:22])
  coverage[i] <- mean(summary$coverage)
  
  # Length
  length_nu[i] <- mean(summary$length[1:22])
  length[i] <- mean(summary$length)
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_coverage = mean(coverage),
                                     mean_length_nu = mean(length_nu),
                                     mean_length = mean(length))
round(coverage_length_result,digits = 3)

# Time
time <- read_cmdstan_csv(files)
round(colMeans(time$time$chains[,2:4])/60,digits=2)

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

p2 <- df %>% 
  filter(startsWith(as.character(variable), "nu_")) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 16.5, linetype = "dashed") +
  geom_vline(xintercept = 19.5, linetype = "dashed") +
  ggtitle("CKRK Numeric Solver") +
  labs(x = "", y = "") + ylim(c(-2,2)) + 
  annotate("text", x = 8.5, y = -1.6, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 18, y = -1.6, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 21, y = -1.6, label = "nu[C]", parse = TRUE, size = 8) +
  annotate("point", x = 0.5, y = 1.9, shape = 1, size = 4, color = "black") +
  annotate("text", x = 0.8, y = 1.9, label = "Posterior mean", hjust = 0, size = 5) +
  annotate("point", x = 0.5, y = 1.6, shape = 8, size = 6, color = "blue") +
  annotate("text", x = 0.8, y = 1.6, label = "True value", hjust = 0, size = 5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 14))

# Plot both
grid.arrange(p1, p2, ncol=1, nrow=2,
             top = textGrob("Estimated Diagonal of A Matrix",
                            gp = gpar(fontface = "bold", fontsize = 16)),
             bottom = textGrob("Parameter",
                               gp = gpar(fontface = "bold", fontsize = 14)),
             left = textGrob("Average Posterior Mean Estimate",
                             gp = gpar(fontface = "bold", fontsize = 14), rot = 90))


