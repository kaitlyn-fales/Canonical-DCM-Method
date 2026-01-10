library(posterior)
library(cmdstanr)

# Simple model
truth <- c(0.4,0.3,0.1,0.1,-0.2,0.3)
nu = truth

# Analytic
# Load in file names
path = "/storage/work/krf5429/Canonical-DCM-Method/Computational_Efficiency_Sim/Output"
temp = list.files(path = path, 
                  pattern="\\.csv$")

# Simple analytic diag var
temp <- temp[grepl("^simple_mdl_analytic_diag_var", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

coverage_analytic <- numeric()
length_analytic <- numeric()
for (i in 1:length(files)){
  result <- read_cmdstan_csv(files[i])
  summary <- posterior::summarise_draws(result$post_warmup_draws)
  summary <- summary[4:9,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  coverage_analytic[i] <- mean(summary$coverage)
  length_analytic[i] <- mean(summary$length)
}

# Time
result_analytic <- read_cmdstan_csv(files)
round(colMeans(result_analytic$time$chains[,2:4])/60,digits=2)

# Plot results from one chain to compare (matching seed)
set.seed(1234)
chain <- sample(1:5, size = 1)

x = c(1:6)
csv_contents <- read_cmdstan_csv(files[chain])
results <- posterior::summarize_draws(csv_contents$post_warmup_draws, 
                                      c(default_summary_measures(),default_mcse_measures()))
df <- data.frame(x,nu,results[4:9,c(2,6:8)])

ggplot(df, aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("Simple Model: Piecewise Analytic Solution (Estimated Diagonal)") +
  labs(x = "Parameter", y = "Estimate") + ylim(c(-0.5,0.5)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 16)) 



# CKRK
# Load in file names
temp = list.files(path = path, 
                  pattern="\\.csv$")

# Simple ckrk diag fix
temp <- temp[grepl("^simple_mdl_ckrk_diag_var", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

coverage_ckrk <- numeric()
length_ckrk <- numeric()
for (i in 1:length(files)){
  result <- read_cmdstan_csv(files[i])
  summary <- posterior::summarise_draws(result$post_warmup_draws)
  summary <- summary[4:9,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  coverage_ckrk[i] <- mean(summary$coverage)
  length_ckrk[i] <- mean(summary$length)
}

result_ckrk <- read_cmdstan_csv(files)
round(colMeans(result_ckrk$time$chains[,2:4])/60,digits=2)

# Plot results from one seed to compare
csv_contents <- read_cmdstan_csv(files[chain]) # seed 4
results <- posterior::summarize_draws(csv_contents$post_warmup_draws, 
                                      c(default_summary_measures(),default_mcse_measures()))
df <- data.frame(x,nu,results[4:9,c(2,6:8)])

ggplot(df, aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("Simple Model: CKRK Numeric Solver (Estimated Diagonal)") +
  labs(x = "Parameter", y = "Estimate") + ylim(c(-0.5,0.5)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 16)) 

# Complex model
truth <- c(c(rep(c(0.4,0.3,0.5,0.1),2),0.4,0.3),rep(0.1,6),rep(-0.2,3),rep(0.3,3))
nu = truth

# Analytic
# Load in file names
temp = list.files(path = path, 
                  pattern="\\.csv$")

# Complex analytic diag var
temp <- temp[grepl("^complex_mdl_analytic_diag_var", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

coverage_analytic <- numeric()
length_analytic <- numeric()
for (i in 1:length(files)){
  result <- read_cmdstan_csv(files[i])
  summary <- posterior::summarise_draws(result$post_warmup_draws)
  summary <- summary[8:29,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  coverage_analytic[i] <- mean(summary$coverage)
  length_analytic[i] <- mean(summary$length)
}

# Time
result_analytic <- read_cmdstan_csv(files)
round(colMeans(result_analytic$time$chains[,2:4])/60,digits=2)

# Make plot
set.seed(12345)
chain <- sample(1:5,1)


x = c(1:22)
csv_contents <- read_cmdstan_csv(files[chain])
results <- posterior::summarize_draws(csv_contents$post_warmup_draws, 
                                      c(default_summary_measures(),default_mcse_measures()))
df <- data.frame(x,nu,results[8:29,c(2,6:8)])

ggplot(df, aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 3) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 4) +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("Complex Model: Piecewise Analytic Solution (Estimated Diagonal)") +
  labs(x = "Parameter", y = "Estimate") + ylim(c(-3,2)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 16))


# CKRK
# Load in file names
temp = list.files(path = path, 
                  pattern="\\.csv$")

# Complex ckrk diag var
temp <- temp[grepl("^complex_mdl_ckrk_diag_var", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

coverage_ckrk <- numeric()
length_ckrk <- numeric()
for (i in 1:length(files)){
  result <- read_cmdstan_csv(files[i])
  summary <- posterior::summarise_draws(result$post_warmup_draws)
  summary <- summary[8:29,]
  summary$coverage <- ifelse(truth >= summary$q5 & truth <= summary$q95, 1, 0)
  summary$length <- abs(summary$q95-summary$q5)
  coverage_ckrk[i] <- mean(summary$coverage)
  length_ckrk[i] <- mean(summary$length)
}

# Time
result_ckrk <- read_cmdstan_csv(files)
round(colMeans(result_ckrk$time$chains[,2:4])/60,digits=2)

# Make plot
set.seed(12345)
chain <- sample(1:5,1)


x = c(1:22)
csv_contents <- read_cmdstan_csv(files[chain])
results <- posterior::summarize_draws(csv_contents$post_warmup_draws, 
                                      c(default_summary_measures(),default_mcse_measures()))
df <- data.frame(x,nu,results[8:29,c(2,6:8)])

ggplot(df, aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 3) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 4) +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("Complex Model: CKRK Numeric Solution (Estimated Diagonal)") +
  labs(x = "Parameter", y = "Estimate") + ylim(c(-3,2)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        plot.title = element_text(size = 16))



