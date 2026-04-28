library(posterior)
library(grid)
library(gridExtra)
library(gtable)
library(dplyr)
library(ggplot2)

get_legend <- function(myplot) {
  plot_grob <- ggplotGrob(myplot)
  legend_index <- which(sapply(plot_grob$grobs, function(x) x$name) == "guide-box")
  plot_grob$grobs[[legend_index]]
}


############### Balloon data generating model #################################

# Simple model with diag A est, zero delays
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7,0,0,0,0)
nu = truth

# Diagonal of A estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("balloon_mdl_diagA_zero", temp)]
temp <- temp[grepl("draws", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
length_nu <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- posterior::summarise_draws(draws_df, mean, sd, ~quantile(.x, probs = c(0.025, 0.975)))
  summary <- summary[4:13,]
  summary$coverage <- ifelse(truth >= summary$`2.5%` & truth <= summary$`97.5%`, 1, 0)
  summary$length <- abs(summary$`97.5%`-summary$`2.5%`)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:6])
  
  # Length
  length_nu[i] <- mean(summary$length[1:6])
  
  # Export post. means 
  post_means <- summary[,1:2]
  save(post_means, file = sprintf("%s/balloon_post_means_diagA_zero_rep%03d.RData", path, i))

  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_length_nu = mean(length_nu))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q2.5   = mean(`2.5%`),
    q97.5  = mean(`97.5%`),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

df_plot <- df %>%
  filter(startsWith(as.character(variable), "nu_"))

posterior_df <- df_plot %>%
  mutate(Type = "Posterior mean")

truth_df <- df_plot %>%
  transmute(
    x = x,
    nu = nu,
    Type = "True value"
  )

p1 <- ggplot() +
  geom_point(data = posterior_df,
             aes(x = x, y = mean, shape = Type, color = Type),
             size = 4) +
  geom_point(data = truth_df,
             aes(x = x, y = nu, shape = Type, color = Type),
             size = 6) +
  geom_errorbar(data = posterior_df,
                aes(x = x, ymin = q2.5, ymax = q97.5),
                linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A, No Delays") +
  labs(x = "", y = "", shape = "", color = "") +
  ylim(c(-1.2, 2.2)) +
  annotate("text", x = 2.5, y = -1, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -1, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -1, label = "nu[C]", parse = TRUE, size = 8) +
  scale_shape_manual(values = c("Posterior mean" = 1, "True value" = 8)) +
  scale_color_manual(values = c("Posterior mean" = "black", "True value" = "blue")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12)
  )

############################################################################
rm(list = setdiff(ls(), c("p1","get_legend")))

# Simple model with diag A and diag B est and zero delays
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,-0.5*-0.5*exp(0.15)*exp(0.05-1),0.7,0,0,0,0)
nu = truth

# Diagonal of A and diag B estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("balloon_mdl_diagAB_zero", temp)]
temp <- temp[grepl("draws", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
length_nu <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- posterior::summarise_draws(draws_df, mean, sd, ~quantile(.x, probs = c(0.025, 0.975)))
  summary <- summary[4:14,]
  summary$coverage <- ifelse(truth >= summary$`2.5%` & truth <= summary$`97.5%`, 1, 0)
  summary$length <- abs(summary$`97.5%`-summary$`2.5%`)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:7])
  
  # Length
  length_nu[i] <- mean(summary$length[1:7])
  
  # Export post. means 
  post_means <- summary[,1:2]
  save(post_means, file = sprintf("%s/balloon_post_means_diagAB_zero_rep%03d.RData", path, i))
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_length_nu = mean(length_nu))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q2.5   = mean(`2.5%`),
    q97.5  = mean(`97.5%`),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

df_plot <- df %>%
  filter(startsWith(as.character(variable), "nu_"))

posterior_df <- df_plot %>%
  mutate(Type = "Posterior mean")

truth_df <- df_plot %>%
  transmute(
    x = x,
    nu = nu,
    Type = "True value"
  )

p2 <- ggplot() +
  geom_point(data = posterior_df,
             aes(x = x, y = mean, shape = Type, color = Type),
             size = 4) +
  geom_point(data = truth_df,
             aes(x = x, y = nu, shape = Type, color = Type),
             size = 6) +
  geom_errorbar(data = posterior_df,
                aes(x = x, ymin = q2.5, ymax = q97.5),
                linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 6.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A and B, No Delays") +
  labs(x = "", y = "", shape = "", color = "") +
  ylim(c(-1.2, 2.2)) +
  annotate("text", x = 2.5, y = -1, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -1, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 7, y = -1, label = "nu[C]", parse = TRUE, size = 8) +
  scale_shape_manual(values = c("Posterior mean" = 1, "True value" = 8)) +
  scale_color_manual(values = c("Posterior mean" = "black", "True value" = "blue")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12)
  )

######################################################################
rm(list = setdiff(ls(), c("p1","p2","get_legend")))

# Simple model with diag A est, nonzero delays
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7,0,0,0,0)
nu = truth

# Diagonal of A estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("balloon_mdl_diagA_nonzero", temp)]
temp <- temp[grepl("draws", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
length_nu <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- posterior::summarise_draws(draws_df, mean, sd, ~quantile(.x, probs = c(0.025, 0.975)))
  summary <- summary[4:13,]
  summary$coverage <- ifelse(truth >= summary$`2.5%` & truth <= summary$`97.5%`, 1, 0)
  summary$length <- abs(summary$`97.5%`-summary$`2.5%`)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:6])
  
  # Length
  length_nu[i] <- mean(summary$length[1:6])
  
  # Export post. means 
  post_means <- summary[,1:2]
  save(post_means, file = sprintf("%s/balloon_post_means_diagA_nonzero_rep%03d.RData", path, i))
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_length_nu = mean(length_nu))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q2.5   = mean(`2.5%`),
    q97.5  = mean(`97.5%`),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

df_plot <- df %>%
  filter(startsWith(as.character(variable), "nu_"))

posterior_df <- df_plot %>%
  mutate(Type = "Posterior mean")

truth_df <- df_plot %>%
  transmute(
    x = x,
    nu = nu,
    Type = "True value"
  )

p3 <- ggplot() +
  geom_point(data = posterior_df,
             aes(x = x, y = mean, shape = Type, color = Type),
             size = 4) +
  geom_point(data = truth_df,
             aes(x = x, y = nu, shape = Type, color = Type),
             size = 6) +
  geom_errorbar(data = posterior_df,
                aes(x = x, ymin = q2.5, ymax = q97.5),
                linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A, 1s Delays") +
  labs(x = "", y = "", shape = "", color = "") +
  ylim(c(-1.2, 2.2)) +
  annotate("text", x = 2.5, y = -1, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -1, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -1, label = "nu[C]", parse = TRUE, size = 8) +
  scale_shape_manual(values = c("Posterior mean" = 1, "True value" = 8)) +
  scale_color_manual(values = c("Posterior mean" = "black", "True value" = "blue")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12)
  )

############################################################################
rm(list = setdiff(ls(), c("p1","p2","p3","get_legend")))

# Simple model with diag A and diag B est and delays
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,-0.5*-0.5*exp(0.15)*exp(0.05-1),0.7,0,0,0,0)
nu = truth

# Diagonal of A and diag B estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("balloon_mdl_diagAB_nonzero", temp)]
temp <- temp[grepl("draws", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
length_nu <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- posterior::summarise_draws(draws_df, mean, sd, ~quantile(.x, probs = c(0.025, 0.975)))
  summary <- summary[4:14,]
  summary$coverage <- ifelse(truth >= summary$`2.5%` & truth <= summary$`97.5%`, 1, 0)
  summary$length <- abs(summary$`97.5%`-summary$`2.5%`)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:7])
  
  # Length
  length_nu[i] <- mean(summary$length[1:7])
  
  # Export post. means 
  post_means <- summary[,1:2]
  save(post_means, file = sprintf("%s/balloon_post_means_diagAB_nonzero_rep%03d.RData", path, i))
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_length_nu = mean(length_nu))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q2.5   = mean(`2.5%`),
    q97.5  = mean(`97.5%`),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

df_plot <- df %>%
  filter(startsWith(as.character(variable), "nu_"))

posterior_df <- df_plot %>%
  mutate(Type = "Posterior mean")

truth_df <- df_plot %>%
  transmute(
    x = x,
    nu = nu,
    Type = "True value"
  )

p4 <- ggplot() +
  geom_point(data = posterior_df,
             aes(x = x, y = mean, shape = Type, color = Type),
             size = 4) +
  geom_point(data = truth_df,
             aes(x = x, y = nu, shape = Type, color = Type),
             size = 6) +
  geom_errorbar(data = posterior_df,
                aes(x = x, ymin = q2.5, ymax = q97.5),
                linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 6.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A and B, 1s Delays") +
  labs(x = "", y = "", shape = "", color = "") +
  ylim(c(-1.2, 2.2)) +
  annotate("text", x = 2.5, y = -1, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -1, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 7, y = -1, label = "nu[C]", parse = TRUE, size = 8) +
  scale_shape_manual(values = c("Posterior mean" = 1, "True value" = 8)) +
  scale_color_manual(values = c("Posterior mean" = "black", "True value" = "blue")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12)
  )

# Plot all 
plot_list <- list(p1, p2, p3, p4)

shared_legend <- get_legend(plot_list[[1]])

plot_list_nolegend <- lapply(plot_list, function(p) {
  p + theme(legend.position = "none")
})

plots_grid <- arrangeGrob(
  grobs = plot_list_nolegend,
  ncol = 2,
  nrow = 2,
  top = textGrob("CDCM Parameter Estimates from SPM Data Generating Model",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Parameter",
                    gp = gpar(fontface = "bold", fontsize = 14)),
  left = textGrob("Average Posterior Mean Estimate",
                  gp = gpar(fontface = "bold", fontsize = 14), rot = 90)
)

final_plot <- arrangeGrob(
  plots_grid,
  shared_legend,
  ncol = 1,
  heights = c(10, 1)
)

grid.newpage()
grid.draw(final_plot)


################################################################################

rm(list = setdiff(ls(), "get_legend"))

############### Canonical data generating model #################################

# Simple model with diag A est
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7,0.1,0.1,0,0)
nu = truth

# Diagonal of A estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("canonical_mdl_diagA_rep", temp)]
temp <- temp[grepl("draws", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
length_nu <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- posterior::summarise_draws(draws_df, mean, sd, ~quantile(.x, probs = c(0.025, 0.975)))
  summary <- summary[4:13,]
  summary$coverage <- ifelse(truth >= summary$`2.5%` & truth <= summary$`97.5%`, 1, 0)
  summary$length <- abs(summary$`97.5%`-summary$`2.5%`)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:6])
  
  # Length
  length_nu[i] <- mean(summary$length[1:6])
  
  # Export post. means 
  post_means <- summary[,1:2]
  save(post_means, file = sprintf("%s/canonical_post_means_diagA_rep%03d.RData", path, i))
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_length_nu = mean(length_nu))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q2.5   = mean(`2.5%`),
    q97.5  = mean(`97.5%`),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

df_plot <- df %>%
  filter(startsWith(as.character(variable), "nu_"))

posterior_df <- df_plot %>%
  mutate(Type = "Posterior mean")

truth_df <- df_plot %>%
  transmute(
    x = x,
    nu = nu,
    Type = "True value"
  )

p1 <- ggplot() +
  geom_point(data = posterior_df,
             aes(x = x, y = mean, shape = Type, color = Type),
             size = 4) +
  geom_point(data = truth_df,
             aes(x = x, y = nu, shape = Type, color = Type),
             size = 6) +
  geom_errorbar(data = posterior_df,
                aes(x = x, ymin = q2.5, ymax = q97.5),
                linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A") +
  labs(x = "", y = "", shape = "", color = "") +
  ylim(c(-1.1, 1.5)) +
  annotate("text", x = 2.5, y = -1, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -1, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -1, label = "nu[C]", parse = TRUE, size = 8) +
  scale_shape_manual(values = c("Posterior mean" = 1, "True value" = 8)) +
  scale_color_manual(values = c("Posterior mean" = "black", "True value" = "blue")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12)
  )

############################################################################
rm(list = setdiff(ls(), c("p1","get_legend")))

# Simple model with diag A and diag B est
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.05,0.7,0.1,0.1,0,0)
nu = truth

# Diagonal of A and diag B estimated
# Load in file names
path = paste0(getwd(),"/Balloon_Simulation/Output")
temp = list.files(path = path, 
                  pattern="\\.RData$")

temp <- temp[grepl("canonical_mdl_diagAB_rep", temp)]
temp <- temp[grepl("draws", temp)]

files <- character()
for (i in 1:length(temp)){
  files[i] <- paste(path,temp[i],sep = "/")
}

# Empty vectors and lists to store results
coverage_nu <- numeric()
length_nu <- numeric()
param_summaries <- list()

# For loop for average results across simulations
for (i in 1:length(files)){
  load(files[i])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[3] and nu_A[4]
  draws_df[,6:7] <- suppressWarnings(-0.5*exp(draws_df[,6:7]))
  
  summary <- posterior::summarise_draws(draws_df, mean, sd, ~quantile(.x, probs = c(0.025, 0.975)))
  summary <- summary[4:14,]
  summary$coverage <- ifelse(truth >= summary$`2.5%` & truth <= summary$`97.5%`, 1, 0)
  summary$length <- abs(summary$`97.5%`-summary$`2.5%`)
  
  # Coverage
  coverage_nu[i] <- mean(summary$coverage[1:7])
  
  # Length
  length_nu[i] <- mean(summary$length[1:7])
  
  # Export post. means 
  post_means <- summary[,1:2]
  save(post_means, file = sprintf("%s/canonical_post_means_diagAB_rep%03d.RData", path, i))
  
  param_summaries[[i]] <- summary
}

# Average HPD length and coverage
coverage_length_result <- data.frame(mean_coverage_nu = mean(coverage_nu),
                                     mean_length_nu = mean(length_nu))
round(coverage_length_result,digits = 3)

# Combine all chain summaries and get average result for plotting
all_summ <- bind_rows(param_summaries)
param_order <- summary$variable

avg_summ <- all_summ %>%
  mutate(variable = factor(variable, levels = param_order)) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(mean),
    q2.5   = mean(`2.5%`),
    q97.5  = mean(`97.5%`),
    .groups = "drop"
  ) %>%
  arrange(variable)

x <- seq_len(nrow(avg_summ))
df <- data.frame(x, nu, avg_summ)

df_plot <- df %>%
  filter(startsWith(as.character(variable), "nu_"))

posterior_df <- df_plot %>%
  mutate(Type = "Posterior mean")

truth_df <- df_plot %>%
  transmute(
    x = x,
    nu = nu,
    Type = "True value"
  )

p2 <- ggplot() +
  geom_point(data = posterior_df,
             aes(x = x, y = mean, shape = Type, color = Type),
             size = 4) +
  geom_point(data = truth_df,
             aes(x = x, y = nu, shape = Type, color = Type),
             size = 6) +
  geom_errorbar(data = posterior_df,
                aes(x = x, ymin = q2.5, ymax = q97.5),
                linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 6.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A and B") +
  labs(x = "", y = "", shape = "", color = "") +
  ylim(c(-1.1, 1.5)) +
  annotate("text", x = 2.5, y = -1, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -1, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 7, y = -1, label = "nu[C]", parse = TRUE, size = 8) +
  scale_shape_manual(values = c("Posterior mean" = 1, "True value" = 8)) +
  scale_color_manual(values = c("Posterior mean" = "black", "True value" = "blue")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12)
  )

# Plot both
plot_list <- list(p1, p2)

shared_legend <- get_legend(plot_list[[1]])

plot_list_nolegend <- lapply(plot_list, function(p) {
  p + theme(legend.position = "none")
})

plots_grid <- arrangeGrob(
  grobs = plot_list_nolegend,
  ncol = 2,
  nrow = 1,
  top = textGrob("CDCM Parameter Estimates from CDCM Data Generating Model",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Parameter",
                    gp = gpar(fontface = "bold", fontsize = 14)),
  left = textGrob("Average Posterior Mean Estimate",
                  gp = gpar(fontface = "bold", fontsize = 14), rot = 90)
)

final_plot <- arrangeGrob(
  plots_grid,
  shared_legend,
  ncol = 1,
  heights = c(10, 1)
)

grid.newpage()
grid.draw(final_plot)


################################################################################


