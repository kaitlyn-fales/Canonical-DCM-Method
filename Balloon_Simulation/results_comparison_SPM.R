library(grid)
library(gridExtra)
library(R.matlab)
library(gtable)
library(dplyr)
library(ggplot2)

############ Functions ##############
# Unstructure params function
unstruct_paramMats = function(params, idxs){
  nu_A = params$A[idxs$A_idxs]
  nu_C = params$C[idxs$C_idxs]
  nu_B = c()
  for(i in 1:nrow(idxs$B_idxs)){
    nu_B[i] = params$B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]]
  }
  c(nu_A = nu_A, nu_B = nu_B, nu_C = nu_C)
}

# Function for summary
get_summary <- function(mu,sigma2){
  
  # Set alpha level
  alpha = 0.05
  
  intervals <- matrix(NA, nrow=length(mu), ncol=2)
  for (i in 1:nrow(intervals)){
    q_lower <- qnorm(alpha/2, mean = mu[i], sd = sqrt(sigma2[i]))
    q_upper <- qnorm(1-alpha/2, mean = mu[i], sd = sqrt(sigma2[i]))
    intervals[i,] <- c(q_lower, q_upper)
  }
  
  df <- data.frame(mu,intervals)
  colnames(df) <- c("mean","q2.5","q97.5")
  
  df$coverage <- ifelse(truth >= df$q2.5 & truth <= df$q97.5, 1, 0)
  df$length <- abs(df$q97.5-df$q2.5)
  
  return(df)
}

get_legend <- function(myplot) {
  plot_grob <- ggplotGrob(myplot)
  legend_index <- which(sapply(plot_grob$grobs, function(x) x$name) == "guide-box")
  plot_grob$grobs[[legend_index]]
}

diagA_reparam <- function(result){
  # Inputs - vectorize
  mu    <- c(result$Ep.A, result$Ep.B, result$Ep.C)
  Sigma <- result$Cp[-c(17:20),-c(17:20)] 
  
  # Indices of the params to be transformed
  idx_A <- c(1,4) 
  
  # simulation settings
  S <- 20000
  set.seed(1234)
  
  # Simulate draws
  draws <- MASS::mvrnorm(n = S, mu = mu, Sigma = Sigma)
  
  # Transform the specified components
  draws[,idx_A] <- -0.5*exp(draws[,idx_A])
  
  # Compute empirical mean and covariance on transformed scale
  mu_emp  <- colMeans(draws)     
  Sigma_emp <- cov(draws)   
  
  # Add in transformed diagonal back to original estimate
  diag(result$Ep.A) <- mu_emp[idx_A]
  diag(result$Cp)[idx_A] <- diag(Sigma_emp)[idx_A]
  
  # Indices of parameters
  idxs = list(A_idxs = matrix(c(2,1,
                                1,2,
                                1,1,
                                2,2), byrow = T, ncol=2),
              B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
              C_idxs = matrix(c(1,1), byrow=T, ncol=2))
  
  # Posterior means
  means <- list(A = result$Ep.A, B = result$Ep.B, C = result$Ep.C)
  
  # Grab posterior means
  means$B <- list(means$B[,,1],means$B[,,2])
  
  # All posterior variances coming from A, B, C (hemodynamic params are last 4)
  Cp <- diag(result$Cp)
  Cp <- Cp[1:(length(Cp)-4)]
  
  # Get rid of nonzero entries - correspond to placeholders not est in A, B, C
  Cp <- Cp[Cp != 0]
  
  # Update ordering of post var to match that of nu - current order is column-major order of vec(A,B,C) matrices
  # Column-major order is (1,1), (2,1), (1,2), (2,2), B(2,1,2), C(1,1)
  vars <- Cp[c(2,3,1,4,5,6)]
  
  # Unstructure means
  mu <- unstruct_paramMats(means,idxs)
  
  # Get results
  summary <- get_summary(mu,vars)
  
  return(summary)
}

diagAB_reparam <- function(result){
  # Inputs - vectorize
  mu    <- c(result$Ep.A, result$Ep.B, result$Ep.C)
  Sigma <- result$Cp[-c(17:20),-c(17:20)] 
  
  # Indices of the params to be transformed
  idx_A <- c(1,4) 
  idx_B <- 12
  
  # simulation settings
  S <- 20000
  set.seed(1234)
  
  # Simulate draws
  draws <- MASS::mvrnorm(n = S, mu = mu, Sigma = Sigma)
  
  # Transform the specified components
  draws[,idx_A] <- -0.5*exp(draws[,idx_A])
  draws[,idx_B] <- -0.5*draws[,idx_A[2]]*exp(draws[,idx_B]-1)
  
  # Compute empirical mean and covariance on transformed scale
  mu_emp  <- colMeans(draws)     
  Sigma_emp <- cov(draws)   
  
  # Add in transformed diagonal back to original estimate
  diag(result$Ep.A) <- mu_emp[idx_A]
  diag(result$Cp)[idx_A] <- diag(Sigma_emp)[idx_A]
  
  result$Ep.B[2,2,2] <- mu_emp[idx_B]
  diag(result$Cp)[idx_B] <- diag(Sigma_emp)[idx_B]
  
  # Indices of parameters
  idxs = list(A_idxs = matrix(c(2,1,
                                1,2,
                                1,1,
                                2,2), byrow = T, ncol=2),
              B_idxs = matrix(c(2,1,2,
                                2,2,2), byrow=T, ncol = 3),
              C_idxs = matrix(c(1,1), byrow=T, ncol=2))
  
  # Posterior means
  means <- list(A = result$Ep.A, B = result$Ep.B, C = result$Ep.C)
  
  # Grab posterior means
  means$B <- list(means$B[,,1],means$B[,,2])
  
  # All posterior variances coming from A, B, C (hemodynamic params are last 4)
  Cp <- diag(result$Cp)
  Cp <- Cp[1:(length(Cp)-4)]
  
  # Get rid of nonzero entries - correspond to placeholders not est in A, B, C
  Cp <- Cp[Cp != 0]
  
  # Update ordering of post var to match that of nu - current order is column-major order of vec(A,B,C) matrices
  # Column-major order is (1,1), (2,1), (1,2), (2,2), B(2,1,2), B(2,2,2), C(1,1)
  vars <- Cp[c(2,3,1,4,5,6,7)]
  
  # Unstructure means
  mu <- unstruct_paramMats(means,idxs)
  
  # Get results
  summary <- get_summary(mu,vars)
  
  return(summary)
}
#####################################

#################### Canonical data generating model ###########################
# Simple model with diag A est
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7)
nu = truth

# Diagonal of A estimated
# Load in file
path = paste0(getwd(),"/Balloon_Simulation/Output")

coverage_nu <- numeric()
length_nu <- numeric()

for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/canonical_diagA_rep%03d_out.mat", k))

  summary <- diagA_reparam(result)
  coverage_nu[k] <- mean(summary$coverage)
  length_nu[k] <- mean(summary$length)
}

coverage_length_result <- data.frame(coverage = mean(coverage_nu),
                                     length = mean(length_nu))
round(coverage_length_result,digits = 3)

# Make df for plotting
x <- seq_len(nrow(summary))
df_plot <- data.frame(x, nu, summary)

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
rm(list = setdiff(ls(), c("p1","get_summary","unstruct_paramMats","get_legend","diagA_reparam","diagAB_reparam")))

# Simple model with diag A and diag B est
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.05,0.7)
nu = truth

# Diagonal of AB estimated
# Load in file
path = paste0(getwd(),"/Balloon_Simulation/Output")

coverage_nu <- numeric()
length_nu <- numeric()

for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/canonical_diagAB_rep%03d_out.mat", k))
  
  summary <- diagAB_reparam(result)
  coverage_nu[k] <- mean(summary$coverage)
  length_nu[k] <- mean(summary$length)
}

coverage_length_result <- data.frame(coverage = mean(coverage_nu),
                                     length = mean(length_nu))
round(coverage_length_result,digits = 3)

# Make df for plotting
x <- seq_len(nrow(summary))
df_plot <- data.frame(x, nu, summary)

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
  top = textGrob("SPM Parameter Estimates from CDCM Data Generating Model",
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

rm(list = setdiff(ls(), c("get_summary","unstruct_paramMats","get_legend","diagA_reparam","diagAB_reparam")))

#################### Balloon data generating model ###########################
# Simple model with diag A est
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7)
nu = truth

# Diagonal of A estimated - zero delays
# Load in file
path = paste0(getwd(),"/Balloon_Simulation/Output")

coverage_nu <- numeric()
length_nu <- numeric()

for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/balloon_diagA_zero_rep%03d_out.mat", k))
  
  summary <- diagA_reparam(result)
  coverage_nu[k] <- mean(summary$coverage)
  length_nu[k] <- mean(summary$length)
}

coverage_length_result <- data.frame(coverage = mean(coverage_nu),
                                     length = mean(length_nu))
round(coverage_length_result,digits = 3)

x <- seq_len(nrow(summary))
df_plot <- data.frame(x, nu, summary)

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
rm(list = setdiff(ls(), c("p1","get_summary","unstruct_paramMats","get_legend","diagA_reparam","diagAB_reparam")))

# Simple model with diag A and diag B est, zero delays
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,-0.5*-0.5*exp(0.15)*exp(0.05-1),0.7)
nu = truth

# Diagonal of A and B estimated
# Load in file
path = paste0(getwd(),"/Balloon_Simulation/Output")

coverage_nu <- numeric()
length_nu <- numeric()

for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/balloon_diagAB_zero_rep%03d_out.mat", k))
  
  summary <- diagAB_reparam(result)
  coverage_nu[k] <- mean(summary$coverage)
  length_nu[k] <- mean(summary$length)
}

coverage_length_result <- data.frame(coverage = mean(coverage_nu),
                                     length = mean(length_nu))
round(coverage_length_result,digits = 3)

x <- seq_len(nrow(summary))
df_plot <- data.frame(x, nu, summary)

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
             
################################################################################

rm(list = setdiff(ls(), c("p1","p2","get_summary","unstruct_paramMats","get_legend","diagA_reparam","diagAB_reparam")))

# Simple model with diag A est
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,0.7)
nu = truth

# Diagonal of A estimated - 1s delays
# Load in file
path = paste0(getwd(),"/Balloon_Simulation/Output")

coverage_nu <- numeric()
length_nu <- numeric()

for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/balloon_diagA_nonzero_rep%03d_out.mat", k))
  
  summary <- diagA_reparam(result)
  coverage_nu[k] <- mean(summary$coverage)
  length_nu[k] <- mean(summary$length)
}

coverage_length_result <- data.frame(coverage = mean(coverage_nu),
                                     length = mean(length_nu))
round(coverage_length_result,digits = 3)

x <- seq_len(nrow(summary))
df_plot <- data.frame(x, nu, summary)

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

################################################################################

rm(list = setdiff(ls(), c("p1","p2","p3","get_summary","unstruct_paramMats","get_legend","diagA_reparam","diagAB_reparam")))

# Simple model with diag A and diag B est, nonzero delays
truth <- c(0.4,0.3,-0.5*exp(-0.1),-0.5*exp(0.15),-0.2,-0.5*-0.5*exp(0.15)*exp(0.05-1),0.7)
nu = truth

# Diagonal of A estimated
# Load in file
path = paste0(getwd(),"/Balloon_Simulation/Output")

coverage_nu <- numeric()
length_nu <- numeric()

for (k in 1:50){
  result <- readMat(sprintf("Balloon_Simulation/Output/balloon_diagAB_nonzero_rep%03d_out.mat", k))
  
  summary <- diagAB_reparam(result)
  coverage_nu[k] <- mean(summary$coverage)
  length_nu[k] <- mean(summary$length)
}

coverage_length_result <- data.frame(coverage = mean(coverage_nu),
                                     length = mean(length_nu))
round(coverage_length_result,digits = 3)

x <- seq_len(nrow(summary))
df_plot <- data.frame(x, nu, summary)

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
  top = textGrob("SPM Parameter Estimates from SPM Data Generating Model",
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
