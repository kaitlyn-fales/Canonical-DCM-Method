library(grid)
library(gridExtra)
library(R.matlab)

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
  alpha = 0.1
  
  intervals <- matrix(NA, nrow=length(mu), ncol=2)
  for (i in 1:nrow(intervals)){
    q_lower <- qnorm(alpha/2, mean = mu[i], sd = sqrt(sigma2[i]))
    q_upper <- qnorm(1-alpha/2, mean = mu[i], sd = sqrt(sigma2[i]))
    intervals[i,] <- c(q_lower, q_upper)
  }
  
  df <- data.frame(mu,intervals)
  colnames(df) <- c("mean","q5","q95")
  
  df$coverage <- ifelse(truth >= df$q5 & truth <= df$q95, 1, 0)
  df$CI_length <- abs(df$q95-df$q5)
  
  return(df)
}
#####################################

# Simple model with diag A est
truth <- c(0.4,0.3,-0.1,0.15,-0.2,0.7)
nu = truth

# Diagonal of A estimated
# Load in file
path = paste0(getwd(),"/Balloon_Simulation/Output")
result <- readMat(paste0(path,"/canonical_diagA_out.mat"))

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

# Make df for plotting
x <- seq_len(nrow(summary))
df <- data.frame(x, nu, summary)

# Plot resulting average estimated nu vs true nu
p1 <- df %>% 
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A") +
  labs(x = "", y = "") + ylim(c(-0.8,1.8)) + 
  annotate("text", x = 2.5, y = -0.7, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5, y = -0.7, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 6, y = -0.7, label = "nu[C]", parse = TRUE, size = 8) +
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
rm(list = setdiff(ls(), c("p1","get_summary","unstruct_paramMats")))

# Simple model with diag A and diag B est
truth <- c(0.4,0.3,-0.1,0.15,-0.2,0.05,0.7)
nu = truth

# Diagonal of A estimated
# Load in file
path = paste0(getwd(),"/Balloon_Simulation/Output")
result <- readMat(paste0(path,"/canonical_diagAB_out.mat"))

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

# Make df for plotting
x <- seq_len(nrow(summary))
df <- data.frame(x, nu, summary)

# Plot resulting average estimated nu vs true nu
p2 <- df %>% 
  ggplot(aes(x = x, y = mean)) +
  geom_point(shape = 1, size = 4) +
  geom_point(aes(x = x, y = nu), shape = 8, color = "blue", size = 6) +
  geom_errorbar(aes(ymax = q95, ymin = q5), linewidth = 0.8, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 4.5, linetype = "dashed") +
  geom_vline(xintercept = 6.5, linetype = "dashed") +
  ggtitle("Diagonal Parameters in A and B") +
  labs(x = "", y = "") + ylim(c(-0.8,1.8)) + 
  annotate("text", x = 2.5, y = -0.7, label = "nu[A]", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -0.7, label = "nu[B]", parse = TRUE, size = 8) +
  annotate("text", x = 7, y = -0.7, label = "nu[C]", parse = TRUE, size = 8) +
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
             top = textGrob("SPM Parameter Estimates from Canonical Hemodynamic Data Generating Model",
                            gp = gpar(fontface = "bold", fontsize = 16)),
             bottom = textGrob("Parameter",
                               gp = gpar(fontface = "bold", fontsize = 14)),
             left = textGrob("Posterior Mean Estimate",
                             gp = gpar(fontface = "bold", fontsize = 14), rot = 90))








