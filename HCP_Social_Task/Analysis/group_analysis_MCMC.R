# PEB for combining subjects into group DCM
library(posterior)
library(cmdstanr)

# Load diagnostics df to use indices for looping through
load("HCP_Social_Task/Analysis/diagnostics_compilation.RData")

# Filter for one phase and mask condition - do for all four conditions
df <- filter(diagnostics_df, phase == "RL" & mask == "R")

# Set up MCMC posterior outputs for meta-analysis
###############################################
y_list <- list()
S_list <- list()

for (i in 1:nrow(df)){
  
    # Load your R posterior draws
    load(paste0("HCP_Social_Task/Output/sub-", df$subject[i],
                "_phase", df$phase[i], "_mask", df$mask[i], "_draws.RData"))
    
    # Get rid of unnecessary columns
    draws <- suppressWarnings(draws_df[,c(4:13)])
    
    # Extract posterior means
    post_means <- summarize_draws(draws, mean)[,2]
    y_list[[i]] <- as.numeric(as.matrix(post_means))
    
    # Extract posterior covariance
    S_list[[i]] <- cov(draws)

}

# Load in HCP covariates
covariates <- read.csv("HCP_Social_Task/HCP_YA_subjects.csv")

# Center and scale PMAT
pmat <- scale(covariates$PMAT24_A_CR)

# Ensure factor levels as desired 
age <- factor(covariates$Age)     
gender  <- factor(covariates$Gender)

# Effect (sum-to-zero) coding for age
contrasts(age) <- contr.sum(length(levels(age)))
age_mat <- model.matrix(~ age, data = data.frame(age = age))[,-1]

# Center gender to -0.5/+0.5
gender <- as.numeric(gender == levels(gender)[2])
gender <- ifelse(gender == 1, 0.5, -0.5)

X <- cbind(gender,age_mat,pmat)
colnames(X)[5] <- "pmat"

# Get data in proper format
K <- length(y_list)         # number of studies
p <- 10                     # number of parameters per subject
q <- 5                      # number of covariates

meta_data <- list(
  K = K,
  p = p,
  q = q,
  y = y_list,
  S = S_list,
  X = as.matrix(X)          
)

# Compile stan model
mod <- cmdstan_model("meta_analysis.stan")  

init_fun <- function() {
  list(
    alpha = rep(0, p),
    B     = matrix(0, p, q),
    Ltau  = diag(p),
    tau   = rep(0.1, p),          
    Lcorr = diag(p)              
  )
}

# Run pathfinder to get good starting MCMC values
pf <- mod$pathfinder(data = meta_data, init = init_fun, num_paths = 1)

# Number of chains
num_chains <- 5

# Column names
param_names <- colnames(pf$draws())

# Extract a single Pathfinder draw (e.g., the last one)
draw_i <- tail(pf$draws(), 1)

# alpha: vector of length p
alpha_vals <- as.numeric(draw_i[grep("^alpha\\[", param_names)])

# B: matrix p x q
B_vals <- as.numeric(draw_i[grep("^B\\[", param_names)])
B_mat <- matrix(B_vals, nrow = p, ncol = q, byrow = TRUE)

# tau: vector of length p
tau_vals <- as.numeric(draw_i[grep("^tau\\[", param_names)])

# Lcorr: matrix p x p
Lcorr_vals <- as.numeric(draw_i[grep("^Lcorr\\[", param_names)])
Lcorr_mat <- matrix(Lcorr_vals, nrow = p, ncol = p, byrow = TRUE)

# Create a single chain init list
single_init <- list(
  alpha = alpha_vals,
  B     = B_mat,
  tau   = tau_vals,
  Lcorr = Lcorr_mat
)

# Replicate the same init for all chains
init_list <- rep(list(single_init), num_chains)

# Fit model
fit <- mod$sample(
  data = meta_data,
  chains = 5,
  parallel_chains = 5,
  init = init_list,
  iter_warmup = 1000,
  iter_sampling = 5000,
  adapt_delta = 0.8,
  max_treedepth = 10,
  seed = 1234
)

# Extract posterior draws for alpha[1]...alpha[p]
alpha_draws <- fit$draws(variables = paste0("alpha[", 1:p, "]"))

# Summarize
alpha_summary <- summarize_draws(alpha_draws)

# Extract tau draws
tau_draws <- fit$draws(variables = paste0("tau[", 1:p, "]"))

# Summarize
tau_summary <- summarize_draws(tau_draws)

# General summary table
summary_table <- tibble(
  parameter = colnames(S_list[[1]]), # original parameter names
  mean_alpha = alpha_summary$mean,
  sd_alpha = alpha_summary$sd,    
  q5_alpha = alpha_summary$q5,
  q95_alpha = alpha_summary$q95,
  between_study_sd = tau_summary$mean  
)

summary_table
assign(paste0("phase",df$phase[1],"_mask",df$mask[1]), summary_table)

# Clear environment except for results
rm(list = setdiff(ls(), c("phaseLR_maskL","phaseLR_maskR",
                          "phaseRL_maskL","phaseRL_maskR")))

# Combine results into one list once done and export
MCMC_results <- list(phaseLR_maskL = phaseLR_maskL,
                     phaseLR_maskR = phaseLR_maskR,
                     phaseRL_maskL = phaseRL_maskL,
                     phaseRL_maskR = phaseRL_maskR)
save(MCMC_results, file = "HCP_Social_Task/Analysis/Results/MCMC_results.RData")
