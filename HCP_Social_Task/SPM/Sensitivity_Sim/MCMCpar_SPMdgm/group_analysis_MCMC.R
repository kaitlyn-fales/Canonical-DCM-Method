# PEB for combining subjects into group DCM
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(tidyverse))

# Get environment variables from Slurm
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Define all combinations
phases <- c("LR", "RL")
masks <- c("L", "R")

# Expand grid
conditions <- expand.grid(phase = phases, mask = masks, stringsAsFactors = FALSE)

# Pick the corresponding row
phase_condition <- conditions$phase[task_id]
mask_condition  <- conditions$mask[task_id]

cat("Running group analysis for Phase =", phase_condition, "and Mask =", mask_condition, "\n")

# Load diagnostics df to use indices for looping through
load("diagnostics_compilation.RData")

# Filter for one phase and mask condition - do for all four conditions
df <- filter(diagnostics_df, phase == phase_condition & mask == mask_condition)

# Set up MCMC posterior outputs for meta-analysis
###############################################
y_list <- list()
S_list <- list()

for (i in 1:nrow(df)){
  
    # Load your R posterior draws
    load(paste0("Output/sub-", df$subject[i],
                "_phase", df$phase[i], "_mask", df$mask[i], "_draws.RData"))
    
    # Get rid of unnecessary columns
    draws <- suppressWarnings(draws_df[,c(4:13)])
    
    # Transform the diag(A) draws according to our reparameterization
    # nu_A[1] and nu_A[4]
    draws[,c(1,4)] <- -0.5*exp(draws[,c(1,4)])
    
    # Extract posterior means
    post_means <- summarize_draws(draws, mean)[,2]
    y_list[[i]] <- as.numeric(as.matrix(post_means))
    
    # Extract posterior covariance
    S_list[[i]] <- cov(draws)

}

# Load in HCP covariates
covariates <- read.csv("../../../HCP_YA_subjects.csv")

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
mod <- cmdstan_model("../../../../meta_analysis.stan")  

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

# Extract draws
draws <- as_draws_df(fit$draws())

# Save
save(draws, file = paste0("Output/phase",phase_condition,"_mask",mask_condition,"_MCMC.RData"))


