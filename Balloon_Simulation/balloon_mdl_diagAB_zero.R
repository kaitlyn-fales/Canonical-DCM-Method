remove(list = ls())

args <- commandArgs(TRUE)
if (length(args) == 0) {
  stop("No arguments supplied.")
} else {
  for (i in seq_along(args)) {
    eval(parse(text = args[[i]]))
  }
}

suppressPackageStartupMessages(library(cdcm))

# index = simulation replicate
rep_id <- index

# Load simulated data
load(sprintf("Data/balloon_sim_diagAB_zero_data_rep%03d.RData", rep_id))

# Output specs
output_dir <- "Output"
basename <- sprintf("balloon_mdl_diagAB_zero_rep%03d", rep_id)

# Set up hypotheses
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2,
                              2,2,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs)

# Convergence criterion - Multivariate ESS (95% HPD, 5% tolerance)
ess_check <- minESS_criterion(stan_dat, alpha = 0.05, eps = 0.05)

# Compile stan program
canonical_dcm <- compile_cdcm()

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(canonical_dcm, stan_dat)

# Run sampler until convergence
mcmc_seed <- 200000 + rep_id
results <- dcm_sample(
  mod = canonical_dcm,
  data = stan_dat,
  inits_list = inits_list,
  output_dir = output_dir,
  basename = basename,
  ess_check = ess_check,
  seed = mcmc_seed
)

