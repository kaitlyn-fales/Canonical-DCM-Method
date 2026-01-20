remove(list=ls())

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mcmcse))
suppressPackageStartupMessages(library(momentLS))

# Load data
load("motion_dat.RData")

# Output specs
output_dir <- "Output"
basename <- paste0("attention_motion")

# Source functions from parent directory
source("../canonical_dcm_functions.R")

########### Get data ready ####################
# Compile stan program from parent directory
canonical_dcm = cmdstanr::cmdstan_model("../canonical_dcm.stan")

# Indices of parameters
A_idxs <- matrix(c(1,2,
                   2,1,
                   2,3,
                   3,2,
                   1,1,
                   2,2,
                   3,3), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,2,1,
                   3,2,1), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

# Put data into form for sampler
stan_dat <- get_stan_dat(motion_dat, idxs)

# Convergence check specs - 90% intervals with 5% tolerance
num_param <- get_num_param(stan_dat)
ess_check <- as.numeric(minESS(num_param, alpha = 0.1, eps = 0.05))

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(canonical_dcm, stan_dat) 

# Run sampler until convergence
dcm_sample(mod = canonical_dcm, 
           data = stan_dat, 
           inits_list = inits_list, 
           output_dir = output_dir,
           basename = basename,
           metric = "dense_e",
           refresh = 100,
           warmup_iter = 5000,
           n_iter_chunk = 1000,
           max_iter = 100000,
           adapt_delta = 0.9,
           seed = 1234,
           chains = 1)




