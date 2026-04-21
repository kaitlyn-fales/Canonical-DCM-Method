remove(list=ls())

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mcmcse))
suppressPackageStartupMessages(library(momentLS))
suppressPackageStartupMessages(library(cdcm))

# Load data
load("motion_dat.RData")

# Output specs
output_dir <- "Output"
basename <- paste0("attention_motion")

########### Get data ready ####################
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

# Put data into form for sampler - warning is for end of scan, no issue
stan_dat <- get_stan_dat(motion_dat, idxs)

# Convergence criterion - Multivariate ESS (95% HPD, 5% tolerance)
ess_check <- minESS_criterion(stan_dat, alpha = 0.05, eps = 0.05)

# Compile stan program
canonical_dcm <- compile_cdcm()

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(canonical_dcm, stan_dat)

# Run sampler until convergence
results <- dcm_sample(mod = canonical_dcm, 
                      data = stan_dat, 
                      inits_list = inits_list, 
                      output_dir = output_dir,
                      basename = basename,
                      ess_check = ess_check,
                      seed = 1234)



