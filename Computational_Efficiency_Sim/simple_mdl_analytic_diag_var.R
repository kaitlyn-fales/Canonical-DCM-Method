remove(list=ls())

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

suppressPackageStartupMessages(library(tidyverse))

# Load in simulated data
load("dat_simple_mdl_diag_var.RData")

# Set up hypotheses
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

# Source functions from parent directory
source("../canonical_dcm_functions.R")

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs, ode_solver_type = 1)

# Compile stan program from parent directory
canonical_dcm = cmdstanr::cmdstan_model("../canonical_dcm.stan")

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(canonical_dcm, stan_dat)

# Running stan program to sample from posterior
fit = canonical_dcm$sample(
  data = stan_dat,
  init = list(inits_list), 
  refresh = 20, # output frequency
  iter_warmup = 3000, # warm-up iterations
  iter_sampling = 5000, # sampling iterations
  seed = index, # seed for reproducibility
  chains = 1,
  adapt_delta = 0.8,
  save_warmup = TRUE, 
  metric = "dense_e") 

# Draws in output csv
fit$save_output_files(dir = "/storage/work/krf5429/Canonical-DCM-Method/Computational_Efficiency_Sim/Output",
                      basename = paste0("simple_mdl_analytic_diag_var_",index), timestamp = F, random = F)