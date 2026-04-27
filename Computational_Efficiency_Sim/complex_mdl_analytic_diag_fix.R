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
suppressPackageStartupMessages(library(cdcm))

# index = simulation replicate
rep_id <- index

# Load in simulated data
load(sprintf("Data/dat_complex_mdl_diag_fix_rep%03d.RData", rep_id))

# Set up hypotheses
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              3,2,
                              2,3,
                              4,3,
                              3,4,
                              5,4,
                              4,5,
                              6,5,
                              5,6), 
                            byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2,
                              2,3,4,
                              2,5,6), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1,
                              3,1,
                              5,1), byrow=T, ncol=2))

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs, ode_solver_type = 1)

# Compile model
mod <- compile_cdcm()

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(mod, stan_dat)

# Running stan program to sample from posterior
mcmc_seed <- 200000 + rep_id
fit = mod$sample(
  data = stan_dat,
  init = list(inits_list), 
  refresh = 20, # output frequency
  iter_warmup = 5000, # warm-up iterations
  iter_sampling = 3000, # sampling iterations
  seed = mcmc_seed, # seed for reproducibility
  chains = 1,
  adapt_delta = 0.9,
  save_warmup = TRUE, 
  metric = "dense_e")

# Draws in output csv
fit$save_output_files(dir = "/storage/work/krf5429/Canonical-DCM-Method/Computational_Efficiency_Sim/Output",
                      basename = sprintf("complex_mdl_analytic_diag_fix_rep%03d", rep_id), 
                      timestamp = F, 
                      random = F)


