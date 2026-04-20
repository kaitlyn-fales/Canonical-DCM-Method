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

# Load in simulated data
load("Data/dat_simple_mdl_diag_fix.RData")

# Set up hypotheses
idxs = list(A_idxs = matrix(c(2,1,
                              1,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs, ode_solver_type = 0)

# Compile model
mod <- compile_cdcm()

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(mod, stan_dat)

# Running stan program to sample from posterior
fit = mod$sample(
  data = stan_dat,
  init = list(inits_list), 
  refresh = 20, # output frequency
  iter_warmup = 5000, # warm-up iterations
  iter_sampling = 3000, # sampling iterations
  seed = index, # seed for reproducibility
  chains = 1,
  adapt_delta = 0.9,
  save_warmup = TRUE, 
  metric = "dense_e") 

# Draws in output csv
fit$save_output_files(dir = "/storage/work/krf5429/Canonical-DCM-Method/Computational_Efficiency_Sim/Output",
                      basename = paste0("simple_mdl_ckrk_diag_fix_",index), timestamp = F, random = F)