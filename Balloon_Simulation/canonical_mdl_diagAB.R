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
suppressPackageStartupMessages(library(mcmcse))
suppressPackageStartupMessages(library(momentLS))
suppressPackageStartupMessages(library(R.matlab))

# Load in simulated data from Matlab and convert to R compatible format
matlab_dat <- readMat("Data/canonical_sim_data_diagAB.mat")
dat <- list(times = seq(0,300,by = 2)[-1], u = matlab_dat$U, y_obs = matlab_dat$Y)

# Output specs
output_dir <- "Output"
basename <- paste0("canonical_mdl_diagAB",index)

# Set up hypotheses
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2,
                              2,2,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

# Source functions from parent directory
source("../canonical_dcm_functions.R")

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs)

# Convergence check specs - 90% intervals with 5% tolerance
num_param <- get_num_param(stan_dat)
ess_check <- as.numeric(minESS(num_param, alpha = 0.1, eps = 0.05))

# Compile stan program from parent directory
canonical_dcm = cmdstanr::cmdstan_model("../canonical_dcm.stan")

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
           seed = index,
           chains = 1)                     
