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
suppressPackageStartupMessages(library(cdcm))

# Load in simulated data from Matlab and convert to R compatible format
matlab_dat <- readMat("Data/canonical_sim_data_diagA.mat")
dat <- list(times = seq(0,300,by = 2)[-1], u = matlab_dat$U, y_obs = matlab_dat$Y)

# Output specs
output_dir <- "Output"
basename <- paste0("canonical_mdl_diagA",index)

# Set up hypotheses
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs)

# Convergence criterion - Multivariate ESS (95% HPD, 5% tolerance)
minESS_check <- minESS_criterion(stan_dat, alpha = 0.05, eps = 0.05)

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
                      seed = index)                     
