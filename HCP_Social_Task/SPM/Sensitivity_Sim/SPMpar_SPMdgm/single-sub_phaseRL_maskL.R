remove(list=ls())

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No data file provided.")
}

data_file <- args[1]

# Load the data
load(data_file)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mcmcse))
suppressPackageStartupMessages(library(momentLS))
suppressPackageStartupMessages(library(cdcm))

######## Change specifications here ###########
# Subject specs
phase <- "RL"
mask_type <- "L"
sub <- sub("sub-(.*)_phaseRL_maskL\\.RData", "\\1", basename(data_file))

# Output specs
output_dir <- "Output"
basename <- paste0("sub-",sub,"_phase",phase,"_mask",mask_type)

###############################################

paste0("DCM for Subject ",sub,", Phase ", phase," Encoding, and Mask ",mask_type)

########### Get data ready ####################
# Compile stan program from parent directory
canonical_dcm <- compile_cdcm()

# Indices of parameters in hypothesis
A_idxs <- matrix(c(1,1,
                   2,1,
                   1,2,
                   2,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,1,1,
                   2,2,1,
                   2,1,2,
                   2,2,2), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1,
                   2,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs)

# Convergence check specs - 95% intervals with 5% tolerance
ess_check <- minESS_criterion(stan_dat, alpha = 0.05, eps = 0.05)
###############################################

########### Initialize sampler ################
# Initial values for pathfinder from group level DCM paper
pathfinder_inits <- list(sigma = c(1,1),
                         nu_A = c(0,0,0,0),
                         nu_B = c(-0.03,1.22,0.16,-0.19),
                         nu_C = c(0.96,-0.06),
                         z0 = c(0.1,0.1),
                         beta = c(0,0))

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(canonical_dcm, stan_dat, pathfinder_init = pathfinder_inits)

# Run sampler until convergence
result <- dcm_sample(mod = canonical_dcm, 
                     data = stan_dat, 
                     inits_list = inits_list, 
                     output_dir = output_dir,
                     basename = basename,
                     refresh = 100,
                     warmup_iter = 5000,
                     n_iter_chunk = 1000,
                     max_iter = 100000,
                     adapt_delta = 0.9,
                     seed = 1234,
                     ess_check = ess_check)
###############################################




