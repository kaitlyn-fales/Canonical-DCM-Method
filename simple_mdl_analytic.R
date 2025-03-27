remove(list=ls())

#args=(commandArgs(TRUE))
#if(length(args)==0){
#  print("No arguments supplied.")
#}else{
#  for(i in 1:length(args)){
#    eval(parse(text=args[[i]]))
#  }
#}

load("simdat_simple_mdl_diag_var.RData")

# Function to get changepoints in u when u(t) is any dimension
get_changepoints <- function(input){
  
  # Get dataframe of indices when changes occur
  changes <- as.data.frame(which(abs(diff(input)) >= 1, arr.ind = T))
  changes$row <- changes$row+1 # add 1 to row index to indicate the row the change starts
  changes <- changes[order(changes$row),] # ascending order
  
  # Make duplicates their own category to reflect when switching from input1 -> input2
  changes <- changes[!duplicated(changes$row),] # remove duplicates
  
  # Add in last row to reflect final block of times
  changes <- rbind(changes, c(nrow(input),0))$row
  
  # Return dataframe
  return(changes)
}

# Obtain dataframe of changes based on input
input <- rbind(0,sim_dat$u)
changes <- get_changepoints(input)

# Set up hypotheses
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

# Compile stan program
dcm_stan7 = cmdstanr::cmdstan_model("dcm_v7.stan")

# Set ODE solver type
ode_solver_type=1 # analytic
tol =10^{-5}

# Getting data in proper list structure for stan program
stan_dat = with(sim_dat,
                list(T=length(times),
                     m = ncol(y_obs),
                     n_u =ncol(u),
                     n_changes = length(changes),
                     change_pts = changes,
                     d_A=nrow(idxs$A_idxs),
                     d_B=nrow(idxs$B_idxs),
                     d_C=nrow(idxs$C_idxs),
                     sigma_nu = 1,
                     sigma_z0 = 0.2,
                     A_idxs = idxs$A_idxs,
                     B_idxs = idxs$B_idxs,
                     C_idxs = idxs$C_idxs,
                     tp= times,
                     u = u,
                     conv = 1, 
                     y_obs = y_obs,
                     rel_tol = tol,
                     abs_tol = tol,
                     max_num_steps = 10^6,
                     ode_solver_type=ode_solver_type))

# Running stan program to sample from posterior
fit = dcm_stan7$sample(
  data=stan_dat,
  init=list(list(nu_A = c(0,0,1,1))), # initialize at 0
  refresh = 1, # output frequency
  iter_warmup = 1000, # warm-up iterations
  iter_sampling = 500, # sampling iterations
  seed = 1234, # seed for reproducibility
  chains = 1,
  adapt_delta = 0.8,
  save_warmup = TRUE) 

fit

# draws in output csv
fit$save_output_files(dir = "/storage/work/krf5429/Comprehensive_Exam/new_DCM")

