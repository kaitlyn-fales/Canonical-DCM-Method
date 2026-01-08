# Comparison of estimated signal using posterior means from each approach with observed signal

setwd("/storage/work/krf5429/attention_motion")

library(deSolve)
library(R.matlab, lib.loc = "/storage/work/krf5429/R_packages")

# Load in posterior means (MCMC)
MCMC <- read.csv("MCMC_summary.csv")
MCMC_means <- MCMC$mean

# Load in posterior means - initial condition (MCMC)
load("MCMC_z0_mean.RData")

# Load in posterior means (SPM)
SPM <- read.csv("SPM_summary.csv")
SPM_means <- SPM$mean

# Load in observed signals (data)
load("motion_dat.RData")
u <- rbind(0,motion_dat$u)
times <- c(0,motion_dat$times)
y_obs <- motion_dat$y_obs

##########################################################################
# Get estimated signal based on posterior mean parameters

# Posterior means - MCMC
nu_A = MCMC_means[1:4]
nu_B = MCMC_means[5:6]
nu_C = MCMC_means[7]

# Initial value - MCMC
z0 <- z0_mean$mean 

# Posterior means - SPM
nu_A = SPM_means[1:4]
nu_B = SPM_means[5:6]
nu_C = SPM_means[7]

# Initial value - SPM
z0 <- rep(0,m)

##########################################################################

# Putting together to form U matrix of experimental inputs
input_u = function(t){
  c(approxfun(x = times,y = u[,1],rule = 2)(t),
    approxfun(x = times,y = u[,2],rule = 2)(t),
    approxfun(x = times,y = u[,3],rule = 2)(t))
}

# 3 nodes, 3 experimental inputs
m = 3; n_u = 3

# Introduce parameters
A = diag(x = rep(-0.5,m))
B = lapply(1:n_u, function(i){matrix(data=0,nrow=m, ncol=m)})
C = matrix(data = 0,nrow = m,ncol = n_u)

A[1,2] = nu_A[1]
A[2,1] = nu_A[2]
A[2,3] = nu_A[3]
A[3,2] = nu_A[4]
B[[2]][2,1] = nu_B[1]
B[[3]][2,1] = nu_B[2]
C[1,1] = nu_C

# ODE for neural activation that needs to be solved
linear <- function(t, z, params, input) {
  # t is the current time point in the integration, 
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]]
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+",lapply(1:nrow(u), function(i) u[i,1]*B[[i]]))
  #B_all = do.call("+", lapply(1:nrow(u), function(i) u[i,1]*B[[i]])) 
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}

# Model params
params <- list(A = A, B = B, C = C)

# Solve the ODE for z
out_z <- ode(y = z0, times = times, func = linear, parms = params, input = input_u)
diagnostics(out_z)

# Get rid of first column of z (same as times)
out_z <- out_z[,-1]

############################################################

# Save individual z separately
MCMC_out_z <- out_z
SPM_out_z <- out_z

titles <- c("V1","V5","SPC")
xlabs <- c("","","Time (Seconds)")

# Plotting resulting z
par(mfrow = c(3,1), mar = c(2,4,4.5,1.5))
for (i in 1:ncol(MCMC_out_z)){
  plot(times[-1],MCMC_out_z[-1,i],type = 'l',xlab = xlabs[i],
       ylab = "Neural Signal", 
       ylim = c(min(c(MCMC_out_z[-1,i],SPM_out_z[-1,i])),
                max(c(MCMC_out_z[-1,i],SPM_out_z[-1,i]))), col = "blue",
       main = titles[i])
  lines(times[-1],SPM_out_z[-1,i],type = 'l',col = "red")
}
mtext("Estimated Neural Signal (z) Comparison Between MCMC and SPM", 
      line = -1.5, outer = T)





