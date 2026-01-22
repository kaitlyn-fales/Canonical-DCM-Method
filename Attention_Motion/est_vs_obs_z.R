# Comparison of estimated signal using posterior means from each approach with observed signal

library(deSolve)
library(R.matlab)

# Load in posterior means (MCMC)
MCMC <- read.csv("Attention_Motion/Output/MCMC_summary.csv")
MCMC_means <- MCMC$mean

# Load in posterior means - initial condition (MCMC)
load("Attention_Motion/Output/MCMC_z0_mean.RData")

# Load in posterior means - diag(A) (MCMC)
load("Attention_Motion/Output/MCMC_diag_A.RData")

# Load in posterior means (SPM)
SPM <- read.csv("Attention_Motion/SPM_VBA_Comparison/SPM_summary.csv")
SPM_means <- SPM$mean

# Load in posterior means - diag(A) (SPM)
SPM_diag_A <- read.csv("Attention_Motion/SPM_VBA_Comparison/SPM_diag_A.csv")

# Load in observed signals (data)
load("Attention_Motion/motion_dat.RData")
u <- rbind(0,motion_dat$u)
times <- c(0,motion_dat$times)
y_obs <- motion_dat$y_obs

# 3 nodes, 3 experimental inputs
m = 3; n_u = 3

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

struct_paramMats = function(m, n_u, idxs, nu){
  
  A = matrix(data = 0,nrow = m, ncol = m)
  for(i in 1:nrow(idxs$A_idxs)){
    A[idxs$A_idxs[i,1], idxs$A_idxs[i,2]] = with(nu,nu_A[i])
  }
  diag(A) = -0.5*exp(diag(A))
  
  B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
  for(i in 1:nrow(idxs$B_idxs)){
    B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
  }
  
  C = matrix(data=0, nrow = m, ncol = n_u)
  for(i in 1:nrow(idxs$C_idxs)){
    C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
  }
  
  return(list(A=A,B=B,C=C))
  
}

# ODE for neural activation that needs to be solved
linear <- function(t, z, params, input) {
  # t is the current time point in the integration, 
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]]
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+",lapply(1:nrow(u), function(i) u[i,1]*B[[i]]))
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}

# Putting together to form U matrix of experimental inputs
input_u = function(t){
  c(approxfun(x = times,y = u[,1],rule = 2)(t),
    approxfun(x = times,y = u[,2],rule = 2)(t),
    approxfun(x = times,y = u[,3],rule = 2)(t))
}

##########################################################################
# Get estimated signal based on posterior mean parameters (MCMC)

# Posterior means - MCMC
nu_A = c(MCMC_means[1:4],diag_A$mean)
nu_B = MCMC_means[5:6]
nu_C = MCMC_means[7]

nu <- list(nu_A, nu_B, nu_C)

# Initial value - MCMC
z0 <- z0_mean$mean 

# Model params
paramMats = struct_paramMats(m = m,n_u = n_u,idxs = idxs,nu = nu)

# Solve the ODE for z
out_z <- ode(y = z0,
             times = times,
             func = linear,
             parms= with(paramMats, list(A=A, B=B, C=C)),
             input = input_u,
             atol = 1e-6, 
             rtol = 1e-6
)
diagnostics(out_z)

# Get rid of first column of z (same as times)
out_z <- out_z[,-1]

assign("MCMC_out_z", out_z)

##########################################################################

##########################################################################
# Get estimated signal based on posterior mean parameters (SPM)

# Posterior means - SPM
nu_A = c(SPM_means[1:4],SPM_diag_A$mean)
nu_B = SPM_means[5:6]
nu_C = SPM_means[7]

nu <- list(nu_A, nu_B, nu_C)

# Initial value - SPM
z0 <- rep(0,m) 

# Model params
paramMats = struct_paramMats(m = m,n_u = n_u,idxs = idxs,nu = nu)

# Solve the ODE for z
out_z <- ode(y = z0,
             times = times,
             func = linear,
             parms= with(paramMats, list(A=A, B=B, C=C)),
             input = input_u,
             atol = 1e-6, 
             rtol = 1e-6
)
diagnostics(out_z)

# Get rid of first column of z (same as times)
out_z <- out_z[,-1]

assign("SPM_out_z", out_z)

##########################################################################


############################################################
# Plotting results

titles <- c("V1","V5","SPC")
xlabs <- c("Time (Seconds)","Time (Seconds)","Time (Seconds)")

# Plotting resulting z
par(mfrow = c(3,1), mar = c(4.5,6,4.5,1.5))
for (i in 1:ncol(MCMC_out_z)){
  plot(times[-1],MCMC_out_z[-1,i],type = 'l',xlab = xlabs[i],
       ylab = "Neural Signal", 
       ylim = c(min(c(MCMC_out_z[-1,i],SPM_out_z[-1,i])),
                max(c(MCMC_out_z[-1,i],SPM_out_z[-1,i]))), col = "blue",
       main = titles[i],
       cex.axis = 1.3,    
       cex.lab = 1.4,    
       cex.main = 1.6,    
       font.axis = 2,      
       font.lab = 2,       
       font.main = 2)
  lines(times[-1],SPM_out_z[-1,i],type = 'l',col = "red")
}





