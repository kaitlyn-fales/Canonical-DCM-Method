# Comparison of estimated signal using posterior means from each approach with observed signal

setwd("/storage/work/krf5429/Canonical-DCM-Method/Balloon_Simulation")

library(deSolve)


######### Diagonal of A only ####################

# Load in posterior means (MCMC)
load("balloon_post_means_diagA.RData")

# Load in observed signals (data)
load("balloon_sim_data_diagA.RData")
u <- rbind(0,dat$u)
times <- c(0,dat$times)
y_obs <- dat$y_obs

# 2 nodes, 2 experimental inputs
m = 2; n_u = 2

# Indices of parameters
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

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
    approxfun(x = times,y = u[,2],rule = 2)(t))
}

##########################################################################
# Get estimated signal based on posterior mean parameters (MCMC)

# Posterior means - MCMC
nu_A = post_means$mean[1:4]
nu_B = post_means$mean[5]
nu_C = post_means$mean[6]

nu <- list(nu_A, nu_B, nu_C)

# Initial value - MCMC
z0 <- post_means$mean[7:8]

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

##########################################################################

######### Hemodynamic model #####################################
# hrf function (using the difference of two gammas - canonical)
HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
HRF_mu = function(mu,tp){
  convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
}

y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
##########################################################################


titles <- c("V1","V2")
xlabs <- c("Time (Seconds)","Time (Seconds)")

par(mfrow = c(1,2))
for (i in 1:ncol(y_pred)){
  plot(times[-1],y_obs[,i],type = 'l',xlab = xlabs[i],
       ylab = "BOLD Response", 
       ylim = c(min(c(y_obs,y_pred)),
                max(c(y_obs,y_pred))), col = "black",
       main = titles[i], lty = "dotted")
  lines(times[-1],y_pred[,i],type = 'l',col = "blue")
}
mtext("MCMC Estimated BOLD Response Using Balloon Hemodynamic Data Generating Model (Diag A)", 
      line = -1.5, outer = T)

# V1
mean((y_obs[,1] - y_pred[,1])^2)

# V2
mean((y_obs[,2] - y_pred[,2])^2)

rm(list = ls())

######### Diagonal of A and B ####################

# Load in posterior means (MCMC)
load("balloon_post_means_diagAB.RData")

# Load in observed signals (data)
load("balloon_sim_data_diagAB.RData")
u <- rbind(0,dat$u)
times <- c(0,dat$times)
y_obs <- dat$y_obs

# 2 nodes, 2 experimental inputs
m = 2; n_u = 2

# Indices of parameters
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2,
                              2,2,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

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
    approxfun(x = times,y = u[,2],rule = 2)(t))
}

##########################################################################
# Get estimated signal based on posterior mean parameters (MCMC)

# Posterior means - MCMC
nu_A = post_means$mean[1:4]
nu_B = post_means$mean[5:6]
nu_C = post_means$mean[7]

nu <- list(nu_A, nu_B, nu_C)

# Initial value - MCMC
z0 <- post_means$mean[8:9]

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

##########################################################################

######### Hemodynamic model #####################################
# hrf function (using the difference of two gammas - canonical)
HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
HRF_mu = function(mu,tp){
  convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
}

y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
##########################################################################


titles <- c("V1","V2")
xlabs <- c("Time (Seconds)","Time (Seconds)")

par(mfrow = c(1,2))
for (i in 1:ncol(y_pred)){
  plot(times[-1],y_obs[,i],type = 'l',xlab = xlabs[i],
       ylab = "BOLD Response", 
       ylim = c(min(c(y_obs,y_pred)),
                max(c(y_obs,y_pred))), col = "black",
       main = titles[i], lty = "dotted")
  lines(times[-1],y_pred[,i],type = 'l',col = "blue")
}
mtext("MCMC Estimated BOLD Response Using Balloon Hemodynamic Data Generating Model (Diag A and B)", 
      line = -1.5, outer = T)

# V1
mean((y_obs[,1] - y_pred[,1])^2)

# V2
mean((y_obs[,2] - y_pred[,2])^2)


