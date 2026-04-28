# Generate predicted curve for each seed using the posterior means

library(deSolve)
library(R.matlab)

############### Balloon data generating model ###########################
######### Diagonal of A only, no delays ####################

pred_signal <- list()

for (i in 1:50){
  # Load in posterior means (MCMC)
  load(sprintf("Balloon_Simulation/Output/balloon_post_means_diagA_zero_rep%03d.RData", i))
  
  # Load in observed signals (data)
  load(sprintf("Balloon_Simulation/Data/balloon_sim_diagA_zero_data_rep%03d.RData", i))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  
  rm(dat)

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
    
    B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
    for(i in 1:nrow(idxs$B_idxs)){
      B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
    }
    
    C = matrix(data=0, nrow = m, ncol = n_u)
    for(i in 1:nrow(idxs$C_idxs)){
      C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
    }
    
    return(list(A=A,B=B,C=C))
    
  } # imported params already have diag(A) transform
  
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
  
  # Get estimated signal based on posterior mean parameters (MCMC)
  
  # Posterior means - MCMC
  nu_A = post_means$mean[1:4]
  nu_B = post_means$mean[5]
  nu_C = post_means$mean[6]
  
  nu <- list(nu_A, nu_B, nu_C)
  
  # Initial value - MCMC
  z0 <- post_means$mean[7:8]
  
  # Constant shift beta - MCMC
  beta <- post_means$mean[9:10]
  
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
  
  # Get rid of first column of z (same as times)
  out_z <- out_z[,-1]
  
  ## Hemodynamic model 
  # hrf function (using the difference of two gammas - canonical)
  HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
  HRF_mu = function(mu,tp){
    convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
  }
  
  y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  
  # Apply constant shift as estimated by beta
  y_pred = sweep(y_pred, 2, beta, "+")
  
  # Add to list
  pred_signal[[i]] <- y_pred
}

save(pred_signal, file = "Balloon_Simulation/Output/balloon_diagA_zero_pred.RData")

rm(list = ls())
##########################################################################

######### Diagonal of A and B, no delays ####################

pred_signal <- list()

for (i in 1:50){
  # Load in posterior means (MCMC)
  load(sprintf("Balloon_Simulation/Output/balloon_post_means_diagAB_zero_rep%03d.RData", i))
  
  # Load in observed signals (data)
  load(sprintf("Balloon_Simulation/Data/balloon_sim_diagAB_zero_data_rep%03d.RData", i))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  
  rm(dat)
  
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
    
    B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
    for(i in 1:nrow(idxs$B_idxs)){
      B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
    }
    
    C = matrix(data=0, nrow = m, ncol = n_u)
    for(i in 1:nrow(idxs$C_idxs)){
      C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
    }
    
    return(list(A=A,B=B,C=C))
    
  } # imported params already have diag(A) transform
  
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
  
  # Get estimated signal based on posterior mean parameters (MCMC)
  
  # Posterior means - MCMC
  nu_A = post_means$mean[1:4]
  nu_B = post_means$mean[5:6]
  nu_C = post_means$mean[7]
  
  nu <- list(nu_A, nu_B, nu_C)
  
  # Initial value - MCMC
  z0 <- post_means$mean[8:9]
  
  # Constant shift beta - MCMC
  beta <- post_means$mean[10:11]
  
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
  
  # Get rid of first column of z (same as times)
  out_z <- out_z[,-1]
  
  ## Hemodynamic model
  # hrf function (using the difference of two gammas - canonical)
  HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
  HRF_mu = function(mu,tp){
    convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
  }
  
  y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  
  # Apply constant shift as estimated by beta
  y_pred = sweep(y_pred, 2, beta, "+")
  
  pred_signal[[i]] <- y_pred
}

save(pred_signal, file = "Balloon_Simulation/Output/balloon_diagAB_zero_pred.RData")

rm(list = ls())
###############################################################################

######### Diagonal of A only, 1s delays ####################

pred_signal <- list()

for (i in 1:50){
  # Load in posterior means (MCMC)
  load(sprintf("Balloon_Simulation/Output/balloon_post_means_diagA_nonzero_rep%03d.RData", i))
  
  # Load in observed signals (data)
  load(sprintf("Balloon_Simulation/Data/balloon_sim_diagA_nonzero_data_rep%03d.RData", i))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  
  rm(dat)
  
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
    
    B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
    for(i in 1:nrow(idxs$B_idxs)){
      B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
    }
    
    C = matrix(data=0, nrow = m, ncol = n_u)
    for(i in 1:nrow(idxs$C_idxs)){
      C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
    }
    
    return(list(A=A,B=B,C=C))
    
  } # imported params already have diag(A) transform
  
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
  
  # Get estimated signal based on posterior mean parameters (MCMC)
  
  # Posterior means - MCMC
  nu_A = post_means$mean[1:4]
  nu_B = post_means$mean[5]
  nu_C = post_means$mean[6]
  
  nu <- list(nu_A, nu_B, nu_C)
  
  # Initial value - MCMC
  z0 <- post_means$mean[7:8]
  
  # Constant shift beta - MCMC
  beta <- post_means$mean[9:10]
  
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
  
  # Get rid of first column of z (same as times)
  out_z <- out_z[,-1]
  
  ## Hemodynamic model 
  # hrf function (using the difference of two gammas - canonical)
  HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
  HRF_mu = function(mu,tp){
    convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
  }
  
  y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  
  # Apply constant shift as estimated by beta
  y_pred = sweep(y_pred, 2, beta, "+")
  
  pred_signal[[i]] <- y_pred
}

save(pred_signal, file = "Balloon_Simulation/Output/balloon_diagA_nonzero_pred.RData")

rm(list = ls())
##########################################################################

######### Diagonal of A and B, 1s delays ####################

pred_signal <- list()

for (i in 1:50){
  # Load in posterior means (MCMC)
  load(sprintf("Balloon_Simulation/Output/balloon_post_means_diagAB_nonzero_rep%03d.RData", i))
  
  # Load in observed signals (data)
  load(sprintf("Balloon_Simulation/Data/balloon_sim_diagAB_nonzero_data_rep%03d.RData", i))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  
  rm(dat)
  
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
    
    B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
    for(i in 1:nrow(idxs$B_idxs)){
      B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
    }
    
    C = matrix(data=0, nrow = m, ncol = n_u)
    for(i in 1:nrow(idxs$C_idxs)){
      C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
    }
    
    return(list(A=A,B=B,C=C))
    
  } # imported params already have diag(A) transform
  
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
  
  # Get estimated signal based on posterior mean parameters (MCMC)
  
  # Posterior means - MCMC
  nu_A = post_means$mean[1:4]
  nu_B = post_means$mean[5:6]
  nu_C = post_means$mean[7]
  
  nu <- list(nu_A, nu_B, nu_C)
  
  # Initial value - MCMC
  z0 <- post_means$mean[8:9]
  
  # Constant shift beta - MCMC
  beta <- post_means$mean[10:11]
  
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
  
  # Get rid of first column of z (same as times)
  out_z <- out_z[,-1]
  
  ## Hemodynamic model
  # hrf function (using the difference of two gammas - canonical)
  HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
  HRF_mu = function(mu,tp){
    convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
  }
  
  y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  
  # Apply constant shift as estimated by beta
  y_pred = sweep(y_pred, 2, beta, "+")
  
  pred_signal[[i]] <- y_pred
}

save(pred_signal, file = "Balloon_Simulation/Output/balloon_diagAB_nonzero_pred.RData")

###############################################################################

rm(list = ls())

############### Canonical data generating model ###########################
######### Diagonal of A only ####################

pred_signal <- list()

for (i in 1:50){
  # Load in posterior means (MCMC)
  load(sprintf("Balloon_Simulation/Output/canonical_post_means_diagA_rep%03d.RData", i))
  
  # Load in observed signals (data)
  load(sprintf("Balloon_Simulation/Data/canonical_sim_data_diagA_rep%03d.RData", i))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  
  rm(dat)
  
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
    
    B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
    for(i in 1:nrow(idxs$B_idxs)){
      B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
    }
    
    C = matrix(data=0, nrow = m, ncol = n_u)
    for(i in 1:nrow(idxs$C_idxs)){
      C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
    }
    
    return(list(A=A,B=B,C=C))
    
  } # imported params already have diag(A) transform
  
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
  
  # Get estimated signal based on posterior mean parameters (MCMC)
  
  # Posterior means - MCMC
  nu_A = post_means$mean[1:4]
  nu_B = post_means$mean[5]
  nu_C = post_means$mean[6]
  
  nu <- list(nu_A, nu_B, nu_C)
  
  # Initial value - MCMC
  z0 <- post_means$mean[7:8]
  
  # Constant shift beta - MCMC
  beta <- post_means$mean[9:10]
  
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
  
  # Get rid of first column of z (same as times)
  out_z <- out_z[,-1]
  
  ## Hemodynamic model 
  # hrf function (using the difference of two gammas - canonical)
  HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
  HRF_mu = function(mu,tp){
    convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
  }
  
  y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  
  # Apply constant shift as estimated by beta
  y_pred = sweep(y_pred, 2, beta, "+")
  
  pred_signal[[i]] <- y_pred
}

save(pred_signal, file = "Balloon_Simulation/Output/canonical_diagA_pred.RData")

rm(list = ls())

######### Diagonal of A and B ####################

pred_signal <- list()

for (i in 1:50){
  # Load in posterior means (MCMC)
  load(sprintf("Balloon_Simulation/Output/canonical_post_means_diagAB_rep%03d.RData", i))
  
  # Load in observed signals (data)
  load(sprintf("Balloon_Simulation/Data/canonical_sim_data_diagAB_rep%03d.RData", i))
  u <- rbind(0,dat$u)
  times <- c(0,dat$times)
  
  rm(dat)
  
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
    
    B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
    for(i in 1:nrow(idxs$B_idxs)){
      B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
    }
    
    C = matrix(data=0, nrow = m, ncol = n_u)
    for(i in 1:nrow(idxs$C_idxs)){
      C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
    }
    
    return(list(A=A,B=B,C=C))
    
  } # imported params already have diag(A) transform
  
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
  
  # Get estimated signal based on posterior mean parameters (MCMC)
  
  # Posterior means - MCMC
  nu_A = post_means$mean[1:4]
  nu_B = post_means$mean[5:6]
  nu_C = post_means$mean[7]
  
  nu <- list(nu_A, nu_B, nu_C)
  
  # Initial value - MCMC
  z0 <- post_means$mean[8:9]
  
  # Constant shift beta - MCMC
  beta <- post_means$mean[10:11]
  
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
  
  # Get rid of first column of z (same as times)
  out_z <- out_z[,-1]
  
  ## Hemodynamic model 
  # hrf function (using the difference of two gammas - canonical)
  HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
  HRF_mu = function(mu,tp){
    convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
  }
  
  y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  
  # Apply constant shift as estimated by beta
  y_pred = sweep(y_pred, 2, beta, "+")
  
  pred_signal[[i]] <- y_pred
}

save(pred_signal, file = "Balloon_Simulation/Output/canonical_diagAB_pred.RData")

###############################################################################

