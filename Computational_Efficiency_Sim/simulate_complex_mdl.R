# Simulating from the complex model with the diag(A) fixed to -0.5
remove(list=ls())

# Packages
library(expm)
library(posterior)

# Generating inputs
##########################
m = 6; n_u = 2; max_time = 150*2
nu = list(nu_A = c(rep(c(0.4,0.3,0.5,0.1),2),0.4,0.3), 
          nu_B = rep(-0.2,3), 
          nu_C = rep(0.3,3))
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
times = NULL
input = NULL

linear <- function(t, z, params, input) {
  # t is the current time point in the integration, 
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]] 
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+", lapply(1:nrow(u), function(i) u[i,1]*B[[i]])) 
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}
gen_input = function(times, start_time, end_time, duration,between_duration){
  input = rep(0,length(times))
  
  current_start_time = start_time
  time_points = matrix(nrow = 0,ncol=2)
  while(1){
    current_end_time = current_start_time + duration 
    if(current_end_time > end_time){current_end_time = end_time}
    time_points = rbind(time_points, c(current_start_time, current_end_time))
    if(current_end_time + between_duration < end_time){
      current_start_time = current_end_time + between_duration 
    }else{
      break
    }
  }
  for(i in 1:nrow(time_points)){
    input[time_points[i,1]<times & times <=time_points[i,2]] = 1
  }
  input
}
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


if(is.null(times)){times <- seq(0, max_time, 2)}

if(is.null(input)){
  input1  = gen_input(times,start_time = 5,end_time = max(times),duration = 20,between_duration=40)
  input2  = gen_input(times,start_time = 25,end_time = max(times),duration = 20,between_duration=40)
  input = cbind(input1, input2)
}else{
  if(length(times)!=nrow(input)){stop("length(input)!=nrow(input)")}
}

input_u = function(t){
  sapply(1:ncol(input), function(j) approxfun(x = times,y = input[,j], rule=2, method="const",f=0)(t))
}

n_u = length(input_u(0))
u = input[-1,]

# set up parameters
paramMats = struct_paramMats(m = m,n_u = n_u,idxs = idxs,nu = nu)

# Simulating data
# Solving piecewise analytically
##########################
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

# Scenario a: u(t) = 0
eval_u_a <- function(t, z0, A, prev_index){
  expm::expAtv(A = A, t = times[t]-times[prev_index], v = z0)$eAtv
}

# Scenario b: u(t) != 0
eval_u_b <- function(t, z0, A, B, C, input, prev_index){
  A0 = A + Reduce("+",lapply(1:ncol(input), function(i) input[t-1,i]*B[[i]]))
  b0 = as.numeric(C%*%input[t-1,])
  zstar = -solve(A0,b0)
  expm::expAtv(A = A0,t = times[t]-times[prev_index],v = z0-zstar)$eAtv + zstar
}

# Function to get piecewise analytic solution to ODE
get_AnalyticSol <- function(input, changes, times, parms, m, v0 = NULL){
  
  # Get initial values if not specified
  if (is.null(v0) == T) {v0 <- rep(0.1,m)}
  
  # Matrix to store solution
  sol <- matrix(NA, nrow = length(times), ncol = m)
  sol[1,] <- v0 # add initial values
  
  # Initialize at u(t)=0
  prev_index <- 1
  
  # For loop to go from one change to the next
  for (i in 1:length(changes)){
    
    # Denote index for which we will work up to
    index <- changes[i]
    
    # Update solution based on change_type
    if (sum(input[index-1,]) == 0){
      sol[(prev_index+1):index,] <- t(sapply((prev_index+1):index, eval_u_a, 
                                             z0 = sol[prev_index,], A = parms$A,
                                             prev_index = prev_index))
    } else {
      if(sum(input[index-1,]) != 0){
        sol[(prev_index+1):index,] <- t(sapply((prev_index+1):index, eval_u_b, 
                                               z0 = sol[prev_index,], A = parms$A, 
                                               B = parms$B, C = parms$C,
                                               input = input, prev_index = prev_index))
      }
    }
    
    # Update indices
    prev_index <- index 
    
  }
  return(sol)
}
##########################

# Running functions to get analytic solution
# Obtain dataframe of changes based on input
changes <- get_changepoints(input)

# Running final function
analytic_sol <- get_AnalyticSol(input, changes, times, parms = paramMats, m)

# Quickly check with forward solving ODE numerically - good
library(deSolve)
numeric_sol <- ode(y = rep(0.1,m),
                   times = times,
                   func = linear,
                   parms= with(paramMats, list(A=A, B=B, C=C)),
                   input = input_u,                   
                   atol = 1e-6, 
                   rtol = 1e-6
)

round(abs(analytic_sol-numeric_sol[,-1]),digits=4)

##########################

# convolution with HRF
# https://bic-berkeley.github.io/psych-214-fall-2016/convolution_background.html
HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
HRF_mu = function(mu,tp){
  convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
}

y_obs_ana = sapply(1:m, function(i) HRF_mu(analytic_sol[-1,i],times[-1]))

par(mfrow = c(3,2))
for (i in 1:m) plot(y_obs_ana[,i], type = 'l', col = i)
#############################


# Add some white Gaussian noise according to approximate SNR obtained from SPM motion data
SNR = 1.679376
variance <-  numeric()
y_obs <- matrix(NA, nrow = length(times)-1, ncol = m)
for (j in 1:m){
  set.seed(j)
  variance[j] <- (var(y_obs_ana[,j])+(mean(y_obs_ana[,j]))^2)/SNR
  y_obs[,j] <- y_obs_ana[,j] + rnorm(length(times)-1, mean = 0, sd = sqrt(variance[j]))
}

variance
sqrt(variance)

par(mfrow=c(3,2))
for (i in 1:m){
  plot(times[-1],y_obs[,i],col=1,
       type='l',lty=2,xlab = "Time",ylab = "Signal")
  lines(times[-1],y_obs_ana[,i],col=2)
}
mtext(paste("Convolved True Signal (Solid) and Noisy Observations (Dashed) at SNR = ",
            SNR,sep = ""), side = 3, line = -2, outer = T)

##########################################

# Save/export data
dat <- list(times = times[-1], u = u, y_obs = y_obs)
save(dat, file = "Computational_Efficiency_Sim/Data/dat_complex_mdl_diag_fix.RData")

############################
# Structural identifiability check
# (A1)-(A4) met
# A%*%s vectors are linearly independent (A5)
z0 <- rep(0.1,m)
mat <- cbind(z0,
             paramMats$A%*%z0,
             paramMats$A%*%paramMats$A%*%z0,
             paramMats$A%*%paramMats$A%*%paramMats$A%*%z0,
             paramMats$A%*%paramMats$A%*%paramMats$A%*%paramMats$A%*%z0,
             paramMats$A%*%paramMats$A%*%paramMats$A%*%paramMats$A%*%paramMats$A%*%z0)
solve(mat,rep(0,m))

# A has d distinct real eigenvalues (A6)
eigen(paramMats$A)$values

# Eigenvalues of matrix exponential of \Tilde(A) are all positive real (A7)
eigen(expm(paramMats$A + paramMats$B[[1]]))$values
eigen(expm(paramMats$A + paramMats$B[[2]]))$values

# Data blocks invertible for each condition (A8)
solve(cbind(dat$y_obs[3:9,],1))
solve(cbind(dat$y_obs[13:19,],1))


################################################################################
# Simulating from the complex model with the diag(A) allowed to vary
remove(list=ls())

# Packages
library(expm)
library(posterior)

# Generating inputs
##########################
m = 6; n_u = 2; max_time = 150*2
nu = list(nu_A = c(rep(c(0.4,0.3,0.5,0.1),2),0.4,0.3,rep(0.1,m)), 
          nu_B = rep(-0.2,3), 
          nu_C = rep(0.3,3))
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              3,2,
                              2,3,
                              4,3,
                              3,4,
                              5,4,
                              4,5,
                              6,5,
                              5,6,
                              1,1,
                              2,2,
                              3,3,
                              4,4,
                              5,5,
                              6,6), 
                            byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2,
                              2,3,4,
                              2,5,6), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1,
                              3,1,
                              5,1), byrow=T, ncol=2))
times = NULL
input = NULL

linear <- function(t, z, params, input) {
  # t is the current time point in the integration, 
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]] 
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+", lapply(1:nrow(u), function(i) u[i,1]*B[[i]])) 
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}
gen_input = function(times, start_time, end_time, duration,between_duration){
  input = rep(0,length(times))
  
  current_start_time = start_time
  time_points = matrix(nrow = 0,ncol=2)
  while(1){
    current_end_time = current_start_time + duration 
    if(current_end_time > end_time){current_end_time = end_time}
    time_points = rbind(time_points, c(current_start_time, current_end_time))
    if(current_end_time + between_duration < end_time){
      current_start_time = current_end_time + between_duration 
    }else{
      break
    }
  }
  for(i in 1:nrow(time_points)){
    input[time_points[i,1]<times & times <=time_points[i,2]] = 1
  }
  input
}
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


if(is.null(times)){times <- seq(0, max_time, 2)}

if(is.null(input)){
  input1  = gen_input(times,start_time = 5,end_time = max(times),duration = 20,between_duration=40)
  input2  = gen_input(times,start_time = 25,end_time = max(times),duration = 20,between_duration=40)
  input = cbind(input1, input2)
}else{
  if(length(times)!=nrow(input)){stop("length(input)!=nrow(input)")}
}

input_u = function(t){
  sapply(1:ncol(input), function(j) approxfun(x = times,y = input[,j], rule=2, method="const",f=0)(t))
}

n_u = length(input_u(0))
u = input[-1,]

# set up parameters
paramMats = struct_paramMats(m = m,n_u = n_u,idxs = idxs,nu = nu)

# Simulating data
# Solving piecewise analytically
##########################
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

# Scenario a: u(t) = 0
eval_u_a <- function(t, z0, A, prev_index){
  expm::expAtv(A = A, t = times[t]-times[prev_index], v = z0)$eAtv
}

# Scenario b: u(t) != 0
eval_u_b <- function(t, z0, A, B, C, input, prev_index){
  A0 = A + Reduce("+",lapply(1:ncol(input), function(i) input[t-1,i]*B[[i]]))
  b0 = as.numeric(C%*%input[t-1,])
  zstar = -solve(A0,b0)
  expm::expAtv(A = A0,t = times[t]-times[prev_index],v = z0-zstar)$eAtv + zstar
}

# Function to get piecewise analytic solution to ODE
get_AnalyticSol <- function(input, changes, times, parms, m, v0 = NULL){
  
  # Get initial values if not specified
  if (is.null(v0) == T) {v0 <- rep(0.1,m)}
  
  # Matrix to store solution
  sol <- matrix(NA, nrow = length(times), ncol = m)
  sol[1,] <- v0 # add initial values
  
  # Initialize at u(t)=0
  prev_index <- 1
  
  # For loop to go from one change to the next
  for (i in 1:length(changes)){
    
    # Denote index for which we will work up to
    index <- changes[i]
    
    # Update solution based on change_type
    if (sum(input[index-1,]) == 0){
      sol[(prev_index+1):index,] <- t(sapply((prev_index+1):index, eval_u_a, 
                                             z0 = sol[prev_index,], A = parms$A,
                                             prev_index = prev_index))
    } else {
      if(sum(input[index-1,]) != 0){
        sol[(prev_index+1):index,] <- t(sapply((prev_index+1):index, eval_u_b, 
                                               z0 = sol[prev_index,], A = parms$A, 
                                               B = parms$B, C = parms$C,
                                               input = input, prev_index = prev_index))
      }
    }
    
    # Update indices
    prev_index <- index 
    
  }
  return(sol)
}
##########################

# Running functions to get analytic solution
# Obtain dataframe of changes based on input
changes <- get_changepoints(input)

# Running final function
analytic_sol <- get_AnalyticSol(input, changes, times, parms = paramMats, m)

# Quickly check with forward solving ODE numerically - good
library(deSolve)
numeric_sol <- ode(y = rep(0.1,m),
                   times = times,
                   func = linear,
                   parms= with(paramMats, list(A=A, B=B, C=C)),
                   input = input_u,                   
                   atol = 1e-6, 
                   rtol = 1e-6
)

round(abs(analytic_sol-numeric_sol[,-1]),digits=4)

##########################

# convolution with HRF
# https://bic-berkeley.github.io/psych-214-fall-2016/convolution_background.html
HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
HRF_mu = function(mu,tp){
  convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
}

y_obs_ana = sapply(1:m, function(i) HRF_mu(analytic_sol[-1,i],times[-1]))

par(mfrow = c(3,2))
for (i in 1:m) plot(y_obs_ana[,i], type = 'l', col = i)
#############################


# Add some white Gaussian noise according to approximate SNR obtained from SPM motion data
SNR = 1.679376
variance <-  numeric()
y_obs <- matrix(NA, nrow = length(times)-1, ncol = m)
for (j in 1:m){
  set.seed(j)
  variance[j] <- (var(y_obs_ana[,j])+(mean(y_obs_ana[,j]))^2)/SNR
  y_obs[,j] <- y_obs_ana[,j] + rnorm(length(times)-1, mean = 0, sd = sqrt(variance[j]))
}

variance
sqrt(variance)

par(mfrow=c(3,2))
for (i in 1:m){
  plot(times[-1],y_obs[,i],col=1,
       type='l',lty=2,xlab = "Time",ylab = "Signal")
  lines(times[-1],y_obs_ana[,i],col=2)
}
mtext(paste("Convolved True Signal (Solid) and Noisy Observations (Dashed) at SNR = ",
            SNR,sep = ""), side = 3, line = -2, outer = T)

##########################################

# Save/export data
dat <- list(times = times[-1], u = u, y_obs = y_obs)
save(dat, file = "Computational_Efficiency_Sim/Data/dat_complex_mdl_diag_var.RData")

############################
# Structural identifiability check
# (A1)-(A4) met
# A%*%s vectors are linearly independent (A5)
z0 <- rep(0.1,m)
mat <- cbind(z0,
             paramMats$A%*%z0,
             paramMats$A%*%paramMats$A%*%z0,
             paramMats$A%*%paramMats$A%*%paramMats$A%*%z0,
             paramMats$A%*%paramMats$A%*%paramMats$A%*%paramMats$A%*%z0,
             paramMats$A%*%paramMats$A%*%paramMats$A%*%paramMats$A%*%paramMats$A%*%z0)
solve(mat,rep(0,m))

# A has d distinct real eigenvalues (A6)
eigen(paramMats$A)$values

# Eigenvalues of matrix exponential of \Tilde(A) are all positive real (A7)
eigen(expm(paramMats$A + paramMats$B[[1]]))$values
eigen(expm(paramMats$A + paramMats$B[[2]]))$values

# Data blocks invertible for each condition (A8)
solve(cbind(dat$y_obs[3:9,],1))
solve(cbind(dat$y_obs[13:19,],1))





