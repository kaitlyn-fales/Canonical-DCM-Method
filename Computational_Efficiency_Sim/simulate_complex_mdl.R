# Simulating from the complex model with the diag(A) fixed to -0.5
remove(list=ls())

# Packages
library(expm)
library(posterior)
library(cdcm)

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
input = NULL

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

times <- seq(0, max_time, 2)

if(is.null(input)){
  input1  = gen_input(times,start_time = 20,end_time = max(times),duration = 20,between_duration=40)
  input2  = gen_input(times,start_time = 40,end_time = max(times),duration = 20,between_duration=40)
  input = cbind(input1, input2)
}else{
  if(length(times)!=nrow(input)){stop("length(input)!=nrow(input)")}
}


# Add some white Gaussian noise according to approximate SNR obtained from SPM motion data
SNR = 1.679376

# Simulate data using cdcm
sim <- simulate_cdcm(
  nu = nu,
  hypothesis_idxs = idxs,
  nscan = 150,
  m = m,
  n_u = n_u,
  TR = 2,
  U = input,
  SNR = SNR,
  seed = c(1:m),
  z0 = rep(0.1,m)
)

dat <- sim$simulated_data

par(mfrow=c(3,2))
for (i in 1:m){
  plot(dat$times,dat$y_obs[,i],col=1,
       type='l',lty=2,xlab = "Time",ylab = "Signal")
}
mtext(paste("Simulated Noisy Observations at SNR = ",
            SNR,sep = ""), side = 3, line = -2, outer = T)

##########################################

# Save/export data
save(dat, file = "Computational_Efficiency_Sim/Data/dat_complex_mdl_diag_fix.RData")

############################


################################################################################
# Simulating from the complex model with the diag(A) allowed to vary
remove(list=ls())

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
input = NULL

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

times <- seq(0, max_time, 2)

if(is.null(input)){
  input1  = gen_input(times,start_time = 20,end_time = max(times),duration = 20,between_duration=40)
  input2  = gen_input(times,start_time = 40,end_time = max(times),duration = 20,between_duration=40)
  input = cbind(input1, input2)
}else{
  if(length(times)!=nrow(input)){stop("length(input)!=nrow(input)")}
}


# Add some white Gaussian noise according to approximate SNR obtained from SPM motion data
SNR = 1.679376

# Simulate data using cdcm
sim <- simulate_cdcm(
  nu = nu,
  hypothesis_idxs = idxs,
  nscan = 150,
  m = m,
  n_u = n_u,
  TR = 2,
  U = input,
  SNR = SNR,
  seed = c(1:m),
  z0 = rep(0.1,m)
)

dat <- sim$simulated_data

par(mfrow=c(3,2))
for (i in 1:m){
  plot(dat$times,dat$y_obs[,i],col=1,
       type='l',lty=2,xlab = "Time",ylab = "Signal")
}
mtext(paste("Simulated Noisy Observations at SNR = ",
            SNR,sep = ""), side = 3, line = -2, outer = T)

##########################################

# Save/export data
save(dat, file = "Computational_Efficiency_Sim/Data/dat_complex_mdl_diag_var.RData")

############################





