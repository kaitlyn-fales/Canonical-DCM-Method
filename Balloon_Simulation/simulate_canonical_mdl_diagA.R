remove(list=ls())

# Packages
library(expm)
library(posterior)
library(R.matlab)
library(cdcm)

# Generating inputs
##########################
m = 2; n_u = 2; max_time = 150*2
nu = list(nu_A = c(0.4,0.3,-0.1,0.15), nu_B = c(-0.2), nu_C=c(0.7)) 
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))
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

par(mfrow=c(1,2))
for (i in 1:m){
  plot(dat$times,dat$y_obs[,i],col=1,
       type='l',lty=2,xlab = "Time",ylab = "Signal")
}
mtext(paste("Simulated Noisy Observations at SNR = ",
            SNR,sep = ""), side = 3, line = -2, outer = T)
##########################################

# Write observed signal to a MATLAB file for SPM
writeMat("Balloon_Simulation/Data/canonical_sim_data_diagA.mat",
         U = dat$u,
         Y = dat$y_obs)

# Export just U matrix for balloon sim in Matlab as well
writeMat("Balloon_Simulation/Data/U_mat.mat",
         U = dat$u)

#########################################

# Simulate data using cdcm, at high SNR to get true curve to use for MSE comparison
sim <- simulate_cdcm(
  nu = nu,
  hypothesis_idxs = idxs,
  nscan = 150,
  m = m,
  n_u = n_u,
  TR = 2,
  U = input,
  SNR = 10000,
  seed = c(1:m),
  z0 = rep(0.1,m)
)

true_signal <- sim$simulated_data$y_obs

save(true_signal, file = "Balloon_Simulation/Data/canonical_sim_signal_diagA.RData")

# Oracle MSE
mse_V1 <- mean((dat$y_obs[,1] - true_signal[,1])^2)
mse_V2 <- mean((dat$y_obs[,2] - true_signal[,2])^2)






