remove(list=ls())

library(R.matlab)

# Generating inputs - hemodynamic params = 0
##########################
m=2; n_u=2; max_time=300
times = NULL
if(is.null(times)){times <- seq(0, max_time, 2)}

# Load in SPM-generated balloon signal
SPM <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagAB_zero.mat")
BOLD <- SPM$y.signal
u <- SPM$U

# Add some white Gaussian noise according to SNR
SNR = 1.679376
variance <-  numeric()
y_obs <- matrix(NA, nrow = length(times)-1, ncol = m)
for (j in 1:m){
  set.seed(j)
  variance[j] <- (var(BOLD[,j])+(mean(BOLD[,j]))^2)/SNR
  y_obs[,j] <- BOLD[,j] + rnorm(length(times)-1, mean = 0, sd = sqrt(variance[j]))
}

variance
sqrt(variance)

# Plot before scaling
par(mfrow=c(1,2))
plot(times[-1],y_obs[,1], ylim=c(min(y_obs),max(y_obs)),col=4,type='l',lty=2,xlab = "Time",ylab = "Signal")
lines(times[-1],BOLD[,1],col=4)
plot(times[-1],y_obs[,2],ylim=c(min(y_obs),max(y_obs)),col=6,type='l',lty=2,xlab = "Time",ylab = "Signal")
lines(times[-1],BOLD[,2],col=6)
mtext(paste("Convolved True Signal (Solid) and Noisy Observations (Dashed) at SNR = ",SNR,sep = ""), side = 3, line = -2, outer = T)

# Enforce scaling like SPM
scale <- max(y_obs) - min(y_obs)
scale <- 4 / max(scale, 4)
y_obs <- y_obs * scale

# Plot after scaling
par(mfrow=c(1,2))
plot(times[-1],y_obs[,1], ylim=c(min(y_obs),max(y_obs)),col=4,type='l',lty=2,xlab = "Time",ylab = "Signal")
plot(times[-1],y_obs[,2],ylim=c(min(y_obs),max(y_obs)),col=6,type='l',lty=2,xlab = "Time",ylab = "Signal")
mtext(paste("Scaled Noisy Observations using SPM12 Hemodynamic Model at SNR = ",SNR,sep = ""), side = 3, line = -2, outer = T)

# Save/export data
dat <- list(times = times[-1], u = u, y_obs = y_obs, true_signal = (BOLD * scale))
save(dat, file = "Balloon_Simulation/Data/balloon_sim_data_diagAB_zero.RData")

# Write observed signal to a MATLAB file for SPM
writeMat("Balloon_Simulation/Data/balloon_sim_data_diagAB_zero.mat",
         U = dat$u,
         Y = dat$y_obs)

##########################

remove(list=ls())

# Generating inputs - hemodynamic params nonzero
##########################
m=2; n_u=2; max_time=300
times = NULL
if(is.null(times)){times <- seq(0, max_time, 2)}

# Load in SPM-generated balloon signal
SPM <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagAB_nonzero.mat")
BOLD <- SPM$y.signal
u <- SPM$U

# Add some white Gaussian noise according to SNR
SNR = 1.679376
variance <-  numeric()
y_obs <- matrix(NA, nrow = length(times)-1, ncol = m)
for (j in 1:m){
  set.seed(j)
  variance[j] <- (var(BOLD[,j])+(mean(BOLD[,j]))^2)/SNR
  y_obs[,j] <- BOLD[,j] + rnorm(length(times)-1, mean = 0, sd = sqrt(variance[j]))
}

variance
sqrt(variance)

# Plot before scaling
par(mfrow=c(1,2))
plot(times[-1],y_obs[,1], ylim=c(min(y_obs),max(y_obs)),col=4,type='l',lty=2,xlab = "Time",ylab = "Signal")
lines(times[-1],BOLD[,1],col=4)
plot(times[-1],y_obs[,2],ylim=c(min(y_obs),max(y_obs)),col=6,type='l',lty=2,xlab = "Time",ylab = "Signal")
lines(times[-1],BOLD[,2],col=6)
mtext(paste("Convolved True Signal (Solid) and Noisy Observations (Dashed) at SNR = ",SNR,sep = ""), side = 3, line = -2, outer = T)

# Enforce scaling like SPM
scale <- max(y_obs) - min(y_obs)
scale <- 4 / max(scale, 4)
y_obs <- y_obs * scale

# Plot after scaling
par(mfrow=c(1,2))
plot(times[-1],y_obs[,1], ylim=c(min(y_obs),max(y_obs)),col=4,type='l',lty=2,xlab = "Time",ylab = "Signal")
plot(times[-1],y_obs[,2],ylim=c(min(y_obs),max(y_obs)),col=6,type='l',lty=2,xlab = "Time",ylab = "Signal")
mtext(paste("Scaled Noisy Observations using SPM12 Hemodynamic Model at SNR = ",SNR,sep = ""), side = 3, line = -2, outer = T)

# Save/export data
dat <- list(times = times[-1], u = u, y_obs = y_obs, true_signal = (BOLD * scale))
save(dat, file = "Balloon_Simulation/Data/balloon_sim_data_diagAB_nonzero.RData")

# Write observed signal to a MATLAB file for SPM
writeMat("Balloon_Simulation/Data/balloon_sim_data_diagAB_nonzero.mat",
         U = dat$u,
         Y = dat$y_obs)

##########################


