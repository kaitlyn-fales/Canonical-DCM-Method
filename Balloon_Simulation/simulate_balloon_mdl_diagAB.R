remove(list=ls())

library(R.matlab)

# Generating inputs - delays = 0
##########################
m=2; n_u=2; max_time=300
times = NULL
if(is.null(times)){times <- seq(0, max_time, 2)}

# Load in SPM-generated balloon signal
SPM <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagAB_zero.mat")
BOLD <- SPM$y.signal
u <- SPM$U

SNR <- 1.679376 # target (real SNR will be larger than this to avoid triggering SPM internal scaling)
y_obs <- matrix(NA, nrow = length(times) - 1, ncol = m)
variance <- numeric()

# Add Gaussian noise but ensure range <= 4
R_signal <- max(BOLD) - min(BOLD)
max_allowed_noise_span <- 4 - R_signal
if (max_allowed_noise_span <= 0)
  stop("Clean signal already has range >= 4: cannot add noise without triggering scaling.")

for (k in 1:m) {
  # target variance for given SNR
  variance[k] <- (var(BOLD[, k]) + (mean(BOLD[, k]))^2) / SNR
  
  # draw Gaussian noise
  noise <- rnorm(length(times) - 1, mean = 0, sd = sqrt(variance[k]))
  
  # compute span of noise and scale it if it would exceed allowed range
  noise_span <- max(noise) - min(noise)
  if (noise_span > max_allowed_noise_span) {
    scale_factor <- max_allowed_noise_span / noise_span
    noise <- noise * scale_factor
  }
  
  # add noise to clean signal
  y_obs[, k] <- BOLD[, k] + noise
}

# Compute achieved SNR (mean across columns) 
# SNR = var(signal) / var(noise)
achieved_SNR <- numeric(m)
for (k in 1:m) {
  achieved_SNR[k] <- var(BOLD[, k]) / var(y_obs[, k] - BOLD[, k])
}

# Final safety scaling like SPM, with flag to indicate when needed
scale <- max(y_obs) - min(y_obs)
final_scale_factor <- 4 / max(scale, 4)
y_obs <- y_obs * final_scale_factor

final_scale_flag <- ifelse(final_scale_factor < 1, 1, 0)
final_scale_vec <- final_scale_factor

# Plot 
par(mfrow=c(1,2))
plot(times[-1],y_obs[,1], ylim=c(min(y_obs),max(y_obs)),col=4,type='l',lty=2,xlab = "Time",ylab = "Signal")
lines(times[-1],BOLD[,1],col=4)
plot(times[-1],y_obs[,2],ylim=c(min(y_obs),max(y_obs)),col=6,type='l',lty=2,xlab = "Time",ylab = "Signal")
lines(times[-1],BOLD[,2],col=6)
mtext(paste("Convolved True Signal (Solid) and Noisy Observations (Dashed) at SNR = ",SNR,sep = ""), side = 3, line = -2, outer = T)

# Save/export data
dat <- list(times = times[-1], u = u, y_obs = y_obs)
save(dat, file = "Balloon_Simulation/Data/balloon_sim_data_diagAB_zero.RData")

# Write observed signal to a MATLAB file for SPM
writeMat("Balloon_Simulation/Data/balloon_sim_data_diagAB_zero.mat",
         U = dat$u,
         Y = dat$y_obs)

##########################

remove(list=ls())

# Generating inputs - delays 1s
##########################
m=2; n_u=2; max_time=300
times = NULL
if(is.null(times)){times <- seq(0, max_time, 2)}

# Load in SPM-generated balloon signal
SPM <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagAB_nonzero.mat")
BOLD <- SPM$y.signal
u <- SPM$U

SNR <- 1.679376 # target (real SNR will be larger than this to avoid triggering SPM internal scaling)
y_obs <- matrix(NA, nrow = length(times) - 1, ncol = m)
variance <- numeric()

# Add Gaussian noise but ensure range <= 4
R_signal <- max(BOLD) - min(BOLD)
max_allowed_noise_span <- 4 - R_signal
if (max_allowed_noise_span <= 0)
  stop("Clean signal already has range >= 4: cannot add noise without triggering scaling.")

for (k in 1:m) {
  # target variance for given SNR
  variance[k] <- (var(BOLD[, k]) + (mean(BOLD[, k]))^2) / SNR
  
  # draw Gaussian noise
  noise <- rnorm(length(times) - 1, mean = 0, sd = sqrt(variance[k]))
  
  # compute span of noise and scale it if it would exceed allowed range
  noise_span <- max(noise) - min(noise)
  if (noise_span > max_allowed_noise_span) {
    scale_factor <- max_allowed_noise_span / noise_span
    noise <- noise * scale_factor
  }
  
  # add noise to clean signal
  y_obs[, k] <- BOLD[, k] + noise
}

# Compute achieved SNR (mean across columns) 
# SNR = var(signal) / var(noise)
achieved_SNR <- numeric(m)
for (k in 1:m) {
  achieved_SNR[k] <- var(BOLD[, k]) / var(y_obs[, k] - BOLD[, k])
}

# Final safety scaling like SPM, with flag to indicate when needed
scale <- max(y_obs) - min(y_obs)
final_scale_factor <- 4 / max(scale, 4)
y_obs <- y_obs * final_scale_factor

final_scale_flag <- ifelse(final_scale_factor < 1, 1, 0)
final_scale_vec <- final_scale_factor

# Plot 
par(mfrow=c(1,2))
plot(times[-1],y_obs[,1], ylim=c(min(y_obs),max(y_obs)),col=4,type='l',lty=2,xlab = "Time",ylab = "Signal")
lines(times[-1],BOLD[,1],col=4)
plot(times[-1],y_obs[,2],ylim=c(min(y_obs),max(y_obs)),col=6,type='l',lty=2,xlab = "Time",ylab = "Signal")
lines(times[-1],BOLD[,2],col=6)
mtext(paste("Convolved True Signal (Solid) and Noisy Observations (Dashed) at SNR = ",SNR,sep = ""), side = 3, line = -2, outer = T)

# Save/export data
dat <- list(times = times[-1], u = u, y_obs = y_obs)
save(dat, file = "Balloon_Simulation/Data/balloon_sim_data_diagAB_nonzero.RData")

# Write observed signal to a MATLAB file for SPM
writeMat("Balloon_Simulation/Data/balloon_sim_data_diagAB_nonzero.mat",
         U = dat$u,
         Y = dat$y_obs)

##########################


