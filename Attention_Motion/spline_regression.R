# Estimation of SNR via spline regression using attention to motion dataset

remove(list=ls())

load("Attention_Motion/motion_dat.RData")

# Plot data
par(mfrow = c(3,1))
for (i in 1:3) plot(x = motion$times, y = motion$y_obs[,i], type = 'l', 
                    xlab = "Time", ylab = "BOLD Signal")


# Fit a smoothing spline to the BOLD signal in each ROI
SNR <- numeric()

for (i in 1:ncol(motion$y_obs)){
  dt <- data.frame(y = motion$y_obs[,i])
  par(mfrow = c(1,1))
  plot(dt$y)
  fit.ss = smooth.spline(x = dt$y)
  muhat = fitted(fit.ss)
  lines(muhat)
  SNR[i] <- sd(muhat)/sd(dt$y - muhat)
}

# Take the mean of the SNR and use as SNR for simulated data
mean(SNR)


