# Examine scaling differences between SPM balloon model and canonical model
library(R.matlab)

# Load in simulated signals
load("Balloon_Simulation/Data/canonical_sim_signal_diagA.RData")
canonical_signal <- true_signal

SPM_signal <- readMat("Balloon_Simulation/Data/balloon_sim_signal_diagA_zero.mat")
balloon_signal <- SPM_signal$y.signal

# Plot them on top of each other
max_time = 150*2
times <- seq(0, max_time, 2)

par(mfrow = c(1,2))

# Region 1
plot(x = times[-1], y = balloon_signal[,1], type = "l", xlab = "Time (s)", ylab = "Simulated BOLD Signal",
     ylim = c(min(cbind(balloon_signal,canonical_signal)),max(cbind(balloon_signal,canonical_signal))*1.15),
     col = "black", lty = 1, lwd = 1.5, main = "Region 1")
lines(x = times[-1], y = canonical_signal[,1], type = "l", col = "blue", lty = 2, lwd = 1.5)
legend("topleft", legend = c("SPM", "Canonical DCM"), col = c("black","blue"), lty = c(1,2), horiz = T,
       bty = "n", cex = 0.9, lwd = c(1.5,1.5))

# Region 2
plot(x = times[-1], y = balloon_signal[,2], type = "l", xlab = "Time (s)", ylab = "Simulated BOLD Signal",
     ylim = c(min(cbind(balloon_signal,canonical_signal)),max(cbind(balloon_signal,canonical_signal))*1.15),
     col = "black", lty = 1, lwd = 1.5, main = "Region 2")
lines(x = times[-1], y = canonical_signal[,2], type = "l", col = "blue", lty = 2, lwd = 1.5)
legend("topleft", legend = c("SPM", "Canonical DCM"), col = c("black","blue"), lty = c(1,2), horiz = T,
       bty = "n", cex = 0.9, lwd = c(1.5,1.5))


# Function to estimate scale between simulated canonical and balloon signals
# using same neuronal parameters (A,B,C)
compute_scales <- function(y_1, y_2) {
  if(!all(dim(y_1) == dim(y_2)))
    stop("y_1 and y_2 must have identical dimensions")
  
  m <- ncol(y_1)
  gamma_ls   <- numeric(m)
  gamma_peak <- numeric(m)
  
  for (k in 1:m) {
    ys <- y_1[, k]
    yi <- y_2[, k]
    
    # 1) Least-squares scale (best multiplicative fit)
    denom <- sum(yi^2)
    gamma_ls[k] <- ifelse(denom == 0, NA, sum(yi * ys) / denom)
    
    # 2) Peak-ratio scale
    max_ys <- max(ys)
    max_yi <- max(yi)
    gamma_peak[k] <- ifelse(max_yi == 0, NA, max_ys / max_yi)
  }
  
  data.frame(
    region = 1:m,
    gamma_ls = gamma_ls,
    gamma_peak = gamma_peak
  )
}

scales <- compute_scales(y_1 = canonical_signal, y_2 = balloon_signal)
print(scales)

# Load in est C from canonical estimation of balloon data C[1,1] = 0.7 truth
load("Balloon_Simulation/Output/balloon_post_means_diagA_zero.RData")
est_C = post_means$mean[6]
est_C*scales$gamma_ls[1] # scaling using least-squares
est_C*scales$gamma_peak[1] # scaling using peaks

# Load in est C from SPM estimation of canonical data C[1,1] = 0.7 truth
SPM <- readMat("Balloon_Simulation/Output/canonical_diagA_out.mat")
est_C = SPM$Ep.C[SPM$Ep.C != 0]
est_C/scales$gamma_ls[1] # reverse scaling using least-squares
est_C/scales$gamma_peak[1] # reverse scaling using peaks
