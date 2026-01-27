# Comparison of estimated signal using posterior means from each approach with observed signal

library(R.matlab)

######### Diagonal of A only ####################

# Load in predicted signal (SPM)
SPM <- readMat("Balloon_Simulation/Output/canonical_diagA_out.mat")
y_pred <- SPM$DCM.y

# Load in observed signals (data)
dat <- readMat("Balloon_Simulation/Data/canonical_sim_data_diagA.mat")
times <- seq(0, 300, 2)
y_obs <- dat$Y

# Compute SPM predicted shift (beta) - not available from program directly
beta_offset <- colMeans(y_obs - y_pred)
y_pred <- sweep(y_pred, 2, beta_offset, "+")

titles <- c("V1","V2")
xlabs <- c("Time (Seconds)","Time (Seconds)")

par(mfrow = c(2,1), mar = c(4.5,4.5,2.5,1.5))
for (i in 1:ncol(y_pred)){
  plot(times[-1],y_obs[,i],type = 'l',xlab = xlabs[i],
       ylab = "BOLD Response", 
       ylim = c(min(c(y_obs,y_pred)),
                max(c(y_obs,y_pred))), col = "black",
       main = titles[i], lty = "dotted",
       cex.axis = 1,    
       cex.lab = 1,    
       cex.main = 1.4,    
       font.axis = 2,      
       font.lab = 2,       
       font.main = 2)
  lines(times[-1],y_pred[,i],type = 'l',col = "blue")
}

# Calculate MSE and boostrapped SE
bootstrap_mse <- function(y_obs, y_pred, B = 10000, seed = 1234) {
  stopifnot(length(y_obs) == length(y_pred))
  set.seed(seed)
  
  # point estimate
  mse <- mean((y_obs - y_pred)^2)
  
  # bootstrap replicates
  n <- length(y_obs)
  mse_boot <- replicate(B, {
    idx <- sample(seq_len(n), replace = TRUE)
    mean((y_obs[idx] - y_pred[idx])^2)
  })
  
  # standard error and CI
  se <- sd(mse_boot)
  ci <- quantile(mse_boot, c(0.025, 0.975))
  
  out <- data.frame(
    MSE = mse,
    SE = se,
    CI_lower = ci[1],
    CI_upper = ci[2]
  )
  
  out[] <- lapply(out, round, digits = 4)
  out
}

results <- list(
  V1 = bootstrap_mse(y_obs[,1], y_pred[,1]),
  V2 = bootstrap_mse(y_obs[,2], y_pred[,2])
)

# Combine all into one tidy data frame
results_df <- do.call(rbind, results)
results_df

rm(list = setdiff(ls(), "bootstrap_mse"))

######### Diagonal of A and B ####################

# Load in predicted signal (SPM)
SPM <- readMat("Balloon_Simulation/Output/canonical_diagAB_out.mat")
y_pred <- SPM$DCM.y

# Load in observed signals (data)
dat <- readMat("Balloon_Simulation/Data/canonical_sim_data_diagAB.mat")
times <- seq(0, 300, 2)
y_obs <- dat$Y

# Compute SPM predicted shift (beta) - not available from program directly
beta_offset <- colMeans(y_obs - y_pred)
y_pred <- sweep(y_pred, 2, beta_offset, "+")

titles <- c("V1","V2")
xlabs <- c("Time (Seconds)","Time (Seconds)")

par(mfrow = c(2,1), mar = c(4.5,4.5,2.5,1.5))
for (i in 1:ncol(y_pred)){
  plot(times[-1],y_obs[,i],type = 'l',xlab = xlabs[i],
       ylab = "BOLD Response", 
       ylim = c(min(c(y_obs,y_pred)),
                max(c(y_obs,y_pred))), col = "black",
       main = titles[i], lty = "dotted",
       cex.axis = 1,    
       cex.lab = 1,    
       cex.main = 1.4,    
       font.axis = 2,      
       font.lab = 2,       
       font.main = 2)
  lines(times[-1],y_pred[,i],type = 'l',col = "blue")
}

results <- list(
  V1 = bootstrap_mse(y_obs[,1], y_pred[,1]),
  V2 = bootstrap_mse(y_obs[,2], y_pred[,2])
)

# Combine all into one tidy data frame
results_df <- do.call(rbind, results)
results_df

