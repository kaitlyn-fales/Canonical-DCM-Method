# Comparison of estimated signal using posterior means from each approach with observed signal

library(deSolve)
library(R.matlab)

# Load in posterior means (MCMC)
MCMC <- read.csv("Attention_Motion/Output/MCMC_summary.csv")
MCMC_means <- MCMC$mean

# Load in posterior means - initial condition (MCMC)
load("Attention_Motion/Output/MCMC_z0_mean.RData")

# Load in posterior means - diag(A) (MCMC)
load("Attention_Motion/Output/MCMC_diag_A.RData")

# Load in beta - (MCMC)
load("Attention_Motion/Output/beta.RData")

# Load in predicted signals (SPM)
y_pred_SPM <- readMat("Attention_Motion/SPM_VBA_Comparison/y_pred_SPM.mat")
y_pred_SPM <- y_pred_SPM$y.pred

# Load in observed signals (data)
load("Attention_Motion/motion_dat.RData")
u <- rbind(0,motion_dat$u)
times <- c(0,motion_dat$times)
y_obs <- motion_dat$y_obs

# 3 nodes, 3 experimental inputs
m = 3; n_u = 3

# Indices of parameters
A_idxs <- matrix(c(1,2,
                   2,1,
                   2,3,
                   3,2,
                   1,1,
                   2,2,
                   3,3), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,2,1,
                   3,2,1), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

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
    approxfun(x = times,y = u[,2],rule = 2)(t),
    approxfun(x = times,y = u[,3],rule = 2)(t))
}

##########################################################################
# Get estimated signal based on posterior mean parameters (MCMC)

# Posterior means - MCMC
nu_A = c(MCMC_means[1:4],diag_A$mean)
nu_B = MCMC_means[5:6]
nu_C = MCMC_means[7]

nu <- list(nu_A, nu_B, nu_C)

# Initial value - MCMC
z0 <- z0_mean$mean 

# Constant shift beta - MCMC
beta <- beta$mean

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

# Convolve with HRF
y_pred = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))

# Apply constant shift as estimated by beta to the data - not predicted signal for ease of viewing
y_obs = sweep(y_obs, 2, beta, "-")
##########################################################################


titles <- c("V1","V5","SPC")
xlabs <- c("Time (Seconds)","Time (Seconds)","Time (Seconds)")

par(mfrow = c(4,1), mar = c(4.5,6,4.5,1.5))
for (i in 1:ncol(y_pred)){
  plot(times[-1],y_obs[,i],type = 'l',xlab = xlabs[i],
       ylab = "BOLD Response", 
       ylim = c(min(c(y_obs[,i],y_pred[,i],y_pred_SPM[,i])),
                max(c(y_obs[,i],y_pred[,i],y_pred_SPM[,i]))), col = "black",
       main = titles[i], lty = "dotted",
       cex.axis = 1.3,    
       cex.lab = 1.4,    
       cex.main = 1.6,    
       font.axis = 2,      
       font.lab = 2,       
       font.main = 2)
  lines(times[-1],y_pred[,i],type = 'l',col = "blue")
  lines(times[-1],y_pred_SPM[,i],type = 'l',col = "red")
}
# Plot experimental design
plot(1, type="n", xlab = "Time (Seconds)", ylab="", yaxt = 'n', frame.plot = F,
     xlim=c(0,times[361]), ylim=c(0,1), main = "Experimental Design",
     cex.axis = 1.3,    
     cex.lab = 1.4,    
     cex.main = 1.6,    
     font.axis = 2,      
     font.lab = 2,       
     font.main = 2)
lines(x = times, y = u[,1], type = 'l',col = 'black', lty = 3)
lines(x = times, y = u[,2]/2, type = 'l',col = 'black', lty = 2)
lines(x = times, y = u[,3]/4, type = 'l',col = 'black')

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
  V1_MCMC = bootstrap_mse(y_obs[,1], y_pred[,1]),
  V1_SPM  = bootstrap_mse(y_obs[,1], y_pred_SPM[,1]),
  V5_MCMC = bootstrap_mse(y_obs[,2], y_pred[,2]),
  V5_SPM  = bootstrap_mse(y_obs[,2], y_pred_SPM[,2]),
  SPC_MCMC = bootstrap_mse(y_obs[,3], y_pred[,3]),
  SPC_SPM  = bootstrap_mse(y_obs[,3], y_pred_SPM[,3])
)

# Combine all into one tidy data frame
results_df <- do.call(rbind, results)
results_df




