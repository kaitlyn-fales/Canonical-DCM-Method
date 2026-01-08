# Comparison of estimated signal using posterior means from each approach with observed signal

setwd("/storage/work/krf5429/attention_motion")

library(deSolve)
library(R.matlab, lib.loc = "/storage/work/krf5429/R_packages")

# Load in posterior means (MCMC)
MCMC <- read.csv("MCMC_summary.csv")
MCMC_means <- MCMC$mean

# Load in posterior means - initial condition (MCMC)
load("MCMC_z0_mean.RData")

# Load in predicted signals (SPM)
y_pred_SPM <- readMat("y_pred_SPM.mat")
y_pred_SPM <- y_pred_SPM$y.pred

# Load in raw observed signals (data pre-detrending and scaling)
#y_obs <- readMat("y_obs_raw.mat")
#y_obs <- cbind(y_obs$V1.BOLD,y_obs$V5.BOLD,y_obs$SPC.BOLD)

# Load in observed signals (data)
load("motion_dat.RData")
u <- rbind(0,motion_dat$u)
times <- c(0,motion_dat$times)
y_obs <- motion_dat$y_obs

##########################################################################
# Get estimated signal based on posterior mean parameters

# Posterior means - canonical
nu_A = MCMC_means[1:4]
nu_B = MCMC_means[5:6]
nu_C = MCMC_means[7]

##########################################################################

# Putting together to form U matrix of experimental inputs
input_u = function(t){
  c(approxfun(x = times,y = u[,1],rule = 2)(t),
    approxfun(x = times,y = u[,2],rule = 2)(t),
    approxfun(x = times,y = u[,3],rule = 2)(t))
}

# 3 nodes, 3 experimental inputs
m = 3; n_u = 3

# Introduce parameters
A = diag(x = rep(-0.5,m))
B = lapply(1:n_u, function(i){matrix(data=0,nrow=m, ncol=m)})
C = matrix(data = 0,nrow = m,ncol = n_u)

A[1,2] = nu_A[1]
A[2,1] = nu_A[2]
A[2,3] = nu_A[3]
A[3,2] = nu_A[4]
B[[2]][2,1] = nu_B[1]
B[[3]][2,1] = nu_B[2]
C[1,1] = nu_C

# ODE for neural activation that needs to be solved
linear <- function(t, z, params, input) {
  # t is the current time point in the integration, 
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]]
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+",lapply(1:nrow(u), function(i) u[i,1]*B[[i]]))
  #B_all = do.call("+", lapply(1:nrow(u), function(i) u[i,1]*B[[i]])) 
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}

# Initial value
z0 <- z0_mean$mean # initial neuronal activity of 0Hz

# Model params
params <- list(A = A, B = B, C = C)

# Solve the ODE for z
out_z <- ode(y = z0, times = times, func = linear, parms = params, input = input_u)
diagnostics(out_z)

# Get rid of first column of z (same as times)
out_z <- out_z[,-1]

# Plotting resulting z
par(mfrow = c(1,1))
plot(times,out_z[,1],type = 'l',xlab = "Time (seconds)",
     ylab = "z", main = "Neural Activation", ylim = c(-4,9))
lines(times,out_z[,2],type = 'l',col = 2)
lines(times,out_z[,3],type = 'l',col = 3)

######### Hemodynamic model #####################################
# hrf function (using the difference of two gammas - canonical)
hrf <- function(t, par = c(6,1,12,1,0.35)){
  result <- dgamma(x = t, shape = par[1],rate = par[2]) - 
    par[5]*dgamma(x = t,shape = par[3],rate = par[4])
  return(result)
}

# Convolution of neaural activation and hrf function
get_y <- function(z,time,par){
  convolution <- convolve(z, rev(hrf(time)), type = "open")[1:length(time)]
  return(convolution)
}

y_pred <- matrix(NA, nrow = length(times)-1, ncol = m)
for (i in 1:m){
  y_pred[,i] <- get_y(out_z[-1,i],times[-1])
}
##########################################################################


titles <- c("V1","V5","SPC")
xlabs <- c("","","Time (Seconds)")

par(mfrow = c(4,1), mar = c(2,4,4.5,1.5))
for (i in 1:ncol(y_pred)){
  plot(times[-1],y_obs[,i],type = 'l',xlab = xlabs[i],
       ylab = "BOLD Response", 
       ylim = c(min(c(y_obs[,i],y_pred[,i],y_pred_SPM[,i])),
                max(c(y_obs[,i],y_pred[,i],y_pred_SPM[,i]))), col = "black",
       main = titles[i], lty = "dotted")
  lines(times[-1],y_pred[,i],type = 'l',col = "blue")
  lines(times[-1],y_pred_SPM[,i],type = 'l',col = "red")
}
# Plot experimental design
plot(1, type="n", xlab = "Time (Seconds)", ylab="", yaxt = 'n', frame.plot = F,
     xlim=c(0,times[361]), ylim=c(0,1), main = "Experimental Design")
lines(x = times, y = u[,1], type = 'l',col = 'black', lty = 3)
lines(x = times, y = u[,2]/2, type = 'l',col = 'black', lty = 2)
lines(x = times, y = u[,3]/4, type = 'l',col = 'black')
mtext("Estimated vs. Observed Signal (y) Comparison Between MCMC and SPM", 
      line = -1.5, outer = T)

# V1
mean((y_obs[,1] - y_pred[,1])^2)
mean((y_obs[,1] - y_pred_SPM[,1])^2)

# V5
mean((y_obs[,2] - y_pred[,2])^2)
mean((y_obs[,2] - y_pred_SPM[,2])^2)

# SPC
mean((y_obs[,3] - y_pred[,3])^2)
mean((y_obs[,3] - y_pred_SPM[,3])^2)



