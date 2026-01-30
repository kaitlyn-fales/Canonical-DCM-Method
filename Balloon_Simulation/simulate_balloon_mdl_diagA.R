remove(list=ls())

# Generating inputs
##########################
m=2; n_u=2; max_time=300
nu = list(nu_A = c(0.4,0.3,-0.1,0.15), nu_B = c(-0.2), nu_C=c(0.7),
          nu_kappa = c(0,0.05), nu_tao = c(-0.05,0), nu_epsilon = c(0,0.02))
#######################################
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))
times = NULL
input = NULL

dcm_balloon <- function(t, x, params, input) {
  # Break x down into corresponding parts for the five hidden states: z, s, f_in, v, q
  z <- x[1:m]
  s <- x[(m+1):(2*m)]
  f_in <- x[(2*m + 1):(3*m)]
  v <- x[(3*m + 1):(4*m)]
  q <- x[(4*m + 1):(5*m)]
  
  # t is the current time point in the integration 
  # z: neural signal
  A = params[["A"]]
  B = params[["B"]] 
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+", lapply(1:nrow(u), function(i) u[i,1]*B[[i]])) 
  A0 <- (A+B_all)
  # Self regulate diagonal of A and B as in SPM12
  diag(A0) <- -0.5*exp(diag(A0))
  dz <- A0%*%z + C%*%u
  
  # Set hemodynamic parameters
  H = params[["H"]]
  nu_kappa = params$nu_kappa
  nu_tao = params$nu_tao
  nu_epsilon = params$nu_epsilon
  
  # Exponentiate hemodynamic state variables (all except s)
  f_in <- exp(f_in)
  v <- exp(v)
  q <- exp(q)
  
  # s: vascular signal
  sd <- H['kappa']*exp(nu_kappa)
  ds <- z - sd*s - H['gamma']*(f_in - 1)
  ds <- as.matrix(ds)
  
  # ln(f_in): log of blood inflow
  dlnf <- s / f_in
  dlnf <- as.matrix(dlnf)
  
  # ln(v): log of venous volume
  tt <- H['tao']*exp(nu_tao) 
  f_out <- v^(1/H['alpha'])
  dlnv <- (f_in - f_out) / (tt*v)
  dlnv <- as.matrix(dlnv)
  
  # ln(q): log of dHb
  f_e <- (1 - (1 - H['E0'])^(1 / f_in)) / H['E0']
  dlnq <- (f_e*f_in - (f_out*q)/v) / (tt*q)
  dlnq <- as.matrix(dlnq)
  
  return(list(rbind(dz,ds,dlnf,dlnv,dlnq)))
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
  
  B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
  for(i in 1:nrow(idxs$B_idxs)){
    B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
  }
  
  C = matrix(data=0, nrow = m, ncol = n_u)
  for(i in 1:nrow(idxs$C_idxs)){
    C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
  }
  
  # Hemodynamic parameters (prespecified)
  H <- c(kappa = 0.64, 
         gamma = 0.32, 
         tao = 2, 
         alpha = 0.32, 
         E0 = 0.4)
  
  # Observation parameters (to be estimated)
  nu_kappa = nu$nu_kappa
  nu_tao = nu$nu_tao
  nu_epsilon = nu$nu_epsilon
  
  return(list(A=A,B=B,C=C,H=H,nu_kappa=nu_kappa,nu_tao=nu_tao,
              nu_epsilon=nu_epsilon))
  
} # modified for balloon model


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

# Forward solve system numerically
library(deSolve)
numeric_sol <- ode(y = c(rep(0.1,m),rep(0,4*m)),
                   times = times,
                   func = dcm_balloon,
                   parms= with(paramMats, list(A=A, B=B, C=C,
                                               H=H, nu_kappa=nu_kappa,
                                               nu_tao=nu_tao,
                                               nu_epsilon=nu_epsilon)),
                   input = input_u,
                   atol = 1e-6, 
                   rtol = 1e-6
)
output <- numeric_sol[-1,-1]

par(mfrow = c(3, 2))
plot(times[-1], output[, "1"], type = "l", main = "Neural state (z)", xlab = "Time", ylab = "")
lines(times[-1], output[, "2"], col = "blue")

plot(times[-1], output[, "3"], type = "l", main = "Signal (s)", xlab = "Time", ylab = "")
lines(times[-1], output[, "4"], col = "blue")

plot(times[-1], exp(output[, "5"]), type = "l", main = "Blood inflow (f_in)", xlab = "Time", ylab = "")
lines(times[-1], exp(output[, "6"]), col = "blue")

plot(times[-1], exp(output[, "7"]), type = "l", main = "Venous volume (v)", xlab = "Time", ylab = "")
lines(times[-1], exp(output[, "8"]), col = "blue")

plot(times[-1], exp(output[, "9"]), type = "l", main = "dHb (q)", xlab = "Time", ylab = "")
lines(times[-1], exp(output[, "10"]), col = "blue")

####### BOLD Signal Model #################################################
# Combine v and q output nonlinearly using the BOLD model to get final output signal

# Set prespecified parameters (SPM12)
TE <- 0.04              # echo time
V0 <- 4                 # resting venous volume (%)
r0 <- 25                # slope of intravascular relaxation rate R_iv as function of oxygen saturation S
nu0 <- 40.3             # frequency offset at the outer surface of magnetized vessels (Hz)
E0 <- paramMats$H['E0'] # resting oxygen extraction fraction

# Free param: estimated ratio of intra- to extra-vascular signal (can be region-specific)
ep <- exp(paramMats$nu_epsilon)

# Coefficients in BOLD signal model
k1 <- 4.3*nu0*E0*TE
k2 <- ep*r0*E0*TE
k3 <- 1 - ep
            
# Output equation of BOLD signal model
v <- exp(output[,(3*m + 1):(4*m)]) 
q <- exp(output[,(4*m + 1):(5*m)]) 
BOLD <- V0*(k1 - k1*q + k2 - k2*q/v + k3 - k3*v)
###########################################

# Plots for output to see results
par(mfrow = c(1,2))
plot(times[-1],BOLD[,1],type = 'l',ylim = c(min(BOLD),max(BOLD)))
plot(times[-1],BOLD[,2],type = 'l',ylim = c(min(BOLD),max(BOLD)))

##########################

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
#############################


# Save/export data
dat <- list(times = times[-1], u = u, y_obs = y_obs)
save(dat, file = "Balloon_Simulation/Data/balloon_sim_data_diagA.RData")

# Write observed signal to a MATLAB file for SPM
library(R.matlab)
writeMat("Balloon_Simulation/Data/balloon_sim_data_diagA.mat",
         U = dat$u,
         Y = dat$y_obs)



# Oracle MSE - also scale true signal 
BOLD <- BOLD * scale
mse_V1 <- mean((y_obs[,1] - BOLD[,1])^2)
mse_V2 <- mean((y_obs[,2] - BOLD[,2])^2)



