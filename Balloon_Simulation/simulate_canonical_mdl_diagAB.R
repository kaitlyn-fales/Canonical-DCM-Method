remove(list = ls())

# Packages
library(R.matlab)
library(cdcm)

# Settings
n_reps <- 50

m <- 2
n_u <- 2
nscan <- 150
TR <- 2
max_time <- nscan * TR
SNR <- 1.679376
z0 <- rep(0.1, m)

out_dir <- "Balloon_Simulation/Data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Model setup
nu = list(nu_A = c(0.4,0.3,-0.1,0.15), 
          nu_B = c(-0.2,0.05), 
          nu_C=c(0.7)) 

idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2,
                              2,2,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

# Input generation
gen_input <- function(times, start_time, end_time, duration, between_duration) {
  input <- rep(0, length(times))
  
  current_start_time <- start_time
  time_points <- matrix(nrow = 0, ncol = 2)
  
  while (TRUE) {
    current_end_time <- current_start_time + duration
    
    if (current_end_time > end_time) {
      current_end_time <- end_time
    }
    
    time_points <- rbind(time_points, c(current_start_time, current_end_time))
    
    if (current_end_time + between_duration < end_time) {
      current_start_time <- current_end_time + between_duration
    } else {
      break
    }
  }
  
  for (i in seq_len(nrow(time_points))) {
    input[time_points[i, 1] < times & times <= time_points[i, 2]] <- 1
  }
  
  input
}

times <- seq(0, max_time, TR)

input1 <- gen_input(
  times = times,
  start_time = 20,
  end_time = max(times),
  duration = 20,
  between_duration = 40
)

input2 <- gen_input(
  times = times,
  start_time = 40,
  end_time = max(times),
  duration = 20,
  between_duration = 40
)

input <- cbind(input1, input2)

# Save true signal once
sim_true <- simulate_cdcm(
  nu = nu,
  hypothesis_idxs = idxs,
  nscan = nscan,
  m = m,
  n_u = n_u,
  TR = TR,
  U = input,
  SNR = 1e6,
  seed = rep(999, m),
  z0 = z0
)

true_signal <- sim_true$simulated_data$y_obs

save(
  true_signal,
  file = file.path(out_dir, "canonical_sim_signal_diagAB.RData")
)

# Generate noisy simulation replicates
set.seed(12345)
sim_seed_bases <- sample.int(1e7, size = n_reps)

for (rep_id in seq_len(n_reps)) {
  
  sim_seed <- sim_seed_bases[rep_id] + seq_len(m)
  
  sim <- simulate_cdcm(
    nu = nu,
    hypothesis_idxs = idxs,
    nscan = nscan,
    m = m,
    n_u = n_u,
    TR = TR,
    U = input,
    SNR = SNR,
    seed = sim_seed,
    z0 = z0
  )
  
  dat <- sim$simulated_data
  
  save(dat,file = file.path(out_dir,sprintf("canonical_sim_data_diagAB_rep%03d.RData", rep_id)))
  
  writeMat(file.path(out_dir,sprintf("canonical_sim_data_diagAB_rep%03d.mat", rep_id)),
           U = dat$u,
           Y = dat$y_obs
  )
}


