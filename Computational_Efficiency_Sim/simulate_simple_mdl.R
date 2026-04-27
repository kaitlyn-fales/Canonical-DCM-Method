remove(list = ls())

library(cdcm)

# Simulation settings
n_reps <- 50  

nscan <- 150
TR <- 2
m <- 2
n_u <- 2
max_time <- nscan * TR
SNR <- 1.679376
z0 <- rep(0.1, m)

out_dir <- "Computational_Efficiency_Sim/Data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
  
  for (i in 1:nrow(time_points)) {
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

# Model specifications
models <- list(
  diag_fix = list(
    nu = list(
      nu_A = c(0.4, 0.3),
      nu_B = c(-0.2),
      nu_C = c(0.3)
    ),
    idxs = list(
      A_idxs = matrix(c(
        2, 1,
        1, 2
      ), byrow = TRUE, ncol = 2),
      B_idxs = matrix(c(2, 1, 2), byrow = TRUE, ncol = 3),
      C_idxs = matrix(c(1, 1), byrow = TRUE, ncol = 2)
    )
  ),
  
  diag_var = list(
    nu = list(
      nu_A = c(0.4, 0.3, rep(0.1, m)),
      nu_B = c(-0.2),
      nu_C = c(0.3)
    ),
    idxs = list(
      A_idxs = matrix(c(
        2, 1,
        1, 2,
        1, 1,
        2, 2
      ), byrow = TRUE, ncol = 2),
      B_idxs = matrix(c(2, 1, 2), byrow = TRUE, ncol = 3),
      C_idxs = matrix(c(1, 1), byrow = TRUE, ncol = 2)
    )
  )
)


# Simulation function
simulate_one_rep <- function(model_spec, rep_id, sim_seed_base) {
  
  # One seed per region for additive noise
  # Example for rep 1 and m = 2: 1001, 1002
  sim_seed <- sim_seed_base + seq_len(m)
  
  sim <- simulate_cdcm(
    nu = model_spec$nu,
    hypothesis_idxs = model_spec$idxs,
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
  
  return(dat)
}

# Generate and save replicates
set.seed(12345)

# One base seed per simulation realization
sim_seed_bases <- sample.int(1e7, size = n_reps)

for (model_name in names(models)) {
  
  for (rep_id in seq_len(n_reps)) {
    
    dat <- simulate_one_rep(
      model_spec = models[[model_name]],
      rep_id = rep_id,
      sim_seed_base = sim_seed_bases[rep_id]
    )
    
    save(dat, file = file.path(out_dir,sprintf("dat_simple_mdl_%s_rep%03d.RData", model_name, rep_id)
      )
    )
  }
}





