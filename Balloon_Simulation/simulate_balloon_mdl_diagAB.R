remove(list = ls())

library(R.matlab)

# Settings
n_reps <- 50

m <- 2
n_u <- 2
max_time <- 300
TR <- 2
times <- seq(0, max_time, TR)

SNR <- 1.679376
out_dir <- "Balloon_Simulation/Data"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function
safe_noise_scale <- function(BOLD, noise_mat, max_range = 4, tol = 1e-10) {
  
  range_with_scale <- function(alpha) {
    y_tmp <- BOLD + alpha * noise_mat
    max(y_tmp) - min(y_tmp)
  }
  
  if (range_with_scale(1) <= max_range - tol) {
    return(1)
  }
  
  lower <- 0
  upper <- 1
  
  for (i in seq_len(100)) {
    mid <- (lower + upper) / 2
    
    if (range_with_scale(mid) <= max_range - tol) {
      lower <- mid
    } else {
      upper <- mid
    }
  }
  
  lower
}

generate_balloon_reps <- function(delay_label, signal_file, out_prefix){
  
  # Load SPM-generated balloon signal
  SPM <- readMat(file.path(out_dir, signal_file))
  BOLD <- SPM$y.signal
  u <- SPM$U
  
  # Clean signal range
  R_signal <- max(BOLD) - min(BOLD)
  
  if (R_signal >= 4) {
    stop("Clean signal already has range >= 4; cannot add noise without exceeding range 4.")
  }
  
  set.seed(12345)
  sim_seed_bases <- sample.int(1e7, size = n_reps)
  
  for (rep_id in seq_len(n_reps)) {
    
    sim_seed <- sim_seed_bases[rep_id] + seq_len(m)
    
    noise_mat <- matrix(NA, nrow = length(times) - 1, ncol = m)
    variance <- numeric(m)
    
    for (k in seq_len(m)) {
      
      set.seed(sim_seed[k])
      
      # SNR is defined as sd(signal) / sd(noise)
      # Therefore Var(noise) = Var(signal) / SNR^2
      variance[k] <- stats::var(BOLD[, k]) / SNR^2
      
      noise_mat[, k] <- stats::rnorm(
        length(times) - 1,
        mean = 0,
        sd = sqrt(variance[k])
      )
    }
    
    # Scale additive noise if needed so that range(BOLD + noise) <= 4
    noisy_unscaled <- BOLD + noise_mat
    unscaled_range <- max(noisy_unscaled) - min(noisy_unscaled)
    
    noise_scale_factor <- safe_noise_scale(
      BOLD = BOLD,
      noise_mat = noise_mat,
      max_range = 4,
      tol = 1e-10
    )
    
    noise_mat_scaled <- noise_scale_factor * noise_mat
    y_obs <- BOLD + noise_mat_scaled
    
    final_range <- max(y_obs) - min(y_obs)
    
    if (final_range > 4) {
      stop(sprintf(
        "Replicate %d exceeded range 4 after noise scaling. Final range = %.12f",
        rep_id, final_range
      ))
    }
    
    dat <- list(times = times[-1],u = u,y_obs = y_obs)
    
    save(dat,file = file.path(out_dir,sprintf("%s_data_rep%03d.RData", out_prefix, rep_id)))
    
    writeMat(file.path(out_dir,sprintf("%s_data_rep%03d.mat", out_prefix, rep_id)),
             U = dat$u,
             Y = dat$y_obs
    )
  }
}


# Generate zero-delay replicates
summary_zero <- generate_balloon_reps(
  delay_label = "zero",
  signal_file = "balloon_sim_signal_diagAB_zero.mat",
  out_prefix = "balloon_sim_diagAB_zero"
)

# Generate nonzero-delay replicates
summary_nonzero <- generate_balloon_reps(
  delay_label = "nonzero",
  signal_file = "balloon_sim_signal_diagAB_nonzero.mat",
  out_prefix = "balloon_sim_diagAB_nonzero"
)


