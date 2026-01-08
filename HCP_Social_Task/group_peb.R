# PEB for combining subjects into group DCM

subjects <- c(100408,101107,101309,101915,103111,103414,103818,105014,
               105115,106016,108828,110411,111312,111716,113619,113922,114419,
               115320,116524,117122,118528,118730,118932,120111,122317,122620,
               123117,123925,124422,125525,126325,127630,127933,128127,128632,
               129028,130013,130316,131217,131722,133019,135225,135932,136833,
               138534,139637,140925,144832,146432,147737,148840,149337,149539,
               149741,151526,151627,153025,154734,156637,159340,160123,161731,
               162733,163129,176542,178950,188347,189450,190031,192540,196750,
               198451,199655,201111,208226,211417,211720,212318,214423,221319,
               239944,245333,280739,298051,366446,397760,414229,499566,654754,
               672756,751348,756055,792564,856766,857263,899885)

phase <- "RL"
mask <- "L"

# Extract posterior means and variances into list

# Empty lists
theta_hat <- list()
Sigma_i <- list()

for (i in 1:length(subjects)){
  load(paste0("Output/sub-",subjects[i],"_phase",phase,"_mask",mask,"_draws.RData"))
  
  summary <- posterior::summarise_draws(draws_df)[-1,]
  
  # Extract posterior means
  posterior_means <- summary$mean
  
  # Extract posterior standard deviations
  posterior_sds <- diag(summary$sd^2)
  
  # Add into master list
  theta_hat[[i]] <- posterior_means
  Sigma_i[[i]] <- posterior_sds
}

# N subjects, p parameters per subject
N <- length(theta_hat) 
p <- length(theta_hat[[1]])

# Stack posterior means into a matrix (p x N)
Theta <- do.call(cbind, theta_hat)  # columns = subjects

# Stack covariance matrices into an array (p x p x N)
Sigma_array <- array(0, dim = c(p, p, N))
for (i in 1:N) Sigma_array[,,i] <- Sigma_i[[i]]

library(Matrix)

# Initialize
mu <- rep(0, p)
Sigma_theta <- diag(1, p)

tol <- 1e-6
max_iter <- 100

for (iter in 1:max_iter) {
  
  # Precompute inverses of subject covariances
  Sigma_inv_list <- lapply(1:N, function(i) solve(Sigma_array[,,i]))
  
  # Step 1: compute group precision-weighted mean
  Sigma_post_inv <- solve(Sigma_theta) + Reduce("+", Sigma_inv_list)
  mu_post_num <- Reduce("+", lapply(1:N, function(i) Sigma_inv_list[[i]] %*% Theta[,i]))
  mu_post <- solve(Sigma_post_inv, mu_post_num)
  
  # Step 2: update between-subject covariance
  Sigma_theta_new <- matrix(0, p, p)
  for (i in 1:N) {
    diff <- Theta[,i] - mu_post
    Sigma_theta_new <- Sigma_theta_new + (diff %*% t(diff) + Sigma_array[,,i])
  }
  Sigma_theta_new <- Sigma_theta_new / N
  
  # Check convergence
  if (max(abs(Sigma_theta_new - Sigma_theta)) < tol) break
  
  mu <- mu_post
  Sigma_theta <- Sigma_theta_new
}


print(data.frame(parameter = summary$variable, 
                 theta_group = round(mu,digits = 3),
                 sigma2_group = round(diag(Sigma_theta),digits = 3)))

