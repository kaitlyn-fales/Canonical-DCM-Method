# Examining results from SPM attention to motion data

setwd("/storage/work/krf5429/attention_motion")

# Packages
library(cmdstanr, lib.loc = "/storage/work/krf5429/R_packages")
library(gridExtra, lib.loc = "/storage/work/krf5429/R_packages")
library(grid, lib.loc = "/storage/work/krf5429/R_packages")
library(bayesplot, lib.loc = "/storage/work/krf5429/R_packages")
library(rstanarm, lib.loc = "/storage/work/krf5429/R_packages")
library(ggridges, lib.loc = "/storage/work/krf5429/R_packages")
library(ggplot2)
library(posterior)
library(abind)


# Load in file names
temp = list.files(pattern="\\.csv$")

# Load in draws
draws5000 <- read_cmdstan_csv(temp[4])
draws65000 <- read_cmdstan_csv(temp[5])

# Chain draws together
draws <- abind(draws5000$post_warmup_draws,draws65000$post_warmup_draws,along = 1)

# Summarize results
results <- posterior::summarize_draws(draws, 
                                      c(default_summary_measures(),
                                        default_mcse_measures()))
results_export <- cbind(results[5:11,c(1:2,6:7)],CI_length = abs(results$q95-results$q5)[5:11])
write.csv(results_export, file = "MCMC_summary.csv")

# To use for est vs obs - initial value means
z0_mean <- results[12:14,1:2]
save(z0_mean, file = "MCMC_z0_mean.RData")

# Trace plots
mcmc_trace(draws[,,5:11])

draws <- matrix(draws, ncol = 14, nrow = 70000)
colnames(draws) <- results$variable
mcmc_areas(draws, pars = results$variable[5:11], prob = 0.9, point_est = "mean")


########## Saving param from first 5000 sampling draws ############

csv_contents <- read_cmdstan_csv(temp[4])
results <- posterior::summarize_draws(csv_contents$post_warmup_draws, 
                                      c(default_summary_measures(),default_mcse_measures()))
results_export <- cbind(results[5:11,c(1:2,6:7)],CI_length = abs(results$q95-results$q5)[5:11])
write.csv(results_export, file = "MCMC_summary.csv")

# Save tuning parameters and final draw to reuse
init <- list(sigma = array(csv_contents$post_warmup_draws[5000,1,2:4]),
             nu_A = array(csv_contents$post_warmup_draws[5000,1,5:8]),
             nu_B = array(csv_contents$post_warmup_draws[5000,1,9:10]),
             nu_C = array(csv_contents$post_warmup_draws[5000,1,11]),
             z0 = array(csv_contents$post_warmup_draws[5000,1,12:14]))
step_size <- csv_contents$step_size$'1'
inv_metric <- csv_contents$inv_metric$'1'

motion_draws <- list(init = init, step_size = step_size, inv_metric = inv_metric)
save(motion_draws, file = "motion_draws.RData")
