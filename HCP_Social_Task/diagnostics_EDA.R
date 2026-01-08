setwd("/storage/work/krf5429/HCP_Social_Task")

library(dplyr)
library(mcmcse)
library(bayesplot)

load("diagnostics_compilation.RData")

diagnostics_df$prop_divergence <- round(diagnostics_df$prop_divergence, digits = 2)
diagnostics_df$prop_max_treedepth <- round(diagnostics_df$prop_max_treedepth, digits = 2)

# Summary
summary(diagnostics_df)

# Summary by subject
subject_summary <- diagnostics_df %>%
  group_by(subject) %>%
  summarise(prop_converged = mean(ess_multi_good, na.rm = T)) %>%
  arrange(desc(prop_converged))

table(subject_summary$prop_converged)

# Phase mask summary
phase_mask_summary <- diagnostics_df %>%
  group_by(phase, mask) %>%
  summarise(prop_converged = mean(ess_multi_good, na.rm = T)) %>%
  arrange(desc(prop_converged))

print(phase_mask_summary)

# Use only diagnostic features for clustering
diag_features <- diagnostics_df %>%
  select(prop_max_treedepth, prop_divergence, e_bfmi, mean_accept_stat, 
         stepsize, mean_mcse, ess_multi)
diag_features[is.na(diag_features)] <- 0

# Scale features
diag_scaled <- scale(diag_features)

# Try k-means clustering (3 groups)
set.seed(123)
km <- kmeans(diag_scaled, centers = 3, nstart = 25)

diagnostics_df$cluster <- factor(km$cluster)

# Summarize clustering to get a sense of chains that should converge, might converge, and likely won't
cluster_summary <- diagnostics_df %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    prop_converged = mean(ess_multi_ok, na.rm = T),
    prop_ess_bulk_ok = mean(ess_bulk_ok),
    prop_ess_tail_ok = mean(ess_tail_ok),
    mean_divergence = mean(prop_divergence),
    mean_treedepth = mean(prop_max_treedepth),
    mean_accept = mean(mean_accept_stat),
    mean_bfmi = mean(e_bfmi)
  ) %>%
  arrange(desc(prop_converged))
cluster_summary

# Chains that converged with no diagnostic issues
converge_chains <- diagnostics_df %>% filter(cluster == 2)

# Chains that converged with some diagnostic issues (max treedepth), or had a larger proportion of divergences
struggling_chains <- diagnostics_df %>% filter(cluster == 3)

# Chains that really struggled and had high proportion of divergences
bad_chains <- diagnostics_df %>% filter(cluster == 1)

######## Examine a random sample of chain draws from each cluster to visually inspect convergence #########
set.seed(1234)

### Converging chains ###
# Set temporary dataframe for ease of use later
temp <- converge_chains
  
# Get sample
cluster1_samp <- sample(1:nrow(temp), 5, replace = F)
samp <- cluster1_samp

# Visual inspection to identify potential issues to flag for
i = 5
load(paste0("Output/sub-",temp$subject[samp[i]],
            "_phase",temp$phase[samp[i]],
            "_mask",temp$mask[samp[i]],"_draws.RData"))

mcmc_trace(draws_df)
mcmc_acf(draws_df)
###

### Possible struggling chains - look okay ###
# Set temporary dataframe for ease of use later
temp <- struggling_chains

# Go through all options
cluster2_samp <- sample(1:nrow(temp), 5, replace = F)
samp <- cluster2_samp

# Visual inspection to identify potential issues to flag for
i = 5
load(paste0("Output/sub-",temp$subject[samp[i]],
            "_phase",temp$phase[samp[i]],
            "_mask",temp$mask[samp[i]],"_draws.RData"))

mcmc_trace(draws_df)
mcmc_acf(draws_df)
###

### Possible struggling chains - ordered by divergence ###
# Set temporary dataframe for ease of use later
temp <- struggling_chains %>% arrange(desc(prop_divergence))

# Go through first 15 (prop divergence >= 0.12)
samp <- c(1:15)

# Visual inspection to identify potential issues to flag for
i = 6 # 6, 7, 8, 9, 12 bad, 11 weird
load(paste0("Output/sub-",temp$subject[samp[i]],
            "_phase",temp$phase[samp[i]],
            "_mask",temp$mask[samp[i]],"_draws.RData"))

load(paste0("Output/sub-",temp$subject[samp[i]],
            "_phase",temp$phase[samp[i]],
            "_mask",temp$mask[samp[i]],"_diagnostics.RData"))
diagnostic_draws <- do.call(rbind, all_diagnostics)
mcmc_parcoord(draws_df[,-1], np = diagnostic_draws)
mcmc_trace(draws_df)
mcmc_acf(draws_df)
###

### Bad chains ###
# Set temporary dataframe for ease of use later
temp <- bad_chains %>% arrange(desc(prop_divergence))

# Go through all
samp <- c(1:nrow(temp))

# Visual inspection
i = 13 # all look bad, except for one
load(paste0("Output/sub-",temp$subject[samp[i]],
            "_phase",temp$phase[samp[i]],
            "_mask",temp$mask[samp[i]],"_draws.RData"))


mcmc_trace(draws_df)
mcmc_acf(draws_df)
###

# Filter to just phase RL mask L, as largest proportion of converged chains (0.94)
df <- diagnostics_df %>% filter(phase == "RL" & mask == "L")

# See if any chains need to obviously be excluded, first sort by cluster
df_cluster <- df %>% arrange(cluster)

# Get rid of chain with high divergences (1), and any chains with multi_ESS < 1000
df <- df %>% filter(prop_divergence <= 0.1 & ess_multi > 1000) # leaves 97/100

paste(df$subject, collapse = ",")

#########################################################################################################
rm(list = ls())

# ACF check on all chains
draw_files <- list.files(
  paste0(getwd(),"/Output"),
  pattern = "_draws\\.RData$",
  full.names = TRUE
)

# Function for computing model ACF
compute_model_acf <- function(file_path, lag_max = 20) {
  e <- new.env()
  load(file_path, envir = e)
  obj_name <- ls(e)[1]
  draws <- e[[obj_name]]
  draws <- draws[,c(2:15)]
  
  # Compute acf at lags 1, 5, 10, and 20 for each parameter
  acf_vals <- map_dfr(draws, ~{
    ac <- acf(.x, plot = FALSE, lag.max = lag_max)
    tibble(
      acf1 = ac$acf[2],
      acf5 = ac$acf[6],
      acf10 = ac$acf[11],
      acf20 = ac$acf[21]
    )
  })
  
  # Take median across parameters
  summary_tbl <- acf_vals %>%
    summarize(
      median_acf1 = median(acf1, na.rm = TRUE),
      median_acf5 = median(acf5, na.rm = TRUE),
      median_acf10 = median(acf10, na.rm = TRUE),
      median_acf20 = median(acf20, na.rm = TRUE)
    )
  
  tibble(
    subject = sub(".*sub-(\\d+).*", "\\1", file_path),
    phase = sub(".*phase([A-Z]+).*", "\\1", file_path),
    mask = sub(".*mask([A-Z]+).*", "\\1", file_path),
    median_acf1 = summary_tbl$median_acf1,
    median_acf5 = summary_tbl$median_acf5,
    median_acf10 = summary_tbl$median_acf10,
    median_acf20 = summary_tbl$median_acf20,
    acf_ok = ifelse(median_acf10 <= 0.5 & median_acf20 <= 0.3, T, F),
    file = basename(file_path)
  )
}

# Apply to all files
acf_summary_models <- map_dfr(draw_files, compute_model_acf)

# Merge with other diagnostics
load("diagnostics_compilation.RData")
diagnostics_df <- merge(diagnostics_df, acf_summary_models, by = c("subject","phase","mask"))

# Check 2x2 table of ESS and ACF
table("ACF OK" = diagnostics_df$acf_ok, "Multivariate ESS OK" = diagnostics_df$ess_multi_ok)

# Filter for ACF OK = F and ESS = F
failing_chains <- diagnostics_df %>% filter(acf_ok == F & ess_multi_ok == F)

# Try thinning these chains and recheck
bad_files <- paste0(getwd(),"/Output/",failing_chains$file)

# Function to thin chain
thin_chain <- function(file_path, thin_factor = 10) {
  e <- new.env()
  load(file_path, envir = e)
  obj_name <- ls(e)[1]
  draws <- e[[obj_name]]
  
  # Thinning: keep every kth row
  thinned_draws <- draws[seq(1, nrow(draws), by = thin_factor), ]
  
  return(thinned_draws)
}

# Function to check thinned chains
check_thinned_chain <- function(thinned_draws) {
  
  thinned_draws <- thinned_draws[,c(2:15)]
  
  # Median ACF per parameter
  acf_vals <- map_dfr(thinned_draws, ~{
    ac <- acf(.x, plot = FALSE, lag.max = 20)
    tibble(
      acf1 = ac$acf[2],
      acf5 = ac$acf[6],
      acf10 = ac$acf[11],
      acf20 = ac$acf[21]
    )
  })
  
  acf_summary <- acf_vals %>%
    summarize(
      median_acf1 = median(acf1, na.rm = TRUE),
      median_acf5 = median(acf5, na.rm = TRUE),
      median_acf10 = median(acf10, na.rm = TRUE),
      median_acf20 = median(acf20, na.rm = TRUE)
    )
  
  # Multivariate ESS 
  ess_multi <- multiESS(thinned_draws)
  
  data.frame(acf_summary,ess_multi)
}

# Apply to all files
thinning_results <- map(bad_files, function(f) {
  thinned <- thin_chain(f)
  check_thinned_chain(thinned)
})
thinning_results <- do.call(rbind, thinning_results)
thinning_results <- cbind(failing_chains[,1:3],thinned_iter = failing_chains$iterations/10,
                          prop_divergence = round(failing_chains[,6], digits = 2),thinning_results)
thinning_results <- thinning_results %>% arrange(median_acf20)

# Filter for ACF OK = F and ESS = T
failing_chains_acf <- diagnostics_df %>% filter(acf_ok == F & ess_multi_ok == T)

# Try thinning these chains and recheck
bad_files <- paste0(getwd(),"/Output/",failing_chains_acf$file)

# Apply to all files
thinning_results_acf <- map(bad_files, function(f) {
  thinned <- thin_chain(f)
  check_thinned_chain(thinned)
})
thinning_results_acf <- do.call(rbind, thinning_results_acf)
thinning_results_acf <- cbind(failing_chains_acf[,1:3],thinned_iter = failing_chains_acf$iterations/10,
                              prop_divergence = round(failing_chains_acf[,6], digits = 2),thinning_results_acf)
thinning_results_acf <- thinning_results_acf %>% arrange(median_acf20)

############# Filter additional chains that might converge to run for longer ###########################

# Filter down to chains who have not converged based on multivariate ESS, but currently have <10% divergences (already done)
rerun_chains <- diagnostics_df %>% filter(ess_multi_ok == F, prop_divergence <= 0.1)

# Make a text file for the files that will continue to run in hopes of converging (already done)
rerun_chains$file_list <- paste0("/storage/work/krf5429/HCP_Social_Task/Data/sub-",
                                 rerun_chains$subject,"_phase",
                                 rerun_chains$phase,"_mask",
                                 rerun_chains$mask,".RData")
# Export (already done)
writeLines(rerun_chains$file_list, "file_list_cont.txt")

#########################################################################################################