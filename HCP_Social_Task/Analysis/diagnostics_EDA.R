library(dplyr)
library(mcmcse)
library(bayesplot)
library(R.matlab)

load("HCP_Social_Task/Analysis/diagnostics_compilation.RData")

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
converge_chains <- diagnostics_df %>% filter(cluster == 1)

# Chains that converged with some diagnostic issues
struggling_chains <- diagnostics_df %>% filter(cluster == 3)

# Chains that really struggled and had high proportion of divergences
bad_chains <- diagnostics_df %>% filter(cluster == 2)

######## Examine a random sample of chain draws from each cluster to visually inspect convergence #########
set.seed(1234)

### Converging chains ###
# Set temporary dataframe for ease of use later
temp <- converge_chains
  
# Get sample
cluster1_samp <- sample(1:nrow(temp), 5, replace = F)
samp <- cluster1_samp

# Visual inspection to identify potential issues to flag for
i = 1
load(paste0("HCP_Social_Task/Output/sub-",temp$subject[samp[i]],
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
load(paste0("HCP_Social_Task/Output/sub-",temp$subject[samp[i]],
            "_phase",temp$phase[samp[i]],
            "_mask",temp$mask[samp[i]],"_draws.RData"))

mcmc_trace(draws_df)
mcmc_acf(draws_df)
###

### Bad chains ###
# Set temporary dataframe for ease of use later
temp <- bad_chains %>% arrange(prop_divergence)

# Go through all
samp <- c(1:nrow(temp))

# Visual inspection
i = 7
load(paste0("HCP_Social_Task/Output/sub-",temp$subject[samp[i]],
            "_phase",temp$phase[samp[i]],
            "_mask",temp$mask[samp[i]],"_draws.RData"))

mcmc_trace(draws_df)
mcmc_acf(draws_df)
###

rm(list = ls())

# Determine which chains are not converged, or have poor diagnostics and should not be included in later analysis
# Impractical to examine all chains, so provide a set of rules to follow
load("HCP_Social_Task/Analysis/diagnostics_compilation.RData")

diagnostics_df$prop_divergence <- round(diagnostics_df$prop_divergence, digits = 2)
diagnostics_df$prop_max_treedepth <- round(diagnostics_df$prop_max_treedepth, digits = 2)

# Set of T/F flags
diagnostics_df$divergence_ok <- ifelse(diagnostics_df$prop_divergence <= 0.25, T, F)
diagnostics_df$accept_stat_ok <- ifelse(diagnostics_df$mean_accept_stat >= 0.5, T, F)

# Convergence flag
diagnostics_df$converged <- ifelse(diagnostics_df$ess_multi_ok & 
                                     diagnostics_df$divergence_ok &
                                     diagnostics_df$accept_stat_ok, T, F)

summary(diagnostics_df)

# Summary by subject
subject_summary <- diagnostics_df %>%
  group_by(subject) %>%
  summarise(prop_converged = mean(converged)) %>%
  arrange(desc(prop_converged))

table(subject_summary$prop_converged)

# Phase mask summary
phase_mask_summary <- diagnostics_df %>%
  group_by(phase, mask) %>%
  summarise(prop_converged = mean(converged)) %>%
  arrange(desc(prop_converged))

print(phase_mask_summary)

save(diagnostics_df, file = "HCP_Social_Task/Analysis/diagnostics_compilation_final.RData")

# Subset for use in prep for group level PEB
diagnostics_df$condition <- paste0("phase",diagnostics_df$phase,"_mask",diagnostics_df$mask)

diagnostics <- diagnostics_df %>% select(subject,condition,converged) %>% filter(converged)

# Export as .mat file for PEB in SPM
save(diagnostics, file = "HCP_Social_Task/Analysis/MCMC_diagnostics.RData")

