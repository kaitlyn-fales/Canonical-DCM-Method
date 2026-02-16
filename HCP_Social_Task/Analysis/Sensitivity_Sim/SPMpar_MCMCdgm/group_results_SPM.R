# Extract and summarize posterior draws from group level model for MCMC
library(posterior)
library(dplyr)
library(R.matlab)

# Number of parameters in group level model of interest
p <- 10

# Names from original notation
par_names <- c("nu_A[1]","nu_A[2]","nu_A[3]","nu_A[4]","nu_B[1]",
               "nu_B[2]","nu_B[3]","nu_B[4]","nu_C[1]","nu_C[2]")

# Conditions
conditions <- c("phaseLR_maskL","phaseLR_maskR",
                "phaseRL_maskL","phaseRL_maskR")

# Loop through conditions
for (i in 1:length(conditions)){
  load(paste0("HCP_Social_Task/Analysis/Sensitivity_Sim/SPMpar_MCMCdgm/Output/",conditions[i],"_SPM.RData"))
  
  # Extract draws for alpha[1] ... alpha[p]
  alpha_vars <- paste0("alpha[", 1:p, "]")
  alpha_draws <- subset_draws(draws, variable = alpha_vars)
  alpha_summary <- summarize_draws(alpha_draws, mean, sd, ~quantile(.x, probs = c(0.025, 0.975))) 
  
  # Extract and summarize tau draws
  tau_vars <- paste0("tau[", 1:p, "]")
  tau_draws <- subset_draws(draws, variable = tau_vars)
  tau_summary <- summarize_draws(tau_draws, mean, sd, ~quantile(.x, probs = c(0.025, 0.975))) 
  
  # General summary table
  summary_table <- tibble(
    parameter = par_names, # original parameter names
    mean_alpha = alpha_summary$mean,
    sd_alpha = alpha_summary$sd,    
    q2.5_alpha = alpha_summary$`2.5%`,
    q97.5_alpha = alpha_summary$`97.5%`,
    between_subject_sd = tau_summary$mean,
    between_subject_q2.5 = tau_summary$`2.5%`,
    between_subject_q97.5 = tau_summary$`97.5%`
  )
  
  assign(conditions[i], summary_table)
}

# Clear environment except for results
rm(list = setdiff(ls(), c("phaseLR_maskL","phaseLR_maskR",
                          "phaseRL_maskL","phaseRL_maskR")))

# Combine results into one list once done and export
SPM_results <- list(phaseLR_maskL = phaseLR_maskL,
                    phaseLR_maskR = phaseLR_maskR,
                    phaseRL_maskL = phaseRL_maskL,
                    phaseRL_maskR = phaseRL_maskR)
save(SPM_results, file = "HCP_Social_Task/Analysis/Sensitivity_Sim/SPMpar_MCMCdgm/SPM_results.RData")


