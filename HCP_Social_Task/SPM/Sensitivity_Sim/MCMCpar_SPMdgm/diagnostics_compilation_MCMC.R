library(dplyr)
library(posterior)
library(mcmcse)
library(momentLS)

path <- paste0(getwd(),"/HCP_Social_Task/SPM/Sensitivity_Sim/SPMpar_SPMdgm/Output")

# Vector of subjects to loop through
subjects <- c(100307,100408,101107,101309,101915,103111,103414,103818,105014,
              105115,106016,108828,110411,111312,111716,113619,113922,114419,
              115320,116524,117122,118528,118730,118932,120111,122317,122620,
              123117,123925,124422,125525,126325,127630,127933,128127,128632,
              129028,130013,130316,131217,131722,133019,133928,135225,135932,
              136833,138534,139637,140925,144832,146432,147737,148335,148840,
              149337,149539,149741,151223,151526,151627,153025,154734,156637,
              159340,160123,161731,162733,163129,176542,178950,188347,189450,
              190031,192540,196750,198451,199655,201111,208226,211417,211720,
              212318,214423,221319,239944,245333,280739,298051,366446,397760,
              414229,499566,654754,672756,751348,756055,792564,856766,857263,
              899885)

# Phase-encoding (run/session)
phase <- c("LR","RL")

# Type of mask (right or left)
mask_type <- c("L","R")

# Combinations of all phases and mask types
combos <- expand.grid(phase = phase,mask = mask_type)

# Empty list to store results
diagnostic_list <- list()

# For loop to make master list of diagnostics across all combinations and subjects
for (i in 1:length(subjects)){

    for (j in 1:nrow(combos)){
    
      # Load output diagnostic file
      load(paste0(path,"/sub-",subjects[i],"_phase",combos$phase[j],"_mask",combos$mask[j],"_diagnostics.RData"))
      
      # Put in large df
      all_diagnostics <- do.call(rbind, do.call(c, lapply(all_diagnostics, function(x) if (is.data.frame(x)) list(x) else x)))
      
      # Record number of iterations
      iterations <- nrow(all_diagnostics)/6
      
      # Collect max treedepths
      max_treedepth <- all_diagnostics %>% filter(Parameter == "treedepth__" & Value == 10) %>% nrow()/iterations
      
      # Collect divergences
      divergence <- all_diagnostics %>% filter(Parameter == "divergent__" & Value == 1) %>% nrow()/iterations
      
      # Collect E-BFMI
      energy <- all_diagnostics %>% filter(Parameter == "energy__")
      e_bfmi <- (sum(diff(energy$Value)^2) / (length(energy$Value)-2)) / var(energy$Value)
      
      # Collect average acceptance prob
      accept_stat <- as.numeric(all_diagnostics %>% filter(Parameter == "accept_stat__") %>% summarise(mean(Value)))
      
      # Collect stepsize
      stepsize <- as.numeric(all_diagnostics %>% filter(Parameter == "stepsize__") %>% head(n = 1) %>% select(Value))
      
      # Load in output draws file
      load(paste0(path,"/sub-",subjects[i],"_phase",combos$phase[j],"_mask",combos$mask[j],"_draws.RData"))
      
      # Get summary of param except lp_
      summary <- summarize_draws(draws_df, mcse_mean, ess_bulk, ess_tail)[-1,-1]
      avg_summary <- colMeans(summary)
      names(avg_summary) <- c("mean_mcse","mean_ess_bulk","mean_ess_tail")
      diag_cutoffs <- c(max_mcse = max(summary$mcse_mean),
                        min_ess_bulk = min(summary$ess_bulk),
                        min_ess_tail = min(summary$ess_tail))
      
      # Get multivariate ESS
      param_draws <- suppressWarnings(as.matrix(draws_df[,2:(ncol(draws_df)-3)]))
      avar <- momentLS::mtvMLSE(param_draws)$cov
      ess_multi <- multiESS(param_draws, covmat = avar)
      
      # Convergence indicators
      convergence <- c(mcse_ok = ifelse(diag_cutoffs[1]<0.01,T,F),
                       ess_bulk_ok = ifelse(diag_cutoffs[2]>100,T,F),
                       ess_tail_ok = ifelse(diag_cutoffs[3]>100,T,F),
                       ess_multi_good = ifelse(ess_multi>7859,T,F),
                       ess_multi_ok = ifelse(ess_multi>1965,T,F))
      
      # Combine into dataframe
      results <- as.data.frame.list(
        c(
          subject = subjects[i],
          phase = as.character(combos$phase[j]),
          mask = as.character(combos$mask[j]),
          iterations = iterations,
          prop_max_treedepth = max_treedepth,
          prop_divergence = divergence,
          e_bfmi = e_bfmi,
          mean_accept_stat = accept_stat,
          stepsize = stepsize,
          avg_summary,
          diag_cutoffs,
          ess_multi = ess_multi,
          convergence
        ),
        stringsAsFactors = FALSE
      )
      
      # Format dataframe properly
      names(results) <- sub("\\..*$", "", names(results))
      results <- type.convert(results, as.is = TRUE)
      
      assign(paste0("phase",combos$phase[j],"_mask",combos$mask[j]),results)
    }
  
  # Bind rows of combination results for subject i
  sub_results <- bind_rows(phaseLR_maskL,phaseLR_maskR,
                           phaseRL_maskL,phaseRL_maskR)

  # Add to result diagnostic list
  diagnostic_list[[i]] <- sub_results

}

diagnostics_df <- bind_rows(diagnostic_list)

save(diagnostics_df, file = paste0(path,"/../diagnostics_compilation.RData"))

# Summarize
summary(diagnostics_df)

