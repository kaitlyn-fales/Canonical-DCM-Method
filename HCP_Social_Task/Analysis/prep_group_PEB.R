# PEB for combining subjects into group DCM
library(R.matlab)

rm(list = ls())

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

# Hypothesis
A_idxs <- matrix(c(2,1,
                   1,2,
                   1,1,
                   2,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,2,1,
                   2,1,2,
                   2,1,1,
                   2,2,2), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1,
                   2,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

# Number of regions and inputs
n_regions <- 2
n_inputs  <- 2

# Initialize empty matrices
A_mat <- matrix(0, n_regions, n_regions)
B_mat <- array(0, dim = c(n_inputs, n_regions, n_regions))
C_mat <- matrix(0, n_regions, n_inputs)

# Fill matrices based on your indices
for (r in 1:nrow(A_idxs)) {
  A_mat[A_idxs[r,1], A_idxs[r,2]] <- 1
}

for (r in 1:nrow(B_idxs)) {
  B_mat[B_idxs[r,3], B_idxs[r,2], B_idxs[r,1]] <- 1
}

for (r in 1:nrow(C_idxs)) {
  C_mat[C_idxs[r,1], C_idxs[r,2]] <- 1
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
  
  return(list(A=A,B=B,C=C))
  
} # diagonal not modified for SPM - need raw values

# Extract posterior means and covariances

###############################################

for (i in 1:length(subjects)){
  
  for (j in 1:nrow(combos)){
    # Load your R posterior draws
    load(paste0("HCP_Social_Task/Output/sub-", subjects[i],
                "_phase", combos$phase[j], "_mask", combos$mask[j], "_draws.RData"))
    
    # Summarise posterior draws
    summary <- posterior::summarise_draws(draws_df)[4:13, ]  
    
    # Do posterior means
    nu = list(nu_A = summary$mean[1:4], 
              nu_B = summary$mean[5:8], 
              nu_C = summary$mean[9:10])
    
    paramMats <- struct_paramMats(m = n_regions, n_u = n_inputs, idxs = idxs, nu = nu)
    
    # Do posterior covariance
    draws_mat <- as_draws_matrix(draws_df)[,4:13]
    
    # Change order of columns of draws matrix to match the order of SPM
    SPM_indices <- c(3,1,2,4,7,5,6,8,9,10)
    draws_mat <- draws_mat[,SPM_indices]
    
    # Add zero columns to pad matrix for the params not est in hypothesis to match SPM
    zero_pad4 <- matrix(0, ncol = 4, nrow = nrow(draws_mat))
    zero_pad2 <- matrix(0, ncol = 2, nrow = nrow(draws_mat))
    
    # Adding 4 zero col for B[1] since there is no B param for input1
    # Adding 2 zero col for final two param in C since only two estimated
    draws_mat <- cbind(draws_mat[,1:4],zero_pad4,draws_mat[,5:10],zero_pad2)
    
    Cp <- round(cov(draws_mat), digits = 8)
    rownames(Cp) <- NULL
    colnames(Cp) <- NULL
    
    # Save as .mat file
    writeMat(paste0("HCP_Social_Task/Output/sub_",  subjects[i], "_phase", combos$phase[j], 
                    "_mask", combos$mask[j], ".mat"),
             A = paramMats$A,
             B1 = paramMats$B[[1]],
             B2 = paramMats$B[[2]],
             C = paramMats$C,
             Cp = Cp)
  }
}

