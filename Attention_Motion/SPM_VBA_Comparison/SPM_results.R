# SPM attention to motion results

library(R.matlab)

############ Functions ##############
# Unstructure params function
unstruct_paramMats = function(params, idxs){
  nu_A = params$A[idxs$A_idxs]
  nu_C = params$C[idxs$C_idxs]
  nu_B = c()
  for(i in 1:nrow(idxs$B_idxs)){
    nu_B[i] = params$B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]]
  }
  c(nu_A = nu_A, nu_B = nu_B, nu_C = nu_C)
}

# Function for summary
get_summary <- function(mu,sigma2){
  
  # Set alpha level
  alpha = 0.1
  
  intervals <- matrix(NA, nrow=length(mu), ncol=2)
  for (i in 1:nrow(intervals)){
    q_lower <- qnorm(alpha/2, mean = mu[i], sd = sqrt(sigma2[i]))
    q_upper <- qnorm(1-alpha/2, mean = mu[i], sd = sqrt(sigma2[i]))
    intervals[i,] <- c(q_lower, q_upper)
  }
  
  df <- data.frame(mu,intervals)
  colnames(df) <- c("mean","q5","q95")
  
  df$CI_length <- abs(df$q95-df$q5)
  
  return(df)
}
#####################################

# Indices of parameters
A_idxs <- matrix(c(1,2,
                   2,1,
                   2,3,
                   3,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,2,1,
                   3,2,1), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

# Read in file
means <- readMat("Attention_Motion/SPM VBA Comparison/means_SPM.mat")
vars <- readMat("Attention_Motion/SPM VBA Comparison/variances_SPM.mat")

# Grab posterior means
means$B <- list(means$B[,,1],means$B[,,2],means$B[,,3])

# Grab posterior variances
vars$B <- list(vars$B[,,1],vars$B[,,2],vars$B[,,3])

# Unstructure means and variances
mu <- unstruct_paramMats(means,idxs)
sigma2 <- unstruct_paramMats(vars,idxs)

# Get results
result <- get_summary(mu,sigma2)
result <- round(result, digits=4)

# Export results
write.csv(result, "Attention_Motion/SPM VBA Comparison/SPM_summary.csv")



