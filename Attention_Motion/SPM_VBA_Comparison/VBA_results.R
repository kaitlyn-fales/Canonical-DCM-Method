# VBA toolbox attention to motion results

library(R.matlab)

# Function
get_summary <- function(mu,Sigma){
  
  alpha = 0.05
  
  intervals <- matrix(NA, nrow=length(mu), ncol=2)
  for (i in 1:nrow(intervals)){
    q_lower <- qnorm(alpha/2, mean = mu[i], sd = sqrt(Sigma[i,i]))
    q_upper <- qnorm(1-alpha/2, mean = mu[i], sd = sqrt(Sigma[i,i]))
    intervals[i,] <- c(q_lower, q_upper)
  }
  
  df <- data.frame(mu,sd = sqrt(diag(Sigma)),intervals)
  colnames(df) <- c("mean","sd","2.5%","97.5%")
  
  return(df)
}


file <- readMat("Attention_Motion/SPM_VBA_Comparison/VBA_results.mat")
mu <- file$p2[[3]][1:7]
Sigma <- file$p2[[4]][1:7,1:7]
result <- get_summary(mu,Sigma)
result <- round(result, digits=4)

write.csv(result, "Attention_Motion/SPM_VBA_Comparison/VBA_summary.csv")



