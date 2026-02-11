# Checking single subject convergence for SPM

# Directory where log .out files are stored
log_dir <- paste(getwd(),"HCP_Social_Task/SPM/Output",sep = "/")

# Get all .out files in that directory
files <- list.files(log_dir, pattern = "\\.out$", full.names = TRUE)

# Function to count "convergence" occurrences in one file
count_convergence <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  sum(grepl("convergence", lines, ignore.case = TRUE))
}

# Apply function to all files
conv_counts <- sapply(files, count_convergence)

# Make a tidy summary
summary_df <- data.frame(
  file = basename(files),
  convergence_count = conv_counts
)

summary_df <- summary_df[order(summary_df$convergence_count), ]
table(summary_df$convergence_count) # all have converged
