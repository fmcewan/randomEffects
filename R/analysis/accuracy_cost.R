# Part B1, B2, B3

# B1



# B2



# B3: Functions and settings for computational cost comparisons

# Default settings (to be able to pass this as one variable later)
b3_settings <- list(  
  cluster_counts = c(50, 200, 800),
  m = 5,                          # Patients per practice
  beta = c(-0.5, 1), sigma = 1.5, # Theta
  K_gh = 10, K_agh = 7,           # Nodes
  R_is = 100, df_is = 4,          # For importance sampling
  min_seconds = 0.2,              # For stable timings
  rounds = 3                      # Num rounds to be timed 
)

# Timed evaluations with system.time()
b3_time_evaluation <- function(f, settings=b3_settings) {
  min_seconds <- settings$min_seconds
  rounds <- settings$rounds
  
  invisible(f())    # Checking the call once without returning or timing
  repetitions <- 1  # Starting with 1 call
  # Defining helper to measure elapsed time
  measure_elapsed <- function(n) {
    system.time({ for (i in seq_len(n)) { invisible(f()) }})[["elapsed"]]
  }
  
  # Repeating until proper batch size is found, then timing whole batch
  repeat {
    elapsed <- measure_elapsed(repetitions)
    if (elapsed >= min_seconds) { break } # break if long enough
    repetitions <- repetitions * 2 # Otherwise increase num of calls and retry
  }
  
  # And repeating batches to show the variability in the timings
  seconds_per_call <- numeric(rounds)
  for (r in seq_len(rounds)) {
    elapsed <- measure_elapsed(repetitions)
    seconds_per_call[r] <- elapsed / repetitions # Finding time per repetition
  }
  
  # Individual batch timings and summaries
  list(
    median_seconds = median(seconds_per_call),
    min_seconds = min(seconds_per_call),
    max_seconds = max(seconds_per_call),
    calls_per_batch = repetitions,
    seconds_per_call = seconds_per_call
  )
}