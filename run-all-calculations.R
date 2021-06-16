# -----------------------------------------------------------------------------
# Simulation Study 1
# -----------------------------------------------------------------------------

# Execute each scenario in Simulation Study 1
# Values for all parameters were held fixed in each calculation below, with
# the exception of the probability that a given decision point is available;
# the value of this parameter is varied from .80 to 1, in increments of .05

# Parameters fixed throughout all scenarios
prob_engaged_given_none <- .15
prob_engaged_given_low_effort <- .21
prob_engaged_given_high_effort <- .16 
rho  <- 0.90
# Parameters varied in the different scenarios
prob_avail <- 0.80 
source("calc-power.R")
write.csv(dat_all_results, 
          file.path("sim_study_01", paste("dat_all_results_",prob_avail,".csv", sep="")), 
          row.names = FALSE)
rm(list = ls())

# Parameters fixed throughout all scenarios
prob_engaged_given_none <- .15
prob_engaged_given_low_effort <- .21
prob_engaged_given_high_effort <- .16 
rho  <- 0.90
# Parameters varied in the different scenarios
prob_avail <- 0.85 
source("calc-power.R")
write.csv(dat_all_results, 
          file.path("sim_study_01", paste("dat_all_results_",prob_avail,".csv", sep="")), 
          row.names = FALSE)
rm(list = ls())

# Parameters fixed throughout all scenarios
prob_engaged_given_none <- .15
prob_engaged_given_low_effort <- .21
prob_engaged_given_high_effort <- .16 
rho  <- 0.90
# Parameters varied in the different scenarios
prob_avail <- 0.90
source("calc-power.R")
write.csv(dat_all_results, 
          file.path("sim_study_01", paste("dat_all_results_",prob_avail,".csv", sep="")), 
          row.names = FALSE)
rm(list = ls())

# Parameters fixed throughout all scenarios
prob_engaged_given_none <- .15
prob_engaged_given_low_effort <- .21
prob_engaged_given_high_effort <- .16 
rho  <- 0.90
# Parameters varied in the different scenarios
prob_avail <- 0.95 
source("calc-power.R")
write.csv(dat_all_results, 
          file.path("sim_study_01", paste("dat_all_results_",prob_avail,".csv", sep="")), 
          row.names = FALSE)
rm(list = ls())

# Parameters fixed throughout all scenarios
prob_engaged_given_none <- .15
prob_engaged_given_low_effort <- .21
prob_engaged_given_high_effort <- .16 
rho  <- 0.90
# Parameters varied in the different scenarios
prob_avail <- 1
source("calc-power.R")
write.csv(dat_all_results, 
          file.path("sim_study_01", paste("dat_all_results_",prob_avail,".csv", sep="")), 
          row.names = FALSE)
rm(list = ls())

# -----------------------------------------------------------------------------
# Simulation Study 2
# -----------------------------------------------------------------------------

prob_engaged_given_none <- .10
prob_engaged_given_low_effort <- .155
prob_engaged_given_high_effort <- .107
rho  <- 0.91
prob_avail <- 0.80 
source("calc-power.R")
write.csv(dat_all_results, 
          file.path("sim_study_02", paste("dat_all_results_",prob_avail,"_",prob_engaged_given_none, ".csv", sep="")), 
          row.names = FALSE)
rm(list = ls())

prob_engaged_given_none <- .20
prob_engaged_given_low_effort <- .267
prob_engaged_given_high_effort <- .21
rho  <- 0.89
prob_avail <- 0.80 
source("calc-power.R")
write.csv(dat_all_results, 
          file.path("sim_study_02", paste("dat_all_results_",prob_avail,"_",prob_engaged_given_none, ".csv", sep="")), 
          row.names = FALSE)
rm(list = ls())




