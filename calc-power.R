library(bindata)
library(geeM)
library(dplyr)
options(warn=-1)

# Note: To enable replicability of results of the Monte Carlo simulations,
# we set the seed using the function clusterSetRNGStream()
# This function was used in tandem with other functions in the package
# 'parallel' which was used towards parallelizing computations

# -----------------------------------------------------------------------------
# Uncomment lines within this section if not using run-all-calculations.R
# -----------------------------------------------------------------------------

# Probability that randomization at a given decision point will occur
# For example, if randomization will not occur if the participant is
# either driving, or has activated "sleep mode" in their phone,
# prob_avail is the probability that the participant is NOT 
# driving or in "sleep mode" at a decision point
#prob_avail <- 0.95 

# Elicited parameters from domain scientist: 
# probability of engaging in self-regulatory activities
#prob_engaged_given_none <- .15
#prob_engaged_given_low_effort <- .21
#prob_engaged_given_high_effort <- .16 

# Parameter which governs the simulated within-person correlation;
# The actual range of values of the simulated within-person correlation
# may differ from the value of rho.
#rho  <- 0.90

# -----------------------------------------------------------------------------
# Determine how many Monte Carlo samples to generate
# and number of observation within each Monte Carlo sample
# -----------------------------------------------------------------------------
N_simulations <- 1500  # No. of Monte Carlo simulations
N_participants <- 100  # No. of participants
D <- 60 # Total no. of decision points per participant over entire duration of the study

# -----------------------------------------------------------------------------
# Coin flip probabilities
# -----------------------------------------------------------------------------
prob_any <- .50
prob_none <- .50
prob_low_effort <- .25
prob_high_effort <- .25
prob_low_effort_given_any <- prob_low_effort/prob_any
prob_high_effort_given_any <- prob_high_effort/prob_any

# -----------------------------------------------------------------------------
# Calculate true value of regression model parameters
# -----------------------------------------------------------------------------
theta0 <- log(prob_engaged_given_none)
theta1 <- log(prob_engaged_given_low_effort/prob_engaged_given_none)
theta2 <- log(prob_engaged_given_high_effort/prob_engaged_given_none)

beta0 <- theta0
beta1 <- log(prob_low_effort_given_any*exp(theta1) + prob_high_effort_given_any*exp(theta2))

theta_vec <- as.matrix(c(theta0, theta1, theta2))
beta_vec <- as.matrix(c(beta0, beta1))

print(theta_vec)
print(beta_vec)

# -----------------------------------------------------------------------------
# Calculate true value of Primary Aim & Seconday Aim effects
# -----------------------------------------------------------------------------
effect_aim1 <- exp(beta1)
effect_aim2 <- exp(theta1 - theta2)

print(effect_aim1)
print(effect_aim2)

# -----------------------------------------------------------------------------
# Define functions to be used to calculate power via Monte Carlo simulation
# -----------------------------------------------------------------------------

GenerateData <- function(sim_number=1, N_participants=100, D=60, prob_none=0.50, prob_A1=0.25, prob_A2=0.25, prob_avail=1, theta_vec, rho=0.60, avail_stochastic=FALSE){
  
  # True value of data generating parameters
  truth_theta0 <- theta_vec[1,1]
  truth_theta1 <- theta_vec[2,1]
  truth_theta2 <- theta_vec[3,1]
  
  # Generate data for each of the N_participants
  list_dat <- list()
  
  for(i in 1:N_participants){
    # Randomization assignment
    mat_intervention_now <- rmultinom(n = D, size = 1, prob = c(prob_none, prob_A1, prob_A2))
    A1_now <- mat_intervention_now[2,]
    A2_now <- mat_intervention_now[3,]
    
    if(avail_stochastic==TRUE){
      avail_now <- rbinom(n = D, size = 1, prob = prob_avail)
    }else{
      which_avail <- sample.int(D, size = ceiling(prob_avail*(D)))
      timepoints <- 1:(D)
      avail_now <- ifelse(timepoints %in% which_avail, 1, 0)
    }
    
    # Construct data frame for participant
    this_participant_data <- data.frame(sim_number = sim_number,
                                        participant_id = i, 
                                        ema_id = 1:(D),
                                        t = 1:(D), 
                                        A1_now = A1_now, 
                                        A2_now = A2_now, 
                                        A_any_now = 1*(A1_now + A2_now>0),
                                        I_now = avail_now, 
                                        mu_next = NA)
    
    # True mean
    truth_prob_engaged_next <- exp(truth_theta0 + truth_theta1*this_participant_data[["A1_now"]] + truth_theta2*this_participant_data[["A2_now"]])
    this_participant_data[["mu_next"]] <- truth_prob_engaged_next
    this_participant_data[["mu_next"]] <- replace(this_participant_data[["mu_next"]], this_participant_data[["I_now"]]==0, NA)
    
    # Identify available decision points
    idx_avail <- which(this_participant_data[["I_now"]]==1)
    mu_next_given_avail <- this_participant_data[["mu_next"]][idx_avail]
    
    # Correlation matrix
    corr_dim <- sum(avail_now)
    corrmat <- diag(corr_dim) + rho*((as.matrix(rep(1,corr_dim))) %*% t(as.matrix(rep(1,corr_dim))) - diag(corr_dim))
    
    # Independent binary outcome
    #Y_next_given_avail <- rbinom(n = rep(1, length(idx_avail)), size = rep(1, length(idx_avail)), prob = mu_next_given_avail)
    
    # Correlated binary outcome
    Y_next_given_avail <- t(rmvbin(1, margprob = mu_next_given_avail, sigma = corrmat))
    
    ema_id_avail <- this_participant_data[["ema_id"]][idx_avail]
    dat_avail <- data.frame(ema_id = ema_id_avail, Y_next = Y_next_given_avail)
    this_participant_data <- left_join(x = this_participant_data, y = dat_avail, by = "ema_id")
    
    # This is a sanity check
    #check <- this_participant_data %>% group_by(mu_next) %>% summarise(mean(Y_next))
    
    # Add current participant data to list_dat
    list_dat <- append(list_dat, list(this_participant_data))
  }
  
  
  all_participant_dat <- do.call(rbind, list_dat)
  
  return(all_participant_dat)
}

FitPrimaryModel <- function(simdat){
  
  curr_sim_number <- simdat$sim_number[1]
  
  fit <- try(geem(Y_next ~ 1 + A_any_now, data = simdat, family = poisson, id = participant_id, waves = ema_id, corstr = "independence"), silent = TRUE)
  
  if(class(fit) == "try-error"){
    all_output <- list(sim_number = curr_sim_number, converged = FALSE, est_beta = NA, vcov_matrix_beta = NA)
  }else{
    est_beta <- as.matrix(fit$beta)
    vcov_matrix_beta <- fit$var
    converged <- fit$converged
    all_output <- list(sim_number = curr_sim_number, converged = converged, est_beta = est_beta, vcov_matrix_beta = vcov_matrix_beta)
  }
  
  return(all_output)
}

FitSecondaryModel <- function(simdat){
  
  curr_sim_number <- simdat$sim_number[1]
  
  fit <- try(geem(Y_next ~ 1 + A1_now + A2_now, data = simdat, family = poisson, id = participant_id, waves = ema_id), silent = TRUE)
  
  if(class(fit) == "try-error"){
    all_output <- list(sim_number = curr_sim_number, converged = FALSE, est_beta = NA, vcov_matrix_beta = NA)
  }else{
    est_beta <- as.matrix(fit$beta)
    vcov_matrix_beta <- fit$var
    converged <- fit$converged
    all_output <- list(sim_number = curr_sim_number, converged = converged, est_beta = est_beta, vcov_matrix_beta = vcov_matrix_beta)
  }
  
  return(all_output)
}

TestHypothesis <- function(curr_list, L){
  
  if(curr_list$converged == TRUE){
    use_alpha <- 0.05
    linear_combo <- L %*% curr_list$est_beta
    S_mat <- L %*% curr_list$vcov_matrix_beta %*% t(L)
    Z_statistic <- linear_combo/sqrt(S_mat)
    is_reject <- (abs(Z_statistic) > qnorm(p = 1 - use_alpha/2))
    is_reject <- (1*is_reject)
  }else{
    is_reject <- NA
  }
  return(is_reject)
}

# -----------------------------------------------------------------------------
# Test out functions for one simulation
# -----------------------------------------------------------------------------

if(FALSE){
  simdat <- GenerateData(sim_number = 1, 
                         N_participants = N_participants, 
                         D = D, 
                         prob_none = prob_none, 
                         prob_A1 = prob_low_effort, 
                         prob_A2 = prob_high_effort, 
                         prob_avail = prob_avail, 
                         theta_vec = theta_vec, 
                         rho = rho,
                         avail_stochastic=FALSE)
  head(simdat)
  all_output_primary <- FitPrimaryModel(simdat = simdat)
  all_output_secondary <- FitSecondaryModel(simdat = simdat)
  print(all_output_primary)
  print(all_output_secondary)
}

# -----------------------------------------------------------------------------
# Monte Carlo simulation
# -----------------------------------------------------------------------------
library(parallel)

ncore <- detectCores()
cl <- makeCluster(ncore - 1)
clusterSetRNGStream(cl, 102399)

clusterEvalQ(cl,{
  library(bindata)
  library(geeM)
  library(dplyr)
})

clusterExport(cl, c("N_simulations", "N_participants", "D", "prob_none", "prob_low_effort", "prob_high_effort", "prob_avail", "theta_vec", "rho"))
list_simdat <- as.list(1:(N_simulations))

# Data Generation
list_simdat <- parLapply(cl, 
                         list_simdat, 
                         GenerateData, 
                         N_participants = N_participants, 
                         D = D, 
                         prob_none = prob_none, 
                         prob_A1 = prob_low_effort, 
                         prob_A2 = prob_high_effort, 
                         prob_avail = prob_avail, 
                         theta_vec = theta_vec, 
                         rho = rho,
                         avail_stochastic=FALSE)

# Estimation
list_all_output_primary <- parLapply(cl, list_simdat, FitPrimaryModel)
list_all_output_secondary <- parLapply(cl, list_simdat, FitSecondaryModel)
stopCluster(cl)

# -----------------------------------------------------------------------------
# Check: How many simulation runs converged?
# -----------------------------------------------------------------------------
get_conv_primary <- lapply(list_all_output_primary, function(x){x$converged})
get_conv_secondary <- lapply(list_all_output_secondary, function(x){x$converged})

sum(do.call(rbind, get_conv_primary))
sum(do.call(rbind, get_conv_secondary))

# -----------------------------------------------------------------------------
# Check: Calculate empirical mean of simulated data across all decision points
# for which when I_now==1
# -----------------------------------------------------------------------------
list_check <- lapply(list_simdat, function(x){
  x %>% filter(I_now==1) %>% group_by(mu_next) %>% summarise(mean(Y_next, na.rm=TRUE))
  })
list_check <- lapply(list_check, function(x){x[,2]})
dat_check <- do.call(cbind, list_check)
rowMeans(dat_check)

# -----------------------------------------------------------------------------
# Check: Calculate bias across all N_simulations
# -----------------------------------------------------------------------------
list_beta_vec_estimates <- lapply(list_all_output_primary, function(x){x$est_beta})
list_theta_vec_estimates <- lapply(list_all_output_secondary, function(x){x$est_beta})

list_beta_vec_estimates <- do.call(cbind, list_beta_vec_estimates)
list_theta_vec_estimates <- do.call(cbind, list_theta_vec_estimates)

mean_beta_vec_estimates <- rowMeans(list_beta_vec_estimates)
mean_theta_vec_estimates <- rowMeans(list_theta_vec_estimates)

dat_check_bias_beta_vec <- data.frame(truth = beta_vec, mean_of_estimates = mean_beta_vec_estimates)
dat_check_bias_theta_vec <- data.frame(truth = theta_vec, mean_of_estimates = mean_theta_vec_estimates)
dat_check_bias_beta_vec[['bias_beta_vec']] <- dat_check_bias_beta_vec[['truth']] - dat_check_bias_beta_vec[['mean_of_estimates']]
dat_check_bias_theta_vec[['bias_theta_vec']] <- dat_check_bias_theta_vec[['truth']] - dat_check_bias_theta_vec[['mean_of_estimates']]

print(dat_check_bias_beta_vec)
print(dat_check_bias_theta_vec)

# -----------------------------------------------------------------------------
# Check: Relationship between beta1 and theta1 and theta2
# -----------------------------------------------------------------------------
mean_beta_vec_estimates[1] - mean_theta_vec_estimates[1]
mean_beta_vec_estimates[2] - log(prob_low_effort_given_any*exp(mean_theta_vec_estimates[2]) + prob_low_effort_given_any*exp(mean_theta_vec_estimates[3]))

# -----------------------------------------------------------------------------
# Check: Simulated correlation
# -----------------------------------------------------------------------------
list_corrmat <- lapply(list_simdat, function(x){
  curr_simdat <- x
  curr_simdat <- curr_simdat %>% select(participant_id, ema_id, Y_next)
  wide_curr_simdat <- reshape(data = curr_simdat, direction = "wide", idvar = "participant_id", timevar = "ema_id")
  wide_curr_simdat <- wide_curr_simdat[,-1]
  curr_corrmat <- cor(wide_curr_simdat, use = "pairwise.complete.obs")
  return(curr_corrmat)
})

empirical_corrmat <- Reduce(`+`, list_corrmat)/N_simulations
lower_tri_corr <- empirical_corrmat[lower.tri(empirical_corrmat)]
min_simulated_corr <- min(lower_tri_corr, na.rm=TRUE)
max_simulated_corr <- max(lower_tri_corr, na.rm=TRUE)
print(min_simulated_corr)
print(max_simulated_corr)

# -----------------------------------------------------------------------------
# Test hypothesis: Primary
# -----------------------------------------------------------------------------
L_primary <- matrix(c(0,1), nrow = 1, byrow = TRUE)
list_test_result_primary <- lapply(list_all_output_primary, TestHypothesis, L = L_primary)
test_result_primary <- do.call(rbind, list_test_result_primary)
power_primary <- mean(as.array(test_result_primary), na.rm = TRUE)

# -----------------------------------------------------------------------------
# Test hypothesis: Secondary
# -----------------------------------------------------------------------------
L_secondary <- matrix(c(0,1,-1), nrow = 1, byrow = TRUE)
list_test_result_secondary <- lapply(list_all_output_secondary, TestHypothesis, L = L_secondary)
test_result_secondary <- do.call(rbind, list_test_result_secondary)
power_secondary <- mean(as.array(test_result_secondary), na.rm = TRUE)

# -----------------------------------------------------------------------------
# Collate all results
# -----------------------------------------------------------------------------
dat_all_results <- data.frame(min_simulated_corr,
                              max_simulated_corr,
                              prob_engaged_given_none = prob_engaged_given_none,
                              prob_engaged_given_low_effort = prob_engaged_given_low_effort,
                              prob_engaged_given_high_effort = prob_engaged_given_high_effort,
                              RR_primary = exp(beta_vec[2]),
                              RR_secondary = exp(theta_vec[2] - theta_vec[3]),
                              prob_avail = prob_avail,
                              power_primary = power_primary,
                              power_secondary = power_secondary)

print(dat_all_results)

