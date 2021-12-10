#Aim: Model Criticism of Super-Spreading Events Model

#SETUP
library(MASS)
library(pracma)
library(tidyverse)
library(tibble)
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")

#Epidemic params
num_days = 50
shape_gamma = 6 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 1 

############# --- INSERT PARAMETERS! --- ######################################
alphaX = 0.8 #0.7 #0.8 #0.7 
betaX = 0.1 #0.05 #0.025 #0.2 #0.1 
gammaX = 10 #8
true_r0 = alphaX + betaX*gammaX
true_r0
model_params = c(alphaX, betaX, gammaX, true_r0)

#MCMC - sigma
sigma_a = 0.4*alphaX
sigma_b = 1.0*betaX 
sigma_g = 0.85*gammaX
sigma_bg = 1.5*gammaX
sigma = c(sigma_a, sigma_b, sigma_g, sigma_bg)


######################################################
#1. MCMC
mcmc_ss_x4 <- function(data, n, sigma, thinning_factor, folder_results, rep, burn_in, x0 = 1) {
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; count_accept2 = 0;
  count_accept3 = 0; count_accept4 = 0;
  prior = TRUE; flag_true = FALSE
  
  #Create folder for mcmc results 
  folder_mcmc = paste0(folder_results, '/mcmc')
  ifelse(!dir.exists(file.path(folder_mcmc)), dir.create(file.path(folder_mcmc), recursive = TRUE), FALSE)
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #ALPHA
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    
    #************************************************************************
    #BETA
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    #loglikelihood
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1]
    }
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
    }
    
    #************************************************************************
    #GAMMA
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma_g) 
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #abs(gamma_dash)
    }
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev 
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - gamma_dash + gamma_vec[i-1]
    }
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    
    #R0
    r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
    
    #*****************************************************
    #GAMMA-BETA
    gamma_dash <- gamma_vec[i] + rnorm(1, sd = sigma_bg)#Alter sigma_bg depending on acceptance rate. 
    #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
    if(gamma_dash < 1){ #If less then 1
      gamma_dash = 2 - gamma_dash #abs(gamma_dash)
    }
    #New beta 
    beta_new = (r0_vec[i] - alpha_vec[i])/gamma_dash #Proposing new Gamma AND Beta. Beta_dash = f(R0 & gamma_dash)
    
    if(beta_new >= 0){ #Only accept values of beta > 0
      
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_new, gamma_dash)
      logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i])
      log_accept_prob = logl_new - logl_prev  
      #Priors
      if (prior){
        log_accept_prob = log_accept_prob - beta_dash + beta_vec[i]
      }
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_new
        count_accept4 = count_accept4 + 1
      } 
    }
  }
  
  #Final stats
  accept_rate1 = 100*count_accept1/n
  accept_rate2 = 100*count_accept2/n
  accept_rate3 = 100*count_accept3/n
  accept_rate4 = 100*count_accept4/n
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4))
}

################
#1i. REPEAT MCMC
run_mcmc_reps <- function(n, n_reps, model_params, sigma, flag_dt, folder_results, burn_in){
  
  'Run mcmc for n_reps iterations and save'
  
  #Get model params
  alphaX = model_params[1]; betaX = model_params[2]
  gammaX = model_params[3]; r0 = model_params[4];
  
  #Data_type
  flag1 = flag_dt[1]; flag2 = flag_dt[2]; flag3 = flag_dt[3] 
  cat('r0 = ', r0, '\n'); 
  
  #Repeat for n reps
  for(rep in 1:n_reps) {
    
    cat('\n rep =', rep, '\n')
    #Folder
    folder_results_rep = paste0(folder_results, '/rep_', rep)
    cat('\n folder_results_rep =', folder_results_rep, '\n')
    ifelse(!dir.exists(file.path(folder_results_rep)), dir.create(file.path(folder_results_rep), recursive = TRUE), FALSE)
    
    #Simulate data
    if (flag1){
      sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      saveRDS(sim_data, file = paste0(folder_results_rep, '/sim_data.rds'))
    } else if (flag2){
      sim_data = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      cat('simulate ss individs')
      saveRDS(sim_data, file = paste0(folder_results_rep, '/sim_data.rds'))
    } else if (flag3) {
      sim_data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
      saveRDS(sim_data, file = paste0(folder_results_rep, '/sim_data.rds'))
      cat('simulate_branching')
    }
    
    #MCMC
    mcmc_params = mcmc_ss_x4(sim_data, n, sigma, thinning_factor, folder_results_rep, rep, burn_in)
    
    #SAVE MCMC PARAMS 
    saveRDS(mcmc_params, file = paste0(folder_results_rep, '/mcmc_params_rep_', rep, '.rds' ))
  
  }
  
}

######################################################
#2.  MODEL CRITICISM - GET SUMMARY STATS
get_mcmc_rep <- function(results_home, model_type, iter, rep){
  
  for(rep in 1:n_reps) {
    
    #Get results
    results_rep = paste0(results_home, model_type, "/iter_", iter, "/rep_", rep, '/')
    print(results_rep)
    mcmc_params <- readRDS(paste0(results_rep, '/mcmc_params_rep_', rep, '.rds' ))
    
    #Get parameters
    alpha_mcmc = mcmc_params[1]; alpha_mcmc = unlist(alpha_mcmc)
    beta_mcmc = mcmc_params[2]; beta_mcmc = unlist(beta_mcmc)
    gamma_mcmc = mcmc_params[3]; gamma_mcmc = unlist(gamma_mcmc)
    r0_mcmc = mcmc_params[4]; r0_mcmc = unlist(r0_mcmc)
    
    #Simulate data using thinned params
    for(i in seq(burn_in, n, by = thinning_factor)){
      
      #Simulate data
      sim_data_model_crit = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alpha_mcmc[i], beta_mcmc[i], gamma_mcmc[i])
      #Save data
      saveRDS(sim_data_model_crit, file = paste0(results_rep, 'mcmc/sim_data_iter_', i, '.rds' ))
      
      #Get summary stats. 
      if (i == burn_in) { #first rep
        #print('CREATE DF')
        flag_create = TRUE
        df_summary_stats = get_summary_stats_sim_dataX(data, rep, flag_create)
        flag_create = FALSE 
      } else {
        df_summary_stats[nrow(df_summary_stats) + 1, ] = get_summary_stats_sim_dataX(data, rep, flag_create)
      }
    }
    
    #Save rep summary stats
    #save df_summary_stats!!!!!!!
  }
}

#Get summary statisitcs
get_summary_stats <- function(data, rep, flag_create){
  
  'Get summary statisitcs of the simulated data'
  if (flag_create){
    
    #Df
    summary_stats_results = data.frame(
      sum_inf_counts = sum(sim_data_params),
      median_inf_count = median(sim_data_params),
      max_inf_count = max(sim_data_params),
      std_inf_counts = std(sim_data_params),
      val_75_infs_counts = quantile(sim_data_params)[4][1][1],
      val_87_5_infs_counts = mean(quantile(sim_data_params)[4][1][1], quantile(sim_data_params)[5][1][1]),
      max_dif = max(abs(diff(sim_data_params))),
      med_dif = median(abs(diff(sim_data_params))),
      mean_upper_dif = mean(c(quantile(abs(diff(sim_data_params)))[4][1][1], quantile(abs(diff(sim_data_params)))[5][1][1])),
      sum_1st_half  = sum(which(sim_data_params < quantile(sim_data_params)[3][1][1])),
      sum_2nd_half =  sum(which(sim_data_params > quantile(sim_data_params)[3][1][1]))
      
    )
    
  } else {
    #List
    summary_stats_results = list(sum(sim_data_params), median(sim_data_params), max(sim_data_params),
                                 std(sim_data_params), quantile(sim_data_params)[4][1][1], 
                                 mean(quantile(sim_data_params)[4][1][1], quantile(sim_data_params)[5][1][1]),
                                 max(abs(diff(sim_data_params))), median(abs(diff(sim_data_params))),
                                 mean(c(quantile(abs(diff(sim_data_params)))[4][1][1], quantile(abs(diff(sim_data_params)))[5][1][1])),
                                 sum_1st_half  = sum(which(sim_data_params < quantile(sim_data_params)[3][1][1])),
                                 sum_2nd_half =  sum(which(sim_data_params > quantile(sim_data_params)[3][1][1]))
    )
  }
  
  summary_stats_results
  
}

#All P VALS
get_summary_stats <- function(data, rep, flag_create){
  
  for(rep in 1:n_reps) {
    SUM_STATS =
    gET P VAL
  }
}

#APPLY MCMC
model_type = 'ss_events' #base_sim_sse_inf' #'ssi_sim_sse_inf'
flags_data_type = c(TRUE, FALSE, FALSE) #1)ss_events, 2) ss_individuals, 3) basline
iter = 1
base_folder_current = paste0('~/PhD_Warwick/Project_Epidemic_Modelling/Results/super_spreading_events/model_criticism_II/', '', model_type, '/iter_', iter)
print(base_folder_current)

#Repitions 
n = 5500
n_reps = 100
burn_in = 500
thinning_factor = 50 #(1/1000)*n;

run_mcmc_reps(n, n_reps, model_params, sigma, flags_data_type, base_folder_current, burn_in)
