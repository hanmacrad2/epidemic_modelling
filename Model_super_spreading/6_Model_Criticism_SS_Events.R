#Model Criticism
#Super spreading Events model- MCMC Inference + Model Criticism
#Setup
library(MASS)
library(pracma)
library(tidyverse)
library(tibble)
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")
source("plotting_functions.R")

#Epidemic params
num_days = 50
#Gamma params for infectiousness curve (lambda) distribution
shape_gamma = 6
scale_gamma = 1 
#seed_count = 1
#par(mar=c(1,1,1,1))

################################################################################
# MCMC - FOUR PARAMETER UPDATES
################################################################################

mcmc_ss_mod_crit <- function(data, n, sigma_a, sigma_b, sigma_g, sigma_bg, prior, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  cat('Prior = ', prior)
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  count_accept1 = 0; count_accept2 = 0; count_accept3 = 0; count_accept4 = 0;
  thinning_factor = (1/1000)*n
  #vec_data_simulated <- vector('numeric', n/thinning_factor)
  count_thin = 1
  #Set up summary stats df
  df_summary_stats <- data.frame(
   
   
  )
  
  vec_sum = vector('numeric', n)
  vec_median = vector('numeric', n)
  vec_mode = vector('numeric', n)
  vec_std = vector('numeric', n)
  vec_med_diff = vector('numeric', n)
  vec_n1 = vector('numeric', n)
  vec_n2 = vector('numeric', n)
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #alpha
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma_a) 
    #cat("alpha dash: ", alpha_dash, "\n")
    
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev  #+ prior1 - prior2
    
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    #if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
    if(log(runif(1)) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    
    #************************************************************************
    #beta
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b) 
    #cat("Beta dash: ", beta_dash, "\n")
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev
    
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - beta_dash + beta_vec[i-1]
    }
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
    }
    
    #************************************************************************
    #gamma
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
    if(log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    
    #R0
    r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
    
    #*****************************************************
    #Gamma-Beta
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
      if(log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_new
        count_accept4 = count_accept4 + 1
      } 
      
    }
    
    #Get Summary Statistics
    if(mod(i, thinning_factor) == 0){
      
      #If df of summary stats doesn't exist - create it
      if (!exists("df_summary_stats")) {
        
        flag_create = TRUE
        df_summary_stats = get_summary_stats(list_summary_stats_i, flag_create)
        flag_create = FALSE
        
      } else {

        df_summary_stats[nrow(df_summary_stats) + 1, ] = get_summary_stats(list_summary_stats_i, flag_create)
      }
      
      count_thin = count_thin + 1
    }
    
    
  }
  #Final stats
  #alpha
  accept_rate1 = 100*count_accept1/n
  cat("Acceptance rate1 = ",accept_rate1, '\n')
  
  #beta
  accept_rate2 = 100*count_accept2/n
  cat("Acceptance rate2 = ", accept_rate2, '\n')
  
  #gamma
  accept_rate3 = 100*count_accept3/n
  cat("Acceptance rate3 = ", accept_rate3, '\n')
  
  #Gamma-Beta
  accept_rate4 = 100*count_accept4/n
  cat("Acceptance rate4 = ", accept_rate4, '\n')
  
  #Create list of p-values
  list_p_vals = create_lst_p_vals(sim_data, df_summary_stats)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3,
              accept_rate4,
              list_p_vals))
}

#Get summary stats
get_summary_stats <- function(sim_data, alpha_vec_i, beta_vec_i, gamma_vec_i, create_df_flag){
  
  'Get summary statisitcs of the simulated data'
  #Summary stats
  sumX = sum(sim_data)
  medianX = median(sim_data)
  modeX = mode(sim_data)
  stdX = std(sim_data)
  med_dif = median(diff(check))
    
    if (create_df_flag){
      #Df
      summary_stats_results = data.frame(
        sumX = sum(sim_data),
        medianX = median(sim_data),
        modeX = mode(sim_data),
        stdX = std(sim_data),
        med_dif = median(diff(check)))
    } else {
      #List
      summary_stats_results = list(sum(sim_data), median(sim_data), mode(sim_data),
                               std(sim_data), median(diff(check)))
    }
  
  summary_stats_results
  
}

#Run for multiple
#Get p values
get_p_values <- function(n_reps, model_params){
  
  'Run model criticism for n_reps iterations to get a sample of p values for a number of
  different summary statistics'
  
  #Get model params
  alphaX = model_params[1]; betaX = model_params[2]
  gammaX = model_params[3]; r0 = model_params[4];
  cat('r0 = ', r0); 
  
  #Initialise vector of p values
  vec_sum = vector('numeric', n)
  vec_median = vector('numeric', n)
  vec_mode = vector('numeric', n)
  vec_std = vector('numeric', n)
  vec_med_diff = vector('numeric', n)
  vec_n1 = vector('numeric', n)
  vec_n2 = vector('numeric', n)
  
  #Repeat for n reps
  for(rep in 1:n_reps) {
    
    #Simulate data
    sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
      
    #MCMC
    mcmc_params = mcmc_ss_mod_crit(sim_data, n, sigma_a, sigma_b, sigma_g, sigma_bg, prior)
    list_p_vals = mcmc_params[9]
    list_p_vals = unlist(list_p_vals)
    
    if (!exists("df_p_vals")) {
      
      #Create df; sum etc
      df_p_vals = data.frame(
        sumX  = list_p_vals[1],
        medianX = list_p_vals[3],
        modeX = list_p_vals[4],
        stdX = list_p_vals[5],
        med_dif = list_p_vals[6],
      )
      
    } else {
      
      df_p_vals[nrow(df_p_vals) + 1, ] = list_p_vals
    }
    
  }
  
  #Plot df values
  #plot 3x3

}


#Model Criticism Function
plot_model_criticism <- function(mcmc_params, sim_data, max_sum_val) { 
  
  #Plot Model Criticism
  vec_mod_crit = mcmc_params[9]
  vec_mod_crit = unlist(vec_mod_crit)
  true_sum_inf = sum(sim_data)
  
  #P value
  lt = length(which(vec_mod_crit < true_sum_inf))
  print(lt)
  gt = length(which(vec_mod_crit > true_sum_inf))
  print(gt)
  min_val = min(lt, gt)
  pvalue = min_val/length(vec_mod_crit)
  
  #Check
  if (lt < gt){
    flag = 'lt (<)'
  } else if (gt < lt){
    flag = 'gt (>)'
  }
  
  #Histogram
  hist(vec_mod_crit[vec_mod_crit < max_sum_val], breaks = 100, #freq = FALSE, 
       #xlim = c(xmin, xmax),
       xlab = paste('Sum of Infecteds <', max_sum_val), ylab = 'Density',
       main = paste('Model criticism, true R0 = ', true_r0, '.',
                    'P value', flag, '=', pvalue),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_sum_inf, col = 'red', lwd = 2)
  
}


############# --- INSERT PARAMETERS! --- ######################################
alphaX = 0.8 #0.7 #0.8 #0.7 #0.8 #0.7 #0.7 #1.1 #0.8 #1.1 #0.8 #1.1 # 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.1 #0.05 #0.025 #0.2 #0.1 #0.2 #0.05 #0.1 #0.05 #0.2 #0.05 #0.2 #0.05 #0.2 #0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10 #8
true_r0 = alphaX + betaX*gammaX
true_r0
#Seed
#seed_count = 13
seed_count = seed_count + 1
seed_count
##---##############################################################---##
set.seed(seed_count)
#set.seed(9)

#Epidemic data - Neg Bin
sim_data2 = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = paste('Daily Infections count, true R0 = ', true_r0))
#sim_data2 = sim_data
#WITH PRIOR
#MCMC 
n = 30000
sigma_a = 0.4*alphaX
sigma_a
sigma_b = 1.0*betaX 
sigma_b
sigma_g = 0.85*gammaX
sigma_g
sigma_bg = 1.5*gammaX
sigma_bg
start_time = Sys.time()
print('Start time:')
print(start_time)
prior = TRUE

######################
#Get p values

#Plot x4x4
get_p_values(n_reps, model_params)

# mcmc_params = mcmc_ss_mod_crit(sim_data, n, sigma_a, sigma_b, sigma_g, sigma_bg, prior)
# end_time = Sys.time()
# time_elap = round(end_time - start_time, 2)
# print('Time elapsed:')
# print(time_elap)
# 
# #Apply
# dist_type = 'Neg Bin,'
# max_sum_val = 5000
# #plot_mcmc_x4_priors(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count, prior)
# plot_mcmc_x4_II(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count, prior, max_sum_val)
# 
# #Model Criticism
# par(mfrow = c(1,1))
# model_criticism(mcmc_params, sim_data, max_sum_val)
# 
# par(mfrow = c(2,1))