#Description
#Super spreading model- Monte Carlo
library(MASS)

#Setup
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")
#par(mar=c(1,1,1,1))

################################################################################
# MCMC - ONE PARAMETERS AT A TIME
################################################################################

#*MCMC - Investigate one parameter at a time
mcmc_ss_r0_total_x1 <- function(data, n, sigma, alphaX, betaX, gammaX, param_infer, true_R0, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Initialise parameters
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_total_vec <- vector('numeric', n)
  r0_total_vec[1] <- x0
  count_accept = 0
  
  #alpha
  if (param_infer == 'alpha') {
    alpha_vec[1] <- x0
    beta_vec[1] <- betaX
    gamma_vec[1] <- gammaX
  }
  
  #beta
  if (param_infer == 'beta') {
    alpha_vec[1] <- alphaX
    beta_vec[1] <- x0
    gamma_vec[1] <- gammaX
  }
  
  #gamma
  if (param_infer == 'gamma') {
    alpha_vec[1] <- alphaX
    beta_vec[1] <- betaX
    gamma_vec[1] <- x0
  }
  
  # #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #alpha
    if (param_infer == 'alpha') {
      #alpha_vec[1] <- alphaX; beta_vec[1] <- betaX; gamma_vec[1] <- x0
      
      alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma)
      #cat("alpha dash: ", alpha_dash, "\n")
      
      if(alpha_dash < 0){
        alpha_dash = abs(alpha_dash)
      }
      
      #log alpha
      logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
      logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
      #prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
      #prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
      log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
      
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        alpha_vec[i] <- alpha_dash
        count_accept = count_accept + 1
      } else {
        alpha_vec[i] <- alpha_vec[i-1]
      }
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }

    #************************************************************************
    #beta
    if (param_infer == 'beta') {
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma)
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }

    logl_new = log_like_ss_lse(data, alpha_vec[i-1], beta_dash, gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    #prior1 = dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    #prior2 = dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2

    if (!(is.na(log_accept_prob)) &&
        log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept = count_accept + 1
    } else {
      beta_vec[i] <- beta_vec[i - 1]
    }
    } else {
      beta_vec[i] <- beta_vec[i - 1]
    }
    
    #************************************************************************
    #gamma
    if (param_infer == 'gamma') {
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma)
    
    if(gamma_dash < 0){
      gamma_dash = abs(gamma_dash)
    }
    
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash)
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept = count_accept + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    } else {
      gamma_vec[i] <- gamma_vec[i - 1]
    }
    
    #Calculate r0 value
    r0_total_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
  }
  #Final stats - Acceptance Rate
  accept_rate = 100*count_accept/n
  data_10_pc = 0.1*n
  
  #Mean of parameter infered
  if (param_infer == 'alpha') {
    param_mean_val = mean(alpha_vec[n-data_10_pc:n])
    param_mcmc = alpha_vec
  } else if (param_infer == 'beta') {
    param_mean_val = mean(beta_vec[n-data_10_pc:n])
    param_mcmc = beta_vec
  } else if (param_infer == 'gamma') {
    param_mean_val = mean(gamma_vec[n-data_10_pc:n])
    param_mcmc = gamma_vec
  } else {
    print('Set the parameter inference flag')
  }
  
  #Results - Create dataframe
  df_results <- data.frame(
    alpha = alphaX,
    beta = betaX,
    gamma = gammaX,
    param_infer = param_infer,
    param_mean = param_mean_val,
    R0 = true_R0, 
    R0_mean_MCMC = mean(r0_total_vec[n-data_10_pc:n]),
    accept_rate = accept_rate) 
  
  print(df_results)
  
  #Return alpha, acceptance rate
  return(list(param_mcmc, r0_total_vec, accept_rate))
}

#*****************************
#Plot results
plot_mcmc_results <- function(sim_data, param_mcmc, r0_total_vec, param_infer, dist_type, true_param, true_r0){
  
  #Set up
  par(mfrow=c(2,3))

  #Sim data
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste("Daily Infections, Super Spreading Events ", dist_type, ", r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Paramter MCMC chain(alpha, beta, gamma)
  plot.ts(param_mcmc, xlab = 'Time', ylab = param_infer,
          main = paste("MCMC chain of ", param_infer),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Paramter mean
  dat_10_pc = 0.1*n
  param_mean_val = mean(param_mcmc[n-dat_10_pc:n])
  param_mean = cumsum(param_mcmc)/seq_along(param_mcmc)
  plot(seq_along(param_mean), param_mean, xlab = 'Time', ylab = param_infer,
       main = paste("Mean of ", param_infer, ". Final Mean = ", round(param_mean_val, 2), ", True value = ", true_param),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #R0
  plot.ts(r0_total_mcmc,  xlab = 'Time', ylab = 'R0 total',
          main = paste("R0 total, MCMC Super spreading model, true total R0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #R0 mean
  r0_tot_mean = cumsum(r0_total_mcmc)/seq_along(r0_total_mcmc)
  r0_mean_val = mean(r0_total_mcmc[n-dat_10_pc:n])
  plot2 = plot(seq_along(r0_tot_mean), r0_tot_mean,
               xlab = 'Time', ylab = 'r0 total',
               main = paste("Mean of R0 total MCMC chain, Final Mean = ", round(r0_mean_val, 2), ", True R0 = ", true_r0),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  
  #Histogram
  hist(r0_total_mcmc, freq = FALSE, xlab = 'R0 total', #ylab = 'Density', 
       main = 'R0 total MCMC samples',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) #Empirical density of 
  
}

#********
#Implement
num_days = 50

#Gamma params for infectiousness curve (lambda) distribution
shape_gamma = 6
scale_gamma = 1 

############# INSERT PARAMETERS #!!---#######################################
alphaX = 1.1 # 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.05 #0.2 #0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10
param_infer = 'gamma'   #!!!!!!!!!!!
true_param_val = gammaX #!!!!!!!!!!!
true_r0 = alphaX + betaX*gammaX
true_r0
##!!---##############################################################!!---##

#Epidemic data
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#sim_data = sim_data2
#sim_data = sim_data3

#Time 
n = 10000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params = mcmc_ss_r0_total_x1(sim_data, n, sigma, alphaX, betaX, gammaX, param_infer, true_r0)
end_time = Sys.time()
time_elap = end_time - start_time
print('Time elapsed:')
print(time_elap)

#Extract params
param_mcmc = mcmc_params[1]
param_mcmc = unlist(param_mcmc)

r0_total_mcmc = mcmc_params[2]
r0_total_mcmc = unlist(r0_total_mcmc)

accept_rate = mcmc_params[3]

#Plotting
#Apply
dist_type = 'Neg Bin dist'
plot_mcmc_results(sim_data, param_mcmc, r0_total_mcmc, param_infer, dist_type, true_param_val, true_r0)


#********
#Implement
num_days = 50

#Gamma params for infectiousness curve (lambda) distribution
shape_gamma = 6
scale_gamma = 1 

############# INSERT PARAMETERS #!!---#######################################
alphaX = 1.1 # 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.05 #0.2 #0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10
param_infer = 'gamma'   #!!!!!!!!!!!
true_param_val = gammaX #!!!!!!!!!!!
true_r0 = alphaX + betaX*gammaX
true_r0
##!!---##############################################################!!---##

#Epidemic data
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#sim_data = sim_data2
#sim_data = sim_data3

#Time 
n = 10000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params = mcmc_ss_r0_total_x2(sim_data, n, sigma, alphaX, betaX, gammaX, param_infer, true_r0)
end_time = Sys.time()
time_elap = end_time - start_time
print('Time elapsed:')
print(time_elap)

################################################################################
# MCMC - TWO PARAMETERS AT A TIME
################################################################################

#*MCMC - Investigate one parameter at a time
mcmc_ss_x2 <- function(data, n, sigma, sigma_b, alphaX, betaX, gammaX, param_infer, true_R0, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Initialise parameters
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  count_accept1 = 0; count_accept2 = 0; count_accept3 = 0;
  
  #alpha & gamma
  if (param_infer == 'alpha_gamma') {
    alpha_vec[1] <- x0
    beta_vec[1] <- betaX
    gamma_vec[1] <- x0
  }
  
  #beta
  if (param_infer == 'alpha_beta') {
    alpha_vec[1] <- x0
    beta_vec[1] <- x0
    gamma_vec[1] <- gammaX
  }
  
  #gamma
  if (param_infer == 'beta_gamma') {
    alpha_vec[1] <- alphaX
    beta_vec[1] <- x0
    gamma_vec[1] <- x0
  }
  
  # #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #alpha
    if ((param_infer == 'alpha_gamma') | (param_infer == 'alpha_beta')) {
      #alpha_vec[1] <- alphaX; beta_vec[1] <- betaX; gamma_vec[1] <- x0
      
      alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma)
      #cat("alpha dash: ", alpha_dash, "\n")
      
      if(alpha_dash < 0){
        alpha_dash = abs(alpha_dash)
      }
      
      #log alpha
      logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
      logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
      #prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
      #prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
      log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
      
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        alpha_vec[i] <- alpha_dash
        count_accept1 = count_accept1 + 1
      } else {
        alpha_vec[i] <- alpha_vec[i-1]
      }
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
    }
    
    #************************************************************************
    #beta
    if ((param_infer == 'alpha_beta') | (param_infer == 'beta_gamma')) {
      beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma_b)
      if(beta_dash < 0){
        beta_dash = abs(beta_dash)
      }
      
      logl_new = log_like_ss_lse(data, alpha_vec[i-1], beta_dash, gamma_vec[i-1])
      logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
      #prior1 = dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
      #prior2 = dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE)
      log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
      
      if (!(is.na(log_accept_prob)) &&
          log(runif(1)) < log_accept_prob) {
        beta_vec[i] <- beta_dash
        count_accept2 = count_accept2 + 1
      } else {
        beta_vec[i] <- beta_vec[i - 1]
      }
    } else {
      beta_vec[i] <- beta_vec[i - 1]
    }
    
    #************************************************************************
    #gamma
    if ((param_infer == 'alpha_gamma') | (param_infer == 'beta_gamma')) {
      gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma)
      
      if(gamma_dash < 0){
        gamma_dash = abs(gamma_dash)
      }
      
      #Acceptance Probability
      logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash)
      logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
      log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
      
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        gamma_vec[i] <- gamma_dash
        count_accept3 = count_accept3 + 1
      } else {
        gamma_vec[i] <- gamma_vec[i-1]
      }
    } else {
      gamma_vec[i] <- gamma_vec[i - 1]
    }
    
    #Calculate r0 value
    r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
  }
  
  #Final stats
  #alpha
  accept_rate1 = 100*count_accept1/n
  cat("Acceptance rate1 = ",accept_rate1)
  
  #beta
  accept_rate2 = 100*count_accept2/n
  cat("Acceptance rate2 = ", accept_rate2)
  
  #gamma
  accept_rate3 = 100*count_accept3/n
  cat("Acceptance rate3 = ", accept_rate3)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3))
}

#*****************************
#Plot results
plot_mcmc_results_total <- function(sim_data, mcmc_params, true_r0, dist_type, total_time){
  
  #Plot Set up
  plot.new()
  par(mfrow=c(2,4))
  
  #Extract params
  alpha_mcmc = mcmc_params[1]
  alpha_mcmc = unlist(alpha_mcmc)
  
  beta_mcmc = mcmc_params[2]
  beta_mcmc = unlist(beta_mcmc)
  
  gamma_mcmc = mcmc_params[3]
  gamma_mcmc = unlist(gamma_mcmc)
  
  r0_mcmc = mcmc_params[4]
  r0_mcmc = unlist(r0_mcmc)
  
  #Stats
  data_10_pc = 0.1*n
  a_mcmc_mean = round(mean(alpha_mcmc[n-data_10_pc:n]), 2)
  b_mcmc_mean = round(mean(beta_mcmc[n-data_10_pc:n]), 2)
  g_mcmc_mean = round(mean(gamma_mcmc[n-data_10_pc:n]), 2)
  r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
  
  #Plots
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste("Day Infts SS Evnts", dist_type, "r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC
  plot.ts(alpha_mcmc, ylab = 'alpha', main = paste("MCMC SS Events, true alpha = ", alphaX))
  plot.ts(beta_mcmc, ylab = 'beta', main = paste("MCMC SS Events, true beta = ", betaX))
  plot.ts(gamma_mcmc,  ylab = 'gamma', main = paste("MCMC SS Events, true gamma = ", gammaX))
  #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
  
  #Mean data
  #r0 Mean
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  plot2 = plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0))
  print(plot2)
  
  #alpha mean
  alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
  plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Alpha MCMC mean, True alpha = ",alphaX))
  print(plot2)
  
  #beta mean
  beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
  plot2 = plot(seq_along(beta_mean), beta_mean, xlab = 'Time', ylab = 'beta', main = paste("Beta MCMC mean, True beta = ",betaX))
  print(plot2)
  
  #gamma Mean
  gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
  plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Gamma MCMC mean, True gamma = ",gammaX))
  print(plot2)
  
  #Results
  df_results <- data.frame(
    alpha = alphaX,
    a_mc = a_mcmc_mean,
    beta = betaX,
    b_mc = b_mcmc_mean,
    gamma = gammaX,
    g_mc = g_mcmc_mean,
    R0 = true_r0, 
    R0_mc = r0_mcmc_mean,
    accept_rate_a = round(mcmc_params[[5]],2),
    a_rte_b = round(mcmc_params[[6]], 2),
    a_rte_g = round(mcmc_params[[7]],2),
    tot_time = total_time) 
  
  print(df_results)
  
}


############# --- INSERT PARAMETERS! --- ######################################
alphaX = 0.8 # 1.1 #0.8 #1.1 # 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.2 #0.05 #0.2 #0.05 #0.2 #0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10
true_r0 = alphaX + betaX*gammaX
true_r0
##!!---##############################################################!!---##

#Epidemic data - Neg Bin
dist_type = 'Neg bin'
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#Epidemic data - Poisson 
#sim_data = simulate_ss_poisson(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
#plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#MCMC 
param_infer = 'alpha_gamma'
n = 10000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
sigma_b = 0.05
mcmc_params_x2 = mcmc_ss_x2(sim_data, n, sigma, sigma_b, alphaX, betaX, gammaX, param_infer, true_R0, x0 = 1)
#mcmc_params = mcmc_super_spreading(sim_data, n, sigma, sigma_b, x0 = 1)
end_time = Sys.time()
time_elap = round(end_time - start_time, 2)
print('Time elapsed:')
print(time_elap)


#Plotting
dist_type = 'Neg Bin,'
plot_mcmc_results_total(sim_data, mcmc_params_x2, true_r0, dist_type, time_elap)


################################################################################
# ALL THREE PARAMETERS AT ONCE
################################################################################

#********************************************************************
#MCMC Super-spreading
mcmc_super_spreading <- function(data, n, sigma,  sigma_b, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  
  count_accept1 = 0; count_accept2 = 0; count_accept3 = 0;
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #alpha
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma) 
    #cat("alpha dash: ", alpha_dash, "\n")
    
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
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
    prior1 = dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2 
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
    }
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma) 
    
    if(gamma_dash < 0){
      gamma_dash = abs(gamma_dash)
    }
    
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    #prior1 = dgamma(gamma_dash, shape = 1, scale = 1, log = TRUE)
    #prior2 = dgamma(gamma_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2 
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    
    r0_vec[i] = alpha_vec[i] + beta_vec[i]*gamma_vec[i]
    
  }
  #Final stats
  #alpha
  accept_rate1 = 100*count_accept1/n
  cat("Acceptance rate1 = ",accept_rate1)
  
  #beta
  accept_rate2 = 100*count_accept2/n
  cat("Acceptance rate2 = ", accept_rate2)
  
  #gamma
  accept_rate3 = 100*count_accept3/n
  cat("Acceptance rate3 = ", accept_rate3)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3))
}


############# --- INSERT PARAMETERS! --- ######################################
alphaX = 1.1 #0.8 #1.1 # 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.2 #0.05 #0.2 #0.05 #0.2 #0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10
true_r0 = alphaX + betaX*gammaX
true_r0
##!!---##############################################################!!---##

#Epidemic data - Neg Bin
#sim_data2 = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
#plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#Epidemic data - Poisson 
sim_data = simulate_ss_poisson(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#MCMC 
n = 30000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
sigma_b = 0.05
mcmc_params = mcmc_super_spreading(sim_data, n, sigma, sigma_b, x0 = 1)
end_time = Sys.time()
time_elap = round(end_time - start_time, 2)
print('Time elapsed:')
print(time_elap)

#Extract params
# r0_mcmc = mcmc_params[4]
# r0_mcmc = unlist(r0_mcmc)

#Plotting
dist_type = 'Poisson,'
plot_mcmc_results_total(sim_data, mcmc_params, true_r0, dist_type, time_elap)
