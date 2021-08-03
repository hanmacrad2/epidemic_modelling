#Description
#Adaptive Monte Carlo targeting r0
library(gridExtra)

#Setup
source("1_simulation.R")
source("functions.R")
par(mar=c(1,1,1,1))

#Params
num_days = 30 # 60 #100
r0 = 3.0
n = 50000
shape_gamma = 6
scale_gamma = 1
#Priors
prior_r0_k = 1
prior_r0_theta = 1
data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)


#***********************************
#Log Likelihood 
log_like <- function(y, r0_dash){
  
  #Params
  num_days = length(y)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda = r0_dash*sum(y[1:t-1]*rev(prob_infect[1:t-1]))
    logl = logl + y[t]*log(lambda) - lambda
    
  }
  
  logl
  
}

#Version 2
log_likeII <- function(x, r0_dash){
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda = r0_dash*sum(x[1:t-1]*rev(prob_infect[1:t-1]))
    logl = logl + x[t]*log(lambda) - lambda + log(1/(factorial(x[t])))
    
  }
  
  logl
  
}

#********************************************************************
#Adaptive MCMC
adaptive_mc_r0 <- function(data, n, sigma, x0 = 1, burn_in = 5000) { #burn_in = 2500
  
  'Returns mcmc samples of R0 & acceptance rate'
  
  #Set up
  r0_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  U <- runif(n)
  count_accept = 0
  count_reject = 0
  sd_sample = 1
  
  #MCMC chain
  for(i in 2:n) {
    
    #New Proposal
    Y <- r0_vec[i-1] + rnorm(1, sd = sigma)
    print(Y)
    if(Y < 0){
      Y = abs(Y)
    }
    
    log_alpha = log_likeII(data, Y) - log_likeII(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) #log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    print('log_alpha')
    print(log_alpha)
    #if (is.na(log_alpha)){
      #print('na value, Y value:')
      #print(Y)
    #}
    #if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
    if(log(U[i]) < log_alpha) {
      r0_vec[i] <- Y
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
      count_reject = count_reject + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma = var(r0_vec[2:i])*(2.38^2)
    }
    
  }
  #Final stats
  total_iters = count_accept + count_reject
  accept_rate = 100*(count_accept/(count_accept+count_reject))
  num_samples = count_accept
  print("Total iterations = ")
  print(total_iters)
  print("Acceptance rate = ")
  print(accept_rate)
  print("Number samples = ")
  print(count_accept)
  
  #Burn-in
  r0_vec = r0_vec[burn_in:n]
  
  #Return r0, acceptance rate
  return(list(r0_vec, accept_rate, num_samples, sigma))
}


#Apply
mcmc_params_ad = adaptive_mc_r0(data, n, sigma)
r0_as = mcmc_params_ad[1]
r0_as = unlist(r0_as)

#Plot
plot.ts(r0_as)


#Plots
mcmc_plotting_adaptive <- function(mcmc_vector, r0_true, folder_dir_ad) {
  
  #Folder save
  pdf(paste(folder_dir_ad, "/", "adpative_mc_r0_true_", r0_true, ".pdf", sep=""))
  
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector, ylab = 'R0', main = paste("Adaptive MC for R0, true R0 = ", r0_true))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  r0_mean = cumsum(mcmc_vector)/seq_along(mcmc_vector)
  plot2 = plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0', main = paste("Mean of R0 MCMC chain, True R0 = ",r0_true))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'r0', ylab = 'Density', 
               main = 'Empirical density of r0 - MCMC chain')
  print(hist3)
  dev.off()
  
}


#********************************************************************
#Apply Adaptive MCMC to a range of r0s

apply_adaptive_mc_range_r0 <- function(list_r0, sigma, folder_dir_ad){
  
  #Create folder
  ifelse(!dir.exists(file.path(folder_dir_ad)), dir.create(file.path(folder_dir_ad)), FALSE)
  list_accept_rate = vector('numeric', length(list_r0))
  list_sd = vector('numeric', length(list_r0))
  list_num_samp = vector('numeric', length(list_r0))
  list_time_taken = vector('numeric', length(list_r0))
  i = 1
  
  for (r0X in list_r0){
    
    #Get simulated data when r0 is r0X
    print(r0X)
    data = simulate_branching(num_days, r0X, shape_gamma, scale_gamma)
    
    #Time
    start_time = Sys.time()
    mcmc_params_ad = adaptive_mc_r0(data, n, sigma)
    end_time = Sys.time()
    time_elap = end_time - start_time
    print('Time elapsed:')
    print(time_elap)
    
    #Extract params
    r0_mcmc = mcmc_params_ad[1]
    r0_mcmc = unlist(r0_mcmc)
    
    accept_rate = mcmc_params_ad[[2]]
    list_accept_rate[i] = round(accept_rate, 2)
    
    num_samples = mcmc_params_ad[[3]]
    list_num_samp[i] = num_samples
    
    sd_final = mcmc_params_ad[[4]]
    list_sd[i] = round(sd_final, 3)
    
    list_time_taken[i] = round(time_elap, 2)
    i = i + 1
    
    #Apply Plotting
    mcmc_plotting_adaptive(r0_mcmc, r0X, folder_dir_ad)
    
  }
  
  #Create dataframe
  df_results <- data.frame(
    r0 = list_r0,
    accept_rate = unlist(list_accept_rate),
    n_samples = unlist(list_num_samp),
    time_sec = unlist(list_time_taken),
    sd_final = unlist(list_sd))
  
  print(df_results)
  
  df_results
}

#Apply
sigma = 0.75
folder_dir_ad = 'Results/Adaptive_MC/adaptive_mc_iter_VIII_burn_in_5k'
list_r0 = c(0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)  #c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0) #c(0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
#list_r0 = c(0.5, 0.65, 0.70, 0.75, 0.8, 0.85, 0.95, 1.05, 2.80, 3.05, 3.55, 4.05, 4.55, 5.05, 8.05, 10.05)
df_ad_results_formI = apply_adaptive_mc_range_r0(list_r0, sigma, folder_dir_ad)
df_ad_results_formI


#**********************************************************************
#*Adaptive Scaling Algorithm

adaptive_scaling_metropolis_r0 <- function(data, n, sigma, alpha_star, x0 = 1, burn_in = 5000) { #burn_in = 2500
  
  'Returns mcmc samples of R0 & acceptance rate'
  
  #Set up
  r0_vec <- vector('numeric', n)
  scaling_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  scaling_vec[1] <- x0
  U <- runif(n)
  count_accept = 0
  count_reject = 0
  #sd_sample = 1
  
  #MCMC chain
  for(i in 2:n) {
    
    #New Proposal
    Y <- r0_vec[i-1] + exp(scaling_vec[i-1])*rnorm(1) #sd = sigma) exp
    
    #if (is.na(Y) || is.nan(Y) || is.infinite(Y)) next
    
    if(Y < 0){
      Y = abs(Y)
    }
    print('y')
    print(Y)
    
    #Prints
    print('log-likelihoods')
    loglike1 = log_like(data, Y)
    print(loglike1)
    loglike2 = log_like(data, r0_vec[i-1])
    print(loglike2)
    dg1 = dgamma(Y, shape = 1, scale = 1, log = TRUE)
    print('priors')
    print(dg1)
    dg2 = dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE)
    print(dg2)
    
    log_alpha = log_likeII(data, Y) - log_likeII(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) #log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    print('log_alpha')
    print(log_alpha)
    
    #if (is.na(log_alpha)){
      #print('na log_alpha value:')
      #print(log_alpha)
      #print('Y value:')
      #print(Y)
    #}
    #if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
    if(log(U[i]) < log_alpha) {
      r0_vec[i] <- Y
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
      count_reject = count_reject + 1
      log_alpha = 0 
    }
    
    #Scaling factor
    #print('log_alpha')
    #print(log_alpha)
    scaling_vec[i] = scaling_vec[i-1] + (1/i)*(exp(log_alpha) - alpha_star)
    print('scaling_vec')
    print(scaling_vec[i])
    
    #Adaptive MC
    #if (i == burn_in){
    #  sigma = var(r0_vec[2:i])*(2.38^2)
    #}
    
  }
  #Final stats
  total_iters = count_accept + count_reject
  accept_rate = 100*(count_accept/(count_accept+count_reject))
  num_samples = count_accept
  print("Total iterations = ")
  print(total_iters)
  print("Acceptance rate = ")
  print(accept_rate)
  print("Number samples = ")
  print(count_accept)
  
  #Burn-in
  r0_vec = r0_vec[burn_in:n]
  
  #Return r0, acceptance rate
  return(list(r0_vec, accept_rate, num_samples))
}

#Apply
alpha_star = 0.40
as_params = adaptive_scaling_metropolis_r0(data, n, sigma, alpha_star, x0 = 1, burn_in = 5000)

r0_as = as_params[1]
r0_as = unlist(r0_as)

#Plot
plot.ts(r0_as)
