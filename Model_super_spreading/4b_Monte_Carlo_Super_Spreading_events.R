#Description
#Adaptive Monte Carlo targeting alpha, beta, gamma
cat("Current working dir: ", 2)

#Setup
setwd("~/GitHub/epidemic_modelling")
source("functions.R")
#par(mar=c(1,1,1,1))

#Params
n = 50000
num_days = 10 #60 #100
shape_gamma = 6
scale_gamma = 1
#Priors
prior_alpha_k = 1
prior_alpha_theta = 1

#********************************************************************
#Adaptive MCMC
adaptive_mc_ss <- function(data, n, sigma1, sigma2, sigma3, x0 = 1, burn_in = 5000) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n)
  beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n)
  
  alpha_vec[1] <- x0
  beta_vec[1] <- x0
  gamma_vec[1] <- x0
  
  U <- runif(n)
  count_accept1 = 0
  count_reject1 = 0
  count_accept2 = 0
  count_reject2 = 0
  count_accept3 = 0
  count_reject3 = 0
  sd_sample = 1
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************
    #alpha
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma1) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    cat("alpha dash: ", alpha_dash)

    log_alpha = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    - log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    + dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    - dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    
    cat("log_alpha: ", log_alpha)
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
      count_reject1 = count_reject1 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma1 = var(alpha_vec[2:i])*(2.38^2)
      print('Burn in reached')
    }
    
    #******************
    #beta
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma2) 
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    log_alpha = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1]) 
    - log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    + dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    - dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE) #Do other priors cancel?
    
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
      count_reject2 = count_reject2 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma2 = var(beta_vec[2:i])*(2.38^2)
    }
    
    #******************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma3) 
    if(gamma_dash < 0){
      gamma_dash = abs(gamma_dash)
    }
    
    log_alpha = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_dash)  
    - log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    + dgamma(gamma_dash, shape = 1, scale = 1, log = TRUE)
    - dgamma(gamma_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
      count_reject3 = count_reject3 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma3 = var(gamma_vec[2:i])*(2.38^2)
    }
    
    
  }
  #Final stats
  #alpha
  total_iters1 = count_accept1 + count_reject1
  accept_rate1 = 100*(count_accept1/(count_accept1+count_reject1))
  num_samples1 = count_accept1
  print("Acceptance rate1 = ")
  print(accept_rate1)
  
  #beta
  total_iters2 = count_accept2 + count_reject2
  accept_rate2 = 100*(count_accept2/(count_accept2+count_reject2))
  num_samples2 = count_accept2
  print("Acceptance rate2 = ")
  print(accept_rate2)
  
  #gamma
  total_iters3 = count_accept3 + count_reject3
  accept_rate3 = 100*(count_accept3/(count_accept3+count_reject3))
  num_samples3 = count_accept3
  print("Acceptance rate3 = ")
  print(accept_rate3)
  
  #Burn-in 
  alpha_vec = alpha_vec[burn_in:n]
  beta_vec = beta_vec[burn_in:n]
  gamma_vec = gamma_vec[burn_in:n]
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, 
              accept_rate1, num_samples1, sigma1, 
              accept_rate2, num_samples2, sigma2, 
              accept_rate3, num_samples3, sigma3))
}

#********
#*Implement
num_days = 30
num_days = 30 #100
shape_gamma = 6
scale_gamma = 1
alphaX = 3 #Without ss event, ~r0. 
betaX = 3
gammaX = 3
data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
data

#Time
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params_ad = adaptive_mc_ss(data, n, sigma, sigma, sigma)
end_time = Sys.time()
time_elap = end_time - start_time
print('Time elapsed:')
print(time_elap)

#Extract params
alpha_mcmc = mcmc_params_ad[1]
alpha_mcmc = unlist(alpha_mcmc)

beta_mcmc = mcmc_params_ad[2]
beta_mcmc = unlist(beta_mcmc)

gamma_mcmc = mcmc_params_ad[3]
gamma_mcmc = unlist(gamma_mcmc)

#Plot
plot.ts(gamma_mcmc)

#**********************************************
#Plots
mcmc_plotting_adaptive_ss <- function(mcmc_vector1, mcmc_vector2, mcmc_vector3, alpha_true, folder_dir_ad) {
  
  #Folder save
  pdf(paste(folder_dir_ad, "/", "ss_adpative_mc_alpha_true_", alpha_true, ".pdf", sep=""))
  
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector1, ylab = 'alpha', main = paste("Superspreader Adaptive MC for alpha, true alpha = ", alpha_true))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  alpha_mean = cumsum(mcmc_vector1)/seq_along(mcmc_vector1)
  plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alpha_true))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector1, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'alpha', ylab = 'Density', 
               main = 'Empirical density of alpha - MCMC chain')
  print(hist3)
  
  #2. beta
  #**************************************************
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector2, ylab = 'beta', main = paste("Superspreader Adaptive MC for beta, true alpha = ", alpha_true))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  alpha_mean = cumsum(mcmc_vector2)/seq_along(mcmc_vector2)
  plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, True alpha = ",alpha_true))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector2, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'alpha', ylab = 'Density', 
               main = 'Empirical density of beta - MCMC chain')
  print(hist3)
  
  #**************************************************
  #*gamma
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector3, ylab = 'gamma', main = paste("Superspreader Adaptive MC for gamma, true alpha = ", alpha_true))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  gamma_mean = cumsum(mcmc_vector3)/seq_along(mcmc_vector3)
  plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True alpha = ",alpha_true))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector3, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'gamma', ylab = 'Density', 
               main = 'Empirical density of gamma - MCMC chain')
  print(hist3)
  
  dev.off()
  
}

#************************
#Other

#********************************************************************
#Apply Adaptive MCMC to a range of alphas

apply_adaptive_mc_range_alpha_ss <- function(list_alpha, sigma1, sigma2, sigma3, betaX, gammaX, folder_dir_ad){
  
  #Create folder
  ifelse(!dir.exists(file.path(folder_dir_ad)), dir.create(file.path(folder_dir_ad)), FALSE)
  list_accept_rate1 = vector('numeric', length(list_alpha))
  list_accept_rate2 = vector('numeric', length(list_alpha))
  list_accept_rate3 = vector('numeric', length(list_alpha))
  list_sd1 = vector('numeric', length(list_alpha))
  list_sd2 = vector('numeric', length(list_alpha))
  list_sd3 = vector('numeric', length(list_alpha))
  list_num_samp = vector('numeric', length(list_alpha))
  list_time_taken = vector('numeric', length(list_alpha))
  i = 1
  
  for (alphaX in list_alpha){
    
    print('alpha:')
    print(alphaX)
    
    #Get simulated data when alpha is alphaX
    data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
    print('Data simulated')
    
    #Time
    start_time = Sys.time()
    print('Start time:')
    print(start_time)
    mcmc_params_ad = adaptive_mc_ss(data, n, sigma1, sigma2, sigma3)
    end_time = Sys.time()
    time_elap = end_time - start_time
    print('Time elapsed:')
    print(time_elap)
    
    #Extract params
    alpha_mcmc = mcmc_params_ad[1]
    alpha_mcmc = unlist(alpha_mcmc)
    
    beta_mcmc = mcmc_params_ad[2]
    beta_mcmc = unlist(beta_mcmc)
    
    gamma_mcmc = mcmc_params_ad[3]
    gamma_mcmc = unlist(gamma_mcmc)
    
    accept_rate_alpha = mcmc_params_ad[[4]]
    list_accept_rate1[i] = round(accept_rate_alpha, 2)
    
    num_samples1 = mcmc_params_ad[[5]]
    list_num_samp[i] = num_samples1
    
    sd_final_alpha = mcmc_params_ad[[6]]
    list_sd1[i] = round(sd_final_alpha, 3)
    
    accept_rate_beta = mcmc_params_ad[[7]]
    list_accept_rate2[i] = round(accept_rate_beta, 2)
    
    sd_final_beta = mcmc_params_ad[[9]]
    list_sd2[i] = round(sd_final_beta, 3)
    
    accept_rate_gamma = mcmc_params_ad[[10]]
    list_accept_rate3[i] = round(accept_rate_gamma, 2)
    
    sd_final_gamma = mcmc_params_ad[[12]]
    list_sd3[i] = round(sd_final_gamma, 3)
    
    list_time_taken[i] = round(time_elap, 2)
    i = i + 1
    
    #Apply Plotting
    mcmc_plotting_adaptive_ss(alpha_mcmc, beta_mcmc, gamma_mcmc, alphaX, folder_dir_ad)
    
  }
  
  #Create dataframe
  df_results <- data.frame(
    alpha = list_alpha,
    accept_rate_alpha = unlist(list_accept_rate1),
    n_samples_alpha = unlist(list_num_samp),
    sd_final_alpha = unlist(list_sd1),
    accept_rate_beta = unlist(list_accept_rate2),
    sd_final_beta = unlist(list_sd2),
    accept_rate_alpha3 = unlist(list_accept_rate3),
    sd_final_alpha3 = unlist(list_sd3),
    time_sec = unlist(list_time_taken))
  
  print(df_results)
  
  df_results
}

#Apply
betaX = 3
gammaX = 3
sigma = 0.75
folder_dir_ad = 'Results/super_spreaders/ss_model_III_iterIII'
list_alphaX = c(1.0, 2.0) #c(0.9, 1.25, 1.75, 2.0, 2.5, 3, 3.5, 4.0, 5.0, 8.0) #c(0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)  #c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0) #c(0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
#list_alpha = c(0.5, 0.65, 0.70, 0.75, 0.8, 0.85, 0.95, 1.05, 2.80, 3.05, 3.55, 4.05, 4.55, 5.05, 8.05, 10.05)
#df_ss_results = apply_adaptive_mc_range_alpha_ss(list_alphaX, sigma, sigma, sigma, betaX, gammaX, folder_dir_ad)


#***************************************************************************************************************************
#* Adaptive Scaling 

#Adaptive MCMC
adaptive_scaling_met_ss <- function(data, n, sigma1, sigma2, sigma3, x0 = 1, burn_in = 5000) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n)
  beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n)
  
  alpha_vec[1] <- x0
  beta_vec[1] <- x0
  gamma_vec[1] <- x0
  
  U <- runif(n)
  count_accept1 = 0
  count_reject1 = 0
  count_accept2 = 0
  count_reject2 = 0
  count_accept3 = 0
  count_reject3 = 0
  sd_sample = 1
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************
    #alpha
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma1) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    log_alpha = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    - log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    + dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    - dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
      count_reject1 = count_reject1 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma1 = var(alpha_vec[2:i])*(2.38^2)
      print('Burn in reached')
    }
    
    #******************
    #beta
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma2) 
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    log_alpha = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1]) 
    - log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    + dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    - dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE) #Do other priors cancel?
    
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
      count_reject2 = count_reject2 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma2 = var(beta_vec[2:i])*(2.38^2)
    }
    
    #******************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma3) 
    if(gamma_dash < 0){
      gamma_dash = abs(gamma_dash)
    }
    
    log_alpha = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_dash)  
    - log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    + dgamma(gamma_dash, shape = 1, scale = 1, log = TRUE)
    - dgamma(gamma_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
      count_reject3 = count_reject3 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma3 = var(gamma_vec[2:i])*(2.38^2)
    }
    
    
  }
  #Final stats
  #alpha
  total_iters1 = count_accept1 + count_reject1
  accept_rate1 = 100*(count_accept1/(count_accept1+count_reject1))
  num_samples1 = count_accept1
  print("Acceptance rate1 = ")
  print(accept_rate1)
  
  #beta
  total_iters2 = count_accept2 + count_reject2
  accept_rate2 = 100*(count_accept2/(count_accept2+count_reject2))
  num_samples2 = count_accept2
  print("Acceptance rate2 = ")
  print(accept_rate2)
  
  #gamma
  total_iters3 = count_accept3 + count_reject3
  accept_rate3 = 100*(count_accept3/(count_accept3+count_reject3))
  num_samples3 = count_accept3
  print("Acceptance rate3 = ")
  print(accept_rate3)
  
  #Burn-in 
  alpha_vec = alpha_vec[burn_in:n]
  beta_vec = beta_vec[burn_in:n]
  gamma_vec = gamma_vec[burn_in:n]
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, 
              accept_rate1, num_samples1, sigma1, 
              accept_rate2, num_samples2, sigma2, 
              accept_rate3, num_samples3, sigma3))
}

#********
#*Implement
num_days = 30
num_days = 30 #100
shape_gamma = 6
scale_gamma = 1
alphaX = 3 #Without ss event, ~r0. 
betaX = 3
gammaX = 3
data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
data

#Time
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params_ad = adaptive_mc_ss(data, n, sigma, sigma, sigma)
end_time = Sys.time()
time_elap = end_time - start_time
print('Time elapsed:')
print(time_elap)

#Extract params
alpha_mcmc = mcmc_params_ad[1]
alpha_mcmc = unlist(alpha_mcmc)

beta_mcmc = mcmc_params_ad[2]
beta_mcmc = unlist(beta_mcmc)

gamma_mcmc = mcmc_params_ad[3]
gamma_mcmc = unlist(gamma_mcmc)

#Plot
plot.ts(alpha_mcmc, title = 'alpha mcmc')

plot.ts(beta_mcmc, title = 'alpha mcmc')

plot.ts(gamma_mcmc, main = 'gamma mcmc')




