#Description
#Adaptive Monte Carlo targeting alpha, beta, gamma
library(MASS)

#Setup
setwd("~/GitHub/epidemic_modelling")
source("functions.R")
#par(mar=c(1,1,1,1))

#Params
n = 100 # 10000 #50000
burn_in = 2
num_days = 10 #60 #100
shape_gamma = 6
scale_gamma = 1
#Priors
prior_alpha_k = 1
prior_alpha_theta = 1

#********************************************************************
#Adaptive MCMC
mcmc_super_spreading <- function(data, n, sigma1, sigma2, sigma3, burn_in, x0 = 1) { #burn_in = 2500
  
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
    
    #******************************************************
    #alpha
    alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma1) 
    #cat("alpha dash: ", alpha_dash, "\n")
    
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev + prior1 - prior2
    #log_accept_prob = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    #- log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    #+ dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    #- dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    
    #cat("logl_new: ", logl_new, "\n")
    #cat("logl_prev: ", logl_prev, "\n")
    #cat("prior1: ", prior1, "\n")
    #cat("prior2: ", prior2, "\n")
    #cat("log_accept_prob: ", log_accept_prob, "\n")
    #cat("log(U[i]): ", log(U[i]), "\n")
    
    if(!(is.na(log_accept_prob)) && log(U[i]) < log_accept_prob) {
      alpha_vec[i] <- alpha_dash
      count_accept1 = count_accept1 + 1
    } else {
      alpha_vec[i] <- alpha_vec[i-1]
      count_reject1 = count_reject1 + 1
    }
    
    ##cat
    #cat("alpha dash: ", alpha_dash, "\n")
    #cat("alpha[i]: ",  alpha_vec[i], "\n")
    
    #Adaptive MC
    #if (i == burn_in){
    #  sigma1 = var(alpha_vec[2:i])*(2.38^2)
      #print('Burn in reached')
    #}
    
    #************************************************************************
    #beta
    #cat("i: ", i, "\n")
    #cat("beta_vec[i-1]: ", beta_vec[i-1], "\n")
    #cat("rnorm(1, sd = sigma2) ", rnorm(1, sd = sigma2), "\n")
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma2) 
    #cat("Beta dash: ", beta_dash, "\n")
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    prior1 = dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev + prior1 - prior2 
      
    #log_accept_prob = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1]) 
    #- log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    #+ dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    #- dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE) #Do other priors cancel?
    
    #Print
    #cat("beta logl_new: ", logl_new, "\n")
    #cat("beta logl_prev: ", logl_prev, "\n")
    #cat("beta prior1: ", prior1, "\n")
    #cat("beta prior2: ", prior2, "\n")
    #cat("beta log_accept_prob: ", log_accept_prob, "\n")
    #cat("beta log(U[i]): ", log(U[i]), "\n")
    
    
    if(!(is.na(log_accept_prob)) && log(U[i]) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
      count_reject2 = count_reject2 + 1
    }
    
    #New value
    #cat("beta dash: ", beta_dash, "\n")
    #cat("beta[i]: ",  beta_vec[i], "\n")
    
    #Adaptive MC
    #if (i == burn_in){
    #  sigma2 = var(beta_vec[2:i])*(2.38^2)
    #}
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma3) 
    #cat("Gamma dash: ", gamma_dash, "\n")
    
    if(gamma_dash < 0){
      gamma_dash = abs(gamma_dash)
    }
    
    #Acceptance Probability
    log1_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    prior1 = dgamma(gamma_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(gamma_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = log1_new - logl_prev + prior1 - prior2 
    
    #Print
    #cat("gamma logl_new: ", logl_new, "\n")
    #cat("gamma logl_prev: ", logl_prev, "\n")
    #cat("gamma prior1: ", prior1, "\n")
    #cat("gamma prior2: ", prior2, "\n")
    #cat("gamma log_accept_prob: ", log_accept_prob, "\n")
    #cat("gamma log(U[i]): ", log(U[i]), "\n")
    
    #log_accept_prob = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash)  
    #- log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    #+ dgamma(gamma_dash, shape = 1, scale = 1, log = TRUE)
    #- dgamma(gamma_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    
    if(!(is.na(log_accept_prob)) && log(U[i]) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
      count_reject3 = count_reject3 + 1
    }
    
    #New value
    #cat("gamma dash: ", beta_dash, "\n")
    #cat("gamma_vec[i]: ",  gamma_vec[i], "\n")
    
    #Adaptive MC
    #if (i == burn_in){
    #  sigma3 = var(gamma_vec[2:i])*(2.38^2)
    #}
    
    
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
  #alpha_vec = alpha_vec[burn_in:n]
  #beta_vec = beta_vec[burn_in:n]
  #gamma_vec = gamma_vec[burn_in:n]
  
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

#Time
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params_ad = mcmc_super_spreading(data, n, sigma, sigma, sigma, burn_in)

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
plot.ts(alpha_mcmc)
plot.ts(beta_mcmc)
plot.ts(gamma_mcmc)


#*************************************************************************************************
#Super spreading mcmc for a range of alpha values 

ss_mcmc_range_alpha  <- function(list_alpha, sigma1, sigma2, sigma3, betaX, gammaX, folder_dir_ad){
  
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
    
    cat('alpha:', alphaX, "\n")
    
    #Get simulated data when alpha is alphaX
    data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
    print('Data simulated')
    
    #Time
    start_time = Sys.time()
    print('Start time:')
    print(start_time)
    mcmc_params_ad = mcmc_super_spreading(data, n, sigma1, sigma2, sigma3, burn_in)
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
    plot_mcmc_super_spreading(alpha_mcmc, beta_mcmc, gamma_mcmc, alphaX, betaX, gammaX, folder_dir_ad)
    
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
folder_dir_ad = 'Results/super_spreaders/ss_model_mcmc_results_II'
list_alphaX = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5) #c(0.9, 1.25, 1.75, 2.0, 2.5, 3, 3.5, 4.0, 5.0, 8.0) #c(0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)  #c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0) #c(0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
#list_alpha = c(0.5, 0.65, 0.70, 0.75, 0.8, 0.85, 0.95, 1.05, 2.80, 3.05, 3.55, 4.05, 4.55, 5.05, 8.05, 10.05)
df_ss_results = ss_mcmc_range_alpha(list_alphaX, sigma, sigma, sigma, betaX, gammaX, folder_dir_ad)


#**************************************************
#*Multi paramater MCMC 

mcmc_super_spreading_multi_var <- function(data, n, burn_in, x0=as.vector(c(0,0, 0)), Sigma = diag(1, nrow = 3)) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Setup
  params <- matrix(0, nrow=3, ncol= n)
  params[,1] <- x0
  U <- runif(n)
  count_accept1 = 0
  count_reject1 = 0
  count_reject3 = 0
  sd_sample = 1
  
  #MCMC chain
  for(i in 2:n) {
    
    #Sample new value
    params_new = params[, i-1] + mvrnorm(mu=as.vector(c(0,0,0)), Sigma=Sigma)
    params_new = abs(params_new)
    
    #Acceptance probability
    log_accept_prob = log_like_ss_lse(data, params_new[1], params_new[2], params_new[3])
    - log_like_ss_lse(data, params[1, i-1], params[2, i-1], params[3, i-1])
    + dgamma(params_new[1], shape = 1, scale = 1, log = TRUE)
    - dgamma( params[1, i-1], shape = 1, scale = 1, log = TRUE) 
    + dgamma(params_new[2], shape = 1, scale = 1, log = TRUE)
    - dgamma( params[2, i-1], shape = 1, scale = 1, log = TRUE) 
    + dgamma(params_new[3], shape = 1, scale = 1, log = TRUE)
    - dgamma( params[3, i-1], shape = 1, scale = 1, log = TRUE) 

    #Accept/reject
    if(!(is.na(log_accept_prob)) && log(U[i]) < log_accept_prob) {
      params[,i] <- params_new
      count_accept = count_accept + 1
    } else {
      params[,i] <- params[,i-1]
      count_reject = count_reject + 1
    }
    
    #Adaptive MC
    #if (i == burn_in){
    #  sigma3 = var(gamma_vec[2:i])*(2.38^2)
    #}
    
  }
  
  #Final stats
  total_iters = count_accept + count_reject
  accept_rate = 100*(count_accept/(count_accept + count_reject))
  num_samples = count_accept
  cat("Acceptance rate = ", accept_rate, '\n')
  
  #Return alpha, acceptance rate
  return(list(params, accept_rate, num_samples))
}

#Apply
#Time
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params_ad = mcmc_super_spreading(data, n, sigma, sigma, sigma, burn_in)

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
plot.ts(alpha_mcmc)
plot.ts(beta_mcmc)
plot.ts(gamma_mcmc)