#Description
#Super spreading model- Monte Carlo
library(MASS)

#Setup
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")
#par(mar=c(1,1,1,1))

#Params
n = 100 # 10000 #50000
burn_in = 2
shape_gamma = 6
scale_gamma = 1
#Priors
prior_alpha_k = 1
prior_alpha_theta = 1

#********************************************************************
#MCMC Super-spreading
mcmc_super_spreading <- function(data, n, sigma,  x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n)
  
  alpha_vec[1] <- x0; beta_vec[1] <- x0; gamma_vec[1] <- x0
  
  count_accept1 = 0; count_reject1 = 0; count_accept2 = 0
  count_reject2 = 0; count_accept3 = 0; count_reject3 = 0
  
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
      count_reject1 = count_reject1 + 1
    }
    
    #Adaptive MC
    #if (i == burn_in){
    #  sigma1 = var(alpha_vec[2:i])*(2.38^2)
      #print('Burn in reached')
    #}
    
    #************************************************************************
    #beta
    beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma) 
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
      count_reject2 = count_reject2 + 1
    }
    
    #Adaptive MC
    #if (i == burn_in){
    #  sigma2 = var(beta_vec[2:i])*(2.38^2)
    #}
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma) 
    
    if(gamma_dash < 0){
      gamma_dash = abs(gamma_dash)
    }
    
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    prior1 = dgamma(gamma_dash, shape = 1, scale = 1, log = TRUE)
    prior2 = dgamma(gamma_vec[i-1], shape = 1, scale = 1, log = TRUE)
    log_accept_prob = logl_new - logl_prev #+ prior1 - prior2 
    
    #print
    #cat("logl_new: ", logl_new, "\n")
    #cat("logl_prev: ", logl_prev, "\n")
    #cat("prior1: ", prior1, "\n")
    #cat("prior2: ", prior2, "\n")
    #cat("log_accept_prob: ", log_accept_prob, "\n")
    #cat("log(U[i]): ", log(U[i]), "\n")
    
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
      count_reject3 = count_reject3 + 1
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
  #alpha_vec = alpha_vec[burn_in:n]
  #beta_vec = beta_vec[burn_in:n]
  #gamma_vec = gamma_vec[burn_in:n]
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, 
              accept_rate1, num_samples1, sigma, 
              accept_rate2, num_samples2, sigma, 
              accept_rate3, num_samples3, sigma))
}

#********
#*Implement
num_days = 50
#lambda params
shape_gamma = 6
scale_gamma = 1
#params
alphaX = 2 #Without ss event, ~r0. 
betaX = 0.05
gammaX = 10
#Epidemic data
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data)

#Time
n = 1000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params_ad = mcmc_super_spreading(sim_data, n, sigma)

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
plot.ts(alpha_mcmc, ylab = 'alpha', main = paste("MCMC Super spreading model, simulated alpha = ", alphaX))
plot.ts(beta_mcmc, ylab = 'beta', main = paste("MCMC Super spreading model, simulated beta = ", betaX))
plot.ts(gamma_mcmc,  ylab = 'gamma', main = paste("MCMC Super spreading model, simulated gamma = ", gammaX))

#alpha mean
alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alphaX))
print(plot2)

#beta mean
beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
plot2 = plot(seq_along(beta_mean), beta_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, True beta = ",betaX))
print(plot2)

#gamma Mean
gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True gamma = ",gammaX))
print(plot2)

#Plot
file_name = 'ss_mcmc_22_10_21_I'
folder_dir_ad = 'Results/super_spreading_events/ss_mcmc_22_10_21_I' #Use date automate
plot_mcmc_super_spreading_to_screen(sim_data, alpha_mcmc, beta_mcmc, gamma_mcmc, alphaX, betaX, gammaX, file_name, folder_dir_ad)

#*************************************************************************************************
#Super spreading mcmc for a range of alpha values 

ss_mcmc_range_alpha  <- function(num_days, list_alpha, sigma, betaX, gammaX, file_name, folder_dir_ad){
  
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
    sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
    print('Data simulated')
    
    #Time
    start_time = Sys.time()
    print('Start time:')
    print(start_time)
    mcmc_params_ad = mcmc_super_spreading(sim_data, n, sigma)
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
    
    list_time_taken[i] = time_elap # round(time_elap, 4)
    i = i + 1
    
    #Apply Plotting
    plot_mcmc_super_spreading(sim_data, alpha_mcmc, beta_mcmc, gamma_mcmc, alphaX, betaX, gammaX, file_name, folder_dir_ad)
    
  }
  
  #Create dataframe
  df_results <- data.frame(
    alpha = list_alpha,
    accept_rate_alpha = unlist(list_accept_rate1),
    accept_rate_beta = unlist(list_accept_rate2),
    accept_rate_gamma = unlist(list_accept_rate3),
    time_taken = unlist(list_time_taken))
  
  print(df_results)
  
  df_results
}

#Apply
num_days = 50
betaX = 0.05
gammaX = 10
sigma = 0.75
file_name = 'ss_model_iter_I_'
folder_dir_ad = 'Results/super_spreaders/ss_model_iter_II_I'
list_alphaX = c(1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5) #c(0.9, 1.25, 1.75, 2.0, 2.5, 3, 3.5, 4.0, 5.0, 8.0) #c(0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)  #c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0) #c(0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
#list_alpha = c(0.5, 0.65, 0.70, 0.75, 0.8, 0.85, 0.95, 1.05, 2.80, 3.05, 3.55, 4.05, 4.55, 5.05, 8.05, 10.05)
df_ss_results = ss_mcmc_range_alpha(num_days, list_alphaX, sigma, betaX, gammaX, file_name, folder_dir_ad)


#*******************************************************************************
#*
#*MCMC - Investigate one parameter at a time
mcmc_ss_investigate <- function(data, n, sigma, alphaX, betaX, gammaX, x0 = 1) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n)
  
  alpha_vec[1] <- alphaX; beta_vec[1] <- betaX; gamma_vec[1] <- x0
  
  count_accept1 = 0; count_reject1 = 0; count_accept2 = 0
  count_reject2 = 0; count_accept3 = 0; count_reject3 = 0
  
  # #MCMC chain
  for(i in 2:n) {

    #******************************************************
    #alpha
    # alpha_dash <- alpha_vec[i-1] + rnorm(1, sd = sigma)
    # #cat("alpha dash: ", alpha_dash, "\n")
    # 
    # if(alpha_dash < 0){
    #   alpha_dash = abs(alpha_dash)
    # }
    # 
    # #log alpha
    # logl_new = log_like_ss_lse(data, alpha_dash, beta_vec[i-1], gamma_vec[i-1])
    # logl_prev = log_like_ss_lse(data, alpha_vec[i-1], beta_vec[i-1], gamma_vec[i-1])
    # prior1 = dgamma(alpha_dash, shape = 1, scale = 1, log = TRUE)
    # prior2 = dgamma(alpha_vec[i-1], shape = 1, scale = 1, log = TRUE)
    # log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
    # 
    # if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
    #   alpha_vec[i] <- alpha_dash
    #   count_accept1 = count_accept1 + 1
    # } else {
    #   alpha_vec[i] <- alpha_vec[i-1]
    #   count_reject1 = count_reject1 + 1
    # }
    alpha_vec[i] <- alpha_vec[i-1]
    
    #************************************************************************
    # #beta
    # beta_dash <- beta_vec[i-1] + rnorm(1, sd = sigma)
    # if(beta_dash < 0){
    #   beta_dash = abs(beta_dash)
    # }
    # 
    # logl_new = log_like_ss_lse(data, alpha_vec[i], beta_dash, gamma_vec[i-1])
    # logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    # prior1 = dgamma(beta_dash, shape = 1, scale = 1, log = TRUE)
    # prior2 = dgamma(beta_vec[i-1], shape = 1, scale = 1, log = TRUE)
    # log_accept_prob = logl_new - logl_prev #+ prior1 - prior2
    # 
    # if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
    #   beta_vec[i] <- beta_dash
    #   count_accept2 = count_accept2 + 1
    # } else {
    #   beta_vec[i] <- beta_vec[i-1]
    #   count_reject2 = count_reject2 + 1
    # }

    beta_vec[i] <- beta_vec[i-1]
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
      count_reject3 = count_reject3 + 1
    }
    
    #gamma_vec[i] <- gamma_vec[i-1]
    
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
              accept_rate1, num_samples1, sigma, 
              accept_rate2, num_samples2, sigma, 
              accept_rate3, num_samples3, sigma))
}

#********
#*Implement
# num_days = 50
# #lambda params
# shape_gamma = 6
# scale_gamma = 1
# #params
# alphaX = 2 #Without ss event, ~r0. 
# betaX = 0.05
# gammaX = 10
# #Epidemic data
# sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#Time 
n = 5000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
mcmc_params_ad = mcmc_ss_investigate(sim_data, n, sigma, alphaX, betaX,gammaX)

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
plot.ts(alpha_mcmc, ylab = 'alpha', main = paste("MCMC Super spreading model, simulated alpha = ", alphaX))
plot.ts(beta_mcmc, ylab = 'beta', main = paste("MCMC Super spreading model, simulated beta = ", betaX))
plot.ts(gamma_mcmc,  ylab = 'gamma', main = paste("MCMC Super spreading model, simulated gamma = ", gammaX))

#alpha mean
alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alphaX))
print(plot2)

#beta mean
beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
plot2 = plot(seq_along(beta_mean), beta_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, True beta = ",betaX))
print(plot2)

#gamma Mean
gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True gamma = ",gammaX))
print(plot2)




#**************************************************
#*Multi paramater MCMC 

mcmc_super_spreading_multi_var <- function(data, n, x0=as.vector(c(1,1,1)), Sigma = diag(1, nrow = 3)) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Setup
  params <- matrix(0, nrow=3, ncol= n)
  params[,1] <- x0
  U <- runif(n)
  count_accept = 0
  count_reject = 0
  
  #MCMC chain
  for(i in 2:n) {
    
    #Sample new value
    params_new = params[, i-1] + mvrnorm(mu=as.vector(c(0,0,0)), Sigma=Sigma)
    params_new = abs(params_new)
    
    #Acceptance probability
    a = log_like_ss_lse(data, params_new[1], params_new[2], params_new[3])
    b = log_like_ss_lse(data, params[1, i-1], params[2, i-1], params[3, i-1])
    c = dgamma(params_new[1], shape = 1, scale = 1, log = TRUE)
    d = dgamma( params[1, i-1], shape = 1, scale = 1, log = TRUE) 
    e = dgamma(params_new[2], shape = 1, scale = 1, log = TRUE)
    f = dgamma( params[2, i-1], shape = 1, scale = 1, log = TRUE)
    g = dgamma(params_new[3], shape = 1, scale = 1, log = TRUE)
    h = dgamma(params[3, i-1], shape = 1, scale = 1, log = TRUE) 
    
    log_accept_prob = a - b + c - d + e -f + g -h
    
    #Print individual
    # cat("logl_new: ", a, "\n")
    # cat("logl_prev: ", b, "\n")
    # cat("prior_a1: ", c, "\n")
    # cat("prior_a2: ", d, "\n")
    # cat("prior_b1: ", e, "\n")
    # cat("prior_b2: ", f, "\n")
    # cat("prior_c1: ", g, "\n")
    # cat("prior_c2: ", f, "\n")
    # cat("log_accept_prob: ", log_accept_prob, "\n")
    # cat("log(U[i]): ", log(U[i]), "\n")
    
    # log_accept_prob = log_like_ss_lse(data, params_new[1], params_new[2], params_new[3])
    # - log_like_ss_lse(data, params[1, i-1], params[2, i-1], params[3, i-1])
    # + dgamma(params_new[1], shape = 1, scale = 1, log = TRUE)
    # - dgamma( params[1, i-1], shape = 1, scale = 1, log = TRUE) 
    # + dgamma(params_new[2], shape = 1, scale = 1, log = TRUE)
    # - dgamma( params[2, i-1], shape = 1, scale = 1, log = TRUE) 
    # + dgamma(params_new[3], shape = 1, scale = 1, log = TRUE)
    # - dgamma( params[3, i-1], shape = 1, scale = 1, log = TRUE) 
    # cat('log accept prob;', log_accept_prob)
    
    
    #Accept/reject
    if(!(is.na(log_accept_prob)) && log(U[i]) < log_accept_prob) {
      params[,i] <- params_new
      count_accept = count_accept + 1
    } else {
      params[,i] <- params[,i-1]
      count_reject = count_reject + 1
    }
    
    #Sigma
    #if (i == burn_in){
    #  sigma3 = var(gamma_vec[2:i])*(2.38^2)
    #}
    
  }
  
  #Params
  alpha_vec = params[1,]
  beta_vec = params[2,]
  gamma_vec = params[3,]
  
  #Final stats
  total_iters = count_accept + count_reject
  accept_rate = 100*(count_accept/(count_accept + count_reject))
  num_samples = count_accept
  cat("Acceptance rate = ", accept_rate, '\n')
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, accept_rate, num_samples))
}

#Apply
n = 50000
start_time = Sys.time()
print('Start time:')
print(start_time)

mcmc_params_mv = mcmc_super_spreading_multi_var(data, n)

end_time = Sys.time()
time_elap = end_time - start_time
print('Time elapsed:')
print(time_elap)

#Extract params
alpha_mcmc2 = mcmc_params_mv[1]
alpha_mcmc2 = unlist(alpha_mcmc2)

beta_mcmc2 = mcmc_params_mv[2]
beta_mcmc2 = unlist(beta_mcmc2)

gamma_mcmc2 = mcmc_params_mv[3]
gamma_mcmc2 = unlist(gamma_mcmc2)

#Plot
plot.ts(alpha_mcmc2)
plot.ts(beta_mcmc2)
plot.ts(gamma_mcmc2)

#**************************************************************************************************
#Apply multiple times

ss_mcmc_mv_range_alpha  <- function(n, list_alpha, betaX, gammaX, file_name, folder_dir_ad){
  
  #Create folder
  ifelse(!dir.exists(file.path(folder_dir_ad)), dir.create(file.path(folder_dir_ad)), FALSE)
  list_accept_rate = vector('numeric', length(list_alpha))
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
    mcmc_params_ad = mcmc_super_spreading_multi_var(data, n)
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
    list_accept_rate[i] = round(accept_rate_alpha, 2)
    
    num_samples1 = mcmc_params_ad[[5]]
    list_num_samp[i] = num_samples1
    
    list_time_taken[i] = round(time_elap, 2)
    i = i + 1
    
    #Apply Plotting
    plot_mcmc_super_spreading(alpha_mcmc, beta_mcmc, gamma_mcmc, alphaX, betaX, gammaX, file_name, folder_dir_ad)
    
  }
  
  #Create dataframe
  df_results <- data.frame(
    alpha = list_alpha,
    accept_rate_alpha = unlist(list_accept_rate),
    n_samples_alpha = unlist(list_num_samp),
    time_sec = unlist(list_time_taken))
  
  print(df_results)
  
  df_results
}

#Apply
betaX = 3
gammaX = 3
file_name = 'ss_model_iter_I'
folder_dir_ad = 'Results/super_spreaders/ss_model_mcmc_mv_II'
list_alphaX = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5) #c(0.9, 1.25, 1.75, 2.0, 2.5, 3, 3.5, 4.0, 5.0, 8.0) #c(0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)  #c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0) #c(0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
#list_alpha = c(0.5, 0.65, 0.70, 0.75, 0.8, 0.85, 0.95, 1.05, 2.80, 3.05, 3.55, 4.05, 4.55, 5.05, 8.05, 10.05)
df_ss_results = ss_mcmc_mv_range_alpha(n, list_alphaX, betaX, gammaX, file_name, folder_dir_ad)

a = 1; b = 2; c = 4
