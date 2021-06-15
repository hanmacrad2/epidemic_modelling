#Description
#Adaptive Monte Carlo targeting r0
library(gridExtra)

#Setup
source("1_simulation.R")
par(mar=c(1,1,1,1))

#Params
num_days = 60 #100
shape_gamma = 6
scale_gamma = 1
#Priors
prior_r0_k = 1
prior_r0_theta = 1


#***********************************
#Log Likelihood 
log_like_ss <- function(x, r01, r02, p){
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda1 = r01*sum(x[1:t-1]*rev(prob_infect[1:t-1]))
    lambda2 = r02*sum(x[1:t-1]*rev(prob_infect[1:t-1]))
    logl = logl + log(p) + log(lambda1) + log(1-p) + log(lamda2) - lambda1 - lambda2
    
  }
  
  logl
  
}


#********************************************************************
#Adaptive MCMC
adaptive_mc_r0_ss <- function(data, n, sigma, x0 = 1, burn_in = 5000) { #burn_in = 2500
  
  'Returns mcmc samples of R0 & acceptance rate'
  
  #Set up
  r01_vec <- vector('numeric', n)
  r02_vec <- vector('numeric', n)
  p_vec <- vector('numeric', n)
  
  r01_vec[1] <- x0
  r02_vec[1] <- x0
  p_vec[1] <- x0
  
  U <- runif(n)
  count_accept = 0
  count_reject = 0
  sd_sample = 1
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************
    #ro1
    ro1_dash <- r01_vec[i-1] + rnorm(1, sd = sigma1) 
    if(Y < 0){
      Y = abs(Y)
    }
    
    log_alpha = log_like_ss(data, ro1_dash, r02_vec[i-1], p_vec[i-1]) - log_like_ss(data, r01_vec[i-1], r02_vec[i-1], p_vec[i-1]) + dgamma(ro1_dash, shape = 1, scale = 1, log = TRUE) - dgamma(r01_vec[i-1], shape = 1, scale = 1, log = TRUE) #log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      r01_vec[i] <- ro1_dash
      count_accept1 = count_accept1 + 1
    } else {
      r01_vec[i] <- r01_vec[i-1]
      count_reject1 = count_reject1 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma1 = var(r01_vec[2:i])*(2.38^2)
    }
    
    #******************
    #ro2
    ro2_dash <- r02_vec[i-1] + rnorm(1, sd = sigma2) 
    if(ro2_dash < 0){
      ro2_dash = abs(ro2_dash)
    }
    
    log_alpha = log_like_ss(data, r01_vec[i], ro2_dash, p_vec[i-1]) - log_like_ss(data, r01_vec[i], r02_vec[i-1], p_vec[i-1]) + dgamma(ro1_dash, shape = 1, scale = 1, log = TRUE) - dgamma(r01_vec[i-1], shape = 1, scale = 1, log = TRUE)  
    
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      r02_vec[i] <- ro2_dash
      count_accept2 = count_accept2 + 1
    } else {
      r02_vec[i] <- r02_vec[i-1]
      count_reject2 = count_reject2 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma2 = var(r02_vec[2:i])*(2.38^2)
    }
    
    #******************
    #p
    p_dash <- p_vec[i-1] + rnorm(1, sd = sigma3) 
    if(p_dash < 0){
      p_dash = abs(p_dash)
    }
    
    log_alpha = log_like(data, r01_vec[i], r02_vec[i-1], p_dash) - log_like(data, r01_vec[i], r02_vec[i-1], p_vec[i-1]) + dgamma(p_dash, shape = 1, scale = 1, log = TRUE) - dgamma(p_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      p_vec[i] <- p_dash
      count_accept3 = count_accept3 + 1
    } else {
      p_vec[i] <- p_vec[i-1]
      count_reject2 = count_reject2 + 1
    }
    
    #Adaptive MC
    if (i == burn_in){
      sigma3 = var(p_vec[2:i])*(2.38^2)
    }
  }
  #Final stats
  total_iters1 = count_accept1 + count_reject1
  accept_rate1 = 100*(count_accept1/(count_accept1+count_reject1))
  num_samples1 = count_accept1
  print("Acceptance rate1 = ")
  print(accept_rate1)
  
  #Burn-in 
  r01_vec = r01_vec[burn_in:n]
  r02_vec = r02_vec[burn_in:n]
  p_vec = p_vec[burn_in:n]
  
  #Return r0, acceptance rate
  return(list(r01_vec, r02_vec, p_vec, accept_rate1, num_samples1, sigma1, accept_rate2, num_samples2, sigma2, accept_rate3, num_samples3, sigma3))
}

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

#Plot





#************************
#Other

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
folder_dir_ad = 'Results/Adaptive_MC/adaptive_mc_iter_VII_burn_in_5k'
list_r0 = c(0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)  #c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0) #c(0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
#list_r0 = c(0.5, 0.65, 0.70, 0.75, 0.8, 0.85, 0.95, 1.05, 2.80, 3.05, 3.55, 4.05, 4.55, 5.05, 8.05, 10.05)
df_ad_results_formI = apply_adaptive_mc_range_r0(list_r0, sigma, folder_dir_ad)

