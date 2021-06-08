#Description
#Metropolis algorithm with target r0
library(gridExtra)

#Setup
#source("1_simulate_branching.R")
source("1_simulation.R")
par(mar=c(1,1,1,1))

#Params
num_days = 60 #100
r0_true = 2.5 #3.1 #2.8 #
shape_gamma = 6
scale_gamma = 1
data = simulate_branching(num_days, r0_true, shape_gamma, scale_gamma)

#Priors
prior_r0_k = 1
prior_r0_theta = 1


#***********************************
#Likelihood (log)
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

#***********************************
#MCMC
MetropolisHastings_r0 <- function(data, n, sigma, x0 = 1, burn_in = 2500) {
  
  'Returns mcmc samples of R0'
  
  #Set up
  r0_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  U <- runif(n)
  count_accept = 0
  count_reject = 0
  
  #MCMC chain
  for(i in 2:n) {
    Y <- r0_vec[i-1] + rnorm(1, sd = sigma) #, mean = 0, sd = sigma_opt)
    if(Y < 0){
      Y = abs(Y)
    }
    
    #Alpha
    log_alpha = log_like(data, Y) - log_like(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    #Should include: Likelihood + prior + propogsal density x2 (Previous time step & Current time step)
    
    if (is.na(log_alpha)){
      print('na value')
      sprintf("Y: %i", Y)
    }
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      r0_vec[i] <- Y
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
      count_reject = count_reject + 1
    }
  }
  #Final stats
  total_iters = count_accept + count_reject
  accept_rate = 100*(count_accept/(count_accept+count_reject))
  print("Total iterations = ")
  print(total_iters)
  print("Acceptance rate = ")
  print(accept_rate)
  
  r0_vec = r0_vec[burn_in:n]
  r0_vec
}

#Apply MCMC
n = 50000
sigma =  0.5 #(2.38^2/dimesnion_paramter)*Posterior or sample variance ~optimal
r0_mcmc = MetropolisHastings_r0(data, n, sigma)

#Plots
mcmc_plotting <- function(mcmc_vector, r0_true) {

  #i. MCMC chain
  ts.plot(mcmc_vector, ylab = 'R0', main = paste("MCMC of R0, true R0 = ", r0_true, "sd of proposal = ", sigma))
  #main = paste("Tau= ", tau, " Sigma= ", sigma),
  
  #ii. Mean
  #Plot mean
  r0_mean = cumsum(mcmc_vector)/seq_along(mcmc_vector)
  plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0', main = paste("Mean of R0 MCMC chain, True R0 = ",r0_true, ", sd of proposal = ", sigma))
  
  #Histogram
  hist(mcmc_vector, prob = TRUE)
  
  #Hist
  hist1 <- hist(mcmc_vector, breaks = 80)
  hist1$counts <- hist1$counts/sum(hist1$counts)
  plot(hist1, xlab = 'r0', ylab = 'Density', 
       main = 'Empirical density of r0 - MCMC chain')

}
#*******************************************************************
#MCMC v2
MetropolisHastings_r0_vII <- function(data, n, sigma, x0 = 1, burn_in = 2500) {
  
  'Returns mcmc samples of R0 & acceptance rate'
  
  #Set up
  r0_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  U <- runif(n)
  count_accept = 0
  count_reject = 0
  
  #MCMC chain
  for(i in 2:n) {
    Y <- r0_vec[i-1] + rnorm(1, sd = sigma) #, mean = 0, sd = sigma_opt)
    if(Y < 0){
      Y = abs(Y)
    }
    
    log_alpha = log_like(data, Y) - log_like(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) #log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    
    if (is.na(log_alpha)){
      print('na value, Y value:')
      print(Y)
    }
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      r0_vec[i] <- Y
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
      count_reject = count_reject + 1
    }
  }
  #Final stats
  total_iters = count_accept + count_reject
  print("Acceptance count = ")
  print(count_accept)
  accept_rate = (count_accept/(count_accept+count_reject))
  print("Total iterations = ")
  print(total_iters)
  print("Acceptance rate = ")
  print(accept_rate)
  
  #Burn-in
  r0_vec = r0_vec[burn_in:n]
  
  #Return r0, acceptance rate
  return(list(r0_vec, accept_rate))
}

#Apply
mcmc_params = MetropolisHastings_r0_vII(data, n, sigma)
r0X = mcmc_params[1]
r0X = unlist(r0X)

#******************************************************************************
# Investigate affect of varying the sd in the Standard Normal proposal Distribution

MCMC_range_sd <- function(list_sd, data, n, r0_true){
  
  #Initialise
  par(mfrow = c(2, 2))
  list_accept_rate <- vector('numeric', length(list_sd))
  i = 1
  
  #Create folder
  folder_dir = paste("Results/mcmc_results_r0_true_", r0_true)
  ifelse(!dir.exists(file.path(folder_dir)), dir.create(file.path(folder_dir)), FALSE)
  
  for (sigma in list_sd) {
    print('sigma =')
    print(sigma)
    #R0 - mcmc
    mcmc_params = MetropolisHastings_r0_vII(data, n, sigma)
    #Extract
    r0_mcmc = mcmc_params[1]
    r0_mcmc = unlist(r0_mcmc)
    accept_rate = mcmc_params[2]
    print('Accept_rate =')
    print(accept_rate)
    list_accept_rate[i] = accept_rate
    
    #Plots
    
    #i.MCMC chain
    pdf(paste(folder_dir, "/", i, "_mcmc_chain_sd_", sigma, "_true_r0_", r0_true, ".pdf", sep="")) #, width=6, height=4)
    plot1 = ts.plot(r0_mcmc, ylab = 'R0', 
            main = paste("MCMC of R0. True R0 = ", r0_true, ". Sd of proposal = ", sigma))
    print(plot1)
    
    #ii.Cumulative Mean
    r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
    plot2 = plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0',
                 main = paste("Mean of R0 MCMC chain. True R0 = ", r0_true, ". sd of proposal = ", sigma))
    print(plot2)
    
    #iii. Histograms
    hist1 <- hist(r0_mcmc, breaks = 80, xlab = 'r0', ylab = 'Density', 
                  main = paste("Histogram of R0. True R0 = ", r0_true, ". sd of proposal = ", sigma))
    print(hist1)
    #iiib. Histogram - density
    hist1$counts <- hist1$counts/sum(hist1$counts)
    hist2 = plot(hist1, xlab = 'r0', ylab = 'Density', 
         main = paste("Empirical density of r0. True R0 = ", r0_true, ". sd of proposal = ", sigma))
    print(hist2)
    dev.off()
    
    #Increment
    i = i + 1
  }
  
  #Create dataframe
  df_sd_results <- data.frame(
    sd = list_sd,
    acceptance_rate = unlist(list_accept_rate))
  
  print(df_sd_results) 
  
  #Save results
  #pdf(paste(folder_dir, "/df_acceptance_rate_r0_true_",r0_true, ".pdf", sep="")) #pdf("test.pdf", height=11, width=10)
  #grid.table(df_sd_results)
  #dev.off()
  
  df_sd_results
}

#Apply
r0_true = 3.2 #2.8
data = simulate_branching(num_days, r0_true, shape_gamma, scale_gamma)
list_sd = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 5, 7.5, 10)
df_sd_mcmc_results = MCMC_range_sd(list_sd, data, n, r0_true)


#********************************************************************
#Adaptive MCMC

adaptive_mc_r0 <- function(data, n, x0 = 1, burn_in = 2500) {
  
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
    Y <- r0_vec[i-1] + rnorm(1, sd = sd_sample) 
    if(Y < 0){
      Y = abs(Y)
    }
    
    log_alpha = log_like(data, Y) - log_like(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) #log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    
    if (is.na(log_alpha)){
      #print('na value, Y value:')
      #print(Y)
    }
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      r0_vec[i] <- Y
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
      count_reject = count_reject + 1
    }
    #Get sample sd
    sd_sample = var(r0_vec)*(2.38^2)
  }
  #Final stats
  total_iters = count_accept + count_reject
  accept_rate = 100*(count_accept/(count_accept+count_reject))
  print("Total iterations = ")
  print(total_iters)
  print("Acceptance rate = ")
  print(accept_rate)
  
  #Burn-in
  r0_vec = r0_vec[burn_in:n]
  
  #Return r0, acceptance rate
  return(list(r0_vec, accept_rate))
}


adaptive_mc_r0I <- function(data, n, x0 = 1, burn_in = 2500) {
  
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
    Y <- r0_vec[i-1] + rnorm(1, sd = sd_sample) 
    if(Y < 0){
      Y = abs(Y)
    }
    
    log_alpha = log_like(data, Y) - log_like(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) #log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    
    if (is.na(log_alpha)){
      print('na value, Y value:')
      print(Y)
    }
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
      r0_vec[i] <- Y
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
      count_reject = count_reject + 1
    }
    #Get sample sd
    sd_sample = sd(r0_vec)
  }
  #Final stats
  total_iters = count_accept + count_reject
  accept_rate = 100*(count_accept/(count_accept+count_reject))
  print("Total iterations = ")
  print(total_iters)
  print("Acceptance rate = ")
  print(accept_rate)
  
  #Burn-in
  r0_vec = r0_vec[burn_in:n]
  
  #Return r0, acceptance rate
  return(list(r0_vec, accept_rate))
}

#Plots
mcmc_plotting_adaptive <- function(mcmc_vector, r0_true, folder_dir_ad) {
  
  #Folder save
  pdf(paste(folder_dir_ad, "/", "adpative_mc_r0_true_", r0_true, ".pdf", sep=""))
  
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector, ylab = 'R0', main = paste("Adaptive MC for R0, true R0 = ", r0_true, ". SD of proposal = ", sigma))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  r0_mean = cumsum(mcmc_vector)/seq_along(mcmc_vector)
  plot2 = plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0', main = paste("Mean of R0 MCMC chain, True R0 = ",r0_true, ". SD of proposal = ", sigma))
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

#Apply
#Simulated data
r0_true = 2.5 #3.1
data = simulate_branching(num_days, r0_true, shape_gamma, scale_gamma)

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
accept_rate = mcmc_params_ad[2]

#Apply Plotting
folder_dir_ad = 'Results/Adaptive_MC'
mcmc_plotting_adaptive(r0_mcmc, r0_true, folder_dir_ad)

#************
#Apply Adaptive MCMC to a range of r0s

apply_adaptive_mc_range_r0 <- function(list_r0, folder_dir_ad){
  
  #Create folder
  ifelse(!dir.exists(file.path(folder_dir_ad)), dir.create(file.path(folder_dir_ad)), FALSE)
  list_accept_rate = vector('numeric', length(list_r0))
  i = 1
  
  for (r0X in list_r0){
    
    #Get simulated data when r0 is r0X
    data = simulate_branching(num_days, r0X, shape_gamma, scale_gamma)
    
    #Time
    start_time = Sys.time()
    mcmc_params_ad = adaptive_mc_r0(data, n)
    end_time = Sys.time()
    time_elap = end_time - start_time
    print('Time elapsed:')
    print(time_elap)
    
    #Extract params
    r0_mcmc = mcmc_params_ad[1]
    r0_mcmc = unlist(r0_mcmc)
    accept_rate = mcmc_params_ad[2]
    list_accept_rate[i] = accept_rate
    i = i + 1
    
    #Apply Plotting
    mcmc_plotting_adaptive(r0_mcmc, r0X, folder_dir_ad)
    
  }
  
  #Create dataframe
  df_results <- data.frame(
    r0 = list_r0,
    acceptance_rate = unlist(list_accept_rate))
  
  print(df_results)
  
  df_results
}

#Apply
folder_dir_ad = 'Results/adaptive_mc_formula_iter_V'
#list_r0 = c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)  #c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0) #c(0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
list_r0 = c(0.65, 0.70, 0.75, 0.8, 0.85, 0.95, 1.05, 2.80, 3.05, 3.55, 4.05, 4.55, 5.05, 8.05, 10.05)
df_ad_results_formI = apply_adaptive_mc_range_r0(list_r0, folder_dir_ad)


#**************
#Try sd = sd(sample)
apply_adaptive_mc_range_r0_I <- function(list_r0, folder_dir_ad){
  
  #Create folder
  ifelse(!dir.exists(file.path(folder_dir_ad)), dir.create(file.path(folder_dir_ad)), FALSE)
  list_accept_rate = vector('numeric', length(list_r0))
  i = 1
  
  for (r0X in list_r0){
    
    #Get simulated data when r0 is r0X
    data = simulate_branching(num_days, r0X, shape_gamma, scale_gamma)
    
    #Time
    start_time = Sys.time()
    mcmc_params_ad = adaptive_mc_r0I(data, n)
    end_time = Sys.time()
    time_elap = end_time - start_time
    print('Time elapsed:')
    print(time_elap)
    
    #Extract params
    r0_mcmc = mcmc_params_ad[1]
    r0_mcmc = unlist(r0_mcmc)
    accept_rate = mcmc_params_ad[2]
    print('Accept rate')
    print(accept_rate)
    list_accept_rate[i] = accept_rate
    i = i + 1
    
    #Apply Plotting
    mcmc_plotting_adaptive(r0_mcmc, r0X, folder_dir_ad)
    
  }
  
  #Create dataframe
  df_results <- data.frame(
    r0 = list_r0,
    acceptance_rate = unlist(list_accept_rate))
  
  print(df_results)
  
  df_results
}

#Apply
folder_dir_ad = 'Results/adaptive_mc_sd_sample_iterIII'
list_r0 = c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)  #c(0.8, 0.9, 1.0, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0) #c(0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
df_ad_results_formI = apply_adaptive_mc_range_r0(list_r0, folder_dir_ad)


#**************
#*Superspreader - new file (New file for Adaptive too :) )
