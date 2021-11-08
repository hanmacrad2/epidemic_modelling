#Description
#Super spreading Events model- MCMC Inference
#Setup
library(MASS)
library(pracma)
setwd("~/GitHub/epidemic_modelling/Model_super_spreading")
source("functions.R")

#Epidemic params
num_days = 50
#Gamma params for infectiousness curve (lambda) distribution
shape_gamma = 6
scale_gamma = 1 

seed_count = 1
#par(mar=c(1,1,1,1))


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
    log_accept_prob = logl_new - logl_prev  
    
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
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
      beta_vec[i] <- beta_dash
      count_accept2 = count_accept2 + 1
    } else {
      beta_vec[i] <- beta_vec[i-1]
    }
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma) 
    
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #abs(gamma_dash)
    }
    
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_dash) 
    logl_prev = log_like_ss_lse(data, alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    log_accept_prob = logl_new - logl_prev 
    
    #Metropolis Acceptance Step
    if(log(runif(1)) < log_accept_prob) {
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

#***********************************************************
#******** MCMC with prior *********** # 
#*#***********************************************************
mcmc_ss_prior <- function(data, n, sigma,  sigma_b, x0 = 1, prior = TRUE) { #burn_in = 2500
  
  'Returns mcmc samples of alpha & acceptance rate'
  
  #Set up
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  alpha_vec[1] <- x0; beta_vec[1] <- x0;
  gamma_vec[1] <- x0; r0_vec[1] <- x0;
  
  #Priors + Likelihoods
  like_a <- vector('numeric', n); like_b <- vector('numeric', n); 
  like_g = vector('numeric', n); prior_gamma <- vector('numeric', n); 
  prior_alpha <- vector('numeric', n); prior_beta <- vector('numeric', n); 
  
  count_accept1 = 0; count_accept2 = 0; count_accept3 = 0;
  
  #MCMC chain
  for(i in 2:n) {
    
    # #Print
    # if(mod(i, 10) == 0){
    #   print(i)
    # }
    # 
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
    #Store
    like_a[i] = log_like_ss_lse(data,  alpha_vec[i], beta_vec[i-1], gamma_vec[i-1])
    prior_alpha[i] = exp(-alpha_vec[i])
    
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
    #Store
    like_b[i] = log_like_ss_lse(data,  alpha_vec[i], beta_vec[i], gamma_vec[i-1])
    prior_beta[i] = exp(-beta_vec[i])
    
    #************************************************************************
    #gamma
    gamma_dash <- gamma_vec[i-1] + rnorm(1, sd = sigma) 
    
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
    #if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
    if(log(runif(1)) < log_accept_prob) {
      gamma_vec[i] <- gamma_dash
      count_accept3 = count_accept3 + 1
    } else {
      gamma_vec[i] <- gamma_vec[i-1]
    }
    #Priors for plotting
    like_g[i] = log_like_ss_lse(data,  alpha_vec[i], beta_vec[i], gamma_vec[i])
    prior_gamma[i] = 1 + exp(-gamma_vec[i])
    
    #R0
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
              accept_rate1, accept_rate2, accept_rate3,
              like_a, like_b, like_g,
              prior_alpha, prior_beta, prior_gamma))
}
#*****************************
#Plot results
plot_mcmc_results_x4 <- function(sim_data, mcmc_params, true_r0, dist_type, total_time, seed_count, prior){
  
  #Plot Set up
  plot.new()
  par(mfrow=c(3,4))
  
  #Extract params
  alpha_mcmc = mcmc_params[1]
  alpha_mcmc = unlist(alpha_mcmc)
  
  beta_mcmc = mcmc_params[2]
  beta_mcmc = unlist(beta_mcmc)
  
  gamma_mcmc = mcmc_params[3]
  gamma_mcmc = unlist(gamma_mcmc)
  
  r0_mcmc = mcmc_params[4]
  r0_mcmc = unlist(r0_mcmc)
  
  #Cumulative means + param sample limits
  #r0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  r0_lim = max(true_r0, max(r0_mcmc))
  r0_lim2 = max(true_r0, r0_mean)
  
  #alpha
  alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
  a_lim =  max(alphaX, max(alpha_mcmc))
  a_lim2 =  max(alphaX, alpha_mean)
  
  #beta
  beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
  b_lim = max(betaX, max(beta_mcmc))
  b_lim2 = max(betaX, beta_mean)
  
  #gamma
  gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
  g_lim =  max(gammaX, max(gamma_mcmc))
  g_lim2 =  max(gammaX, gamma_mean) 
  
  #***********
  #* Plots *
  
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(seed_count, "Day Infts SS Evnts", dist_type, "r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC Trace Plots
  plot.ts(alpha_mcmc, ylab = 'alpha', ylim=c(0, a_lim),
          main = paste("MCMC SS Events, true alpha = ", alphaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = alphaX, col = 'red', lwd = 2)
  
  plot.ts(beta_mcmc, ylab = 'beta', ylim=c(0, b_lim),
          main = paste("MCMC SS Events, true beta = ", betaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = betaX, col = 'blue', lwd = 2)
  
  plot.ts(gamma_mcmc,  ylab = 'gamma', ylim=c(0,g_lim),
          main = paste("MCMC SS Events, true gamma = ", gammaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = gammaX, col = 'chartreuse4', lwd = 2)
  #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
  
  #iii. Cumulative mean plots
  #r0 Mean
  plot2 = plot(seq_along(r0_mean), r0_mean,
               ylim=c(0, r0_lim),
               xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #alpha mean
  plot2 = plot(seq_along(alpha_mean), alpha_mean,
               ylim=c(0, a_lim),
               xlab = 'Time', ylab = 'alpha', main = paste("Alpha MCMC mean, True alpha = ",alphaX),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = alphaX, col = 'red', lwd = 2)
  
  #beta mean
  plot2 = plot(seq_along(beta_mean), beta_mean,
               ylim=c(0, b_lim),
               xlab = 'Time', ylab = 'beta', main = paste("Beta MCMC mean, True beta = ",betaX),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = betaX, col = 'blue', lwd = 2)
  
  #gamma Mean
  plot2 = plot(seq_along(gamma_mean), gamma_mean,
               xlab = 'Time', ylab = 'gamma', main = paste("Gamma MCMC mean, True gamma = ",gammaX),
               ylim=c(0, g_lim),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = gammaX, col = 'green', lwd = 2)
  
  #iv. Param Histograms (Plots 9,11,12)
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = paste('R0 total MCMC samples. Prior = ', prior),
       xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_r0, col = 'orange', lwd = 2)
  
  #v. Beta vs gamma
  plot(beta_mcmc, gamma_mcmc,
       xlab = 'beta', ylab = 'gamma', main = 'Beta vs Gamma',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Hist Beta 
  hist(beta_mcmc, freq = FALSE, breaks = 100,
       xlab = 'beta', #ylab = 'Density', 
       main = paste("Beta, True beta = ", betaX), 
       xlim=c(0, b_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = betaX, col = 'blue', lwd = 2)
  
  #Hist Gamma 
  hist(gamma_mcmc, freq = FALSE, breaks = 100,
       xlab = 'gamma', #ylab = 'Density', 
       main = paste("Gamma, True gamma = ", gammaX),
       xlim=c(0, g_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = gammaX, col = 'green', lwd = 2)
  
  #Final Mean Stats
  data_10_pc = 0.1*n
  a_mcmc_mean = round(mean(alpha_mcmc[n-data_10_pc:n]), 2)
  b_mcmc_mean = round(mean(beta_mcmc[n-data_10_pc:n]), 2)
  g_mcmc_mean = round(mean(gamma_mcmc[n-data_10_pc:n]), 2)
  r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
  
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

#PLOT WITH PRIORS 
plot_mcmc_x4_priors <- function(sim_data, mcmc_params, true_r0, dist_type, total_time, seed_count, prior){
    
    #Plot Set up
    plot.new()
    par(mfrow=c(3,4))
    
    #Extract params
    alpha_mcmc = mcmc_params[1]
    alpha_mcmc = unlist(alpha_mcmc)
    like_a = mcmc_params[8]
    like_a = unlist(like_a)
    prior_a = mcmc_params[11]
    prior_a = unlist(prior_a)
      
    beta_mcmc = mcmc_params[2]
    beta_mcmc = unlist(beta_mcmc)
    like_b = mcmc_params[9]
    like_b = unlist(like_b)
    prior_b = mcmc_params[12]
    prior_b = unlist(prior_b)
    
    gamma_mcmc = mcmc_params[3]
    gamma_mcmc = unlist(gamma_mcmc)
    like_g = mcmc_params[10]
    like_g = unlist(like_g)
    prior_g = mcmc_params[13]
    prior_g = unlist(prior_g)
    
    r0_mcmc = mcmc_params[4]
    r0_mcmc = unlist(r0_mcmc)
    
    #Cumulative means + param sample limits
    #r0
    r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
    r0_lim = max(true_r0, max(r0_mcmc))
    r0_lim2 = max(true_r0, r0_mean)
    
    #alpha
    alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
    a_lim =  max(alphaX, max(alpha_mcmc))
    a_lim2 =  max(alphaX, alpha_mean)
    
    #beta
    beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
    b_lim = max(betaX, max(beta_mcmc))
    b_lim2 = max(betaX, beta_mean)
    
    #gamma
    gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
    g_lim =  max(gammaX, max(gamma_mcmc))
    g_lim2 =  max(gammaX, gamma_mean) 
    
    #***********
    #* Plots *
    
    #i.Infections
    plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
            main = paste(seed_count, "Day Infts SS Evnts", dist_type, "r0 = ", true_r0),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #ii. MCMC Trace Plots
    plot.ts(alpha_mcmc, ylab = 'alpha', ylim=c(0, a_lim),
            main = paste("SS Evnts, true a = ", alphaX, '(green), ll(red), prior(blue)'),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(h = alphaX, col = 'green', lwd = 2) #True = green
    lines(seq_along(like_a), like_a, col = 'red')
    lines(seq_along(prior_a), prior_a, col = 'blue')
    
    plot.ts(beta_mcmc, ylab = 'beta', ylim=c(0, b_lim),
            main = paste("MCMC SS Events, true beta = ", betaX),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(h = betaX, col = 'green', lwd = 2) #True = green
    lines(seq_along(like_b), like_b, col = 'red')
    lines(seq_along(prior_b), prior_b, col = 'blue')
    
    plot.ts(gamma_mcmc,  ylab = 'gamma', ylim=c(0,g_lim),
            main = paste("MCMC SS Events, true gamma = ", gammaX),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(h = gammaX, col = 'green', lwd = 2) #True = green
    lines(seq_along(like_g), like_g, col = 'red')
    lines(seq_along(prior_g), prior_g, col = 'blue')
    
    #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
    
    #iii. Cumulative mean plots
    #r0 Mean
    plot2 = plot(seq_along(r0_mean), r0_mean,
                 ylim=c(0, r0_lim),
                 xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
                 cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    print(plot2)
    abline(h = true_r0, col = 'orange', lwd = 2)
    
    #alpha mean
    plot2 = plot(seq_along(alpha_mean), alpha_mean,
                 ylim=c(0, a_lim),
                 xlab = 'Time', ylab = 'alpha', main = paste("Alpha MCMC mean, True alpha = ",alphaX),
                 cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    print(plot2)
    abline(h = alphaX, col = 'red', lwd = 2)
    
    #beta mean
    plot2 = plot(seq_along(beta_mean), beta_mean,
                 ylim=c(0, b_lim),
                 xlab = 'Time', ylab = 'beta', main = paste("Beta MCMC mean, True beta = ",betaX),
                 cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    print(plot2)
    abline(h = betaX, col = 'blue', lwd = 2)
    
    #gamma Mean
    plot2 = plot(seq_along(gamma_mean), gamma_mean,
                 xlab = 'Time', ylab = 'gamma', main = paste("Gamma MCMC mean, True gamma = ",gammaX),
                 ylim=c(0, g_lim),
                 cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    print(plot2)
    abline(h = gammaX, col = 'green', lwd = 2)
    
    #iv. Param Histograms (Plots 9,11,12)
    hist(r0_mcmc, freq = FALSE, breaks = 100,
         xlab = 'R0 total', #ylab = 'Density', 
         main = paste('R0 total MCMC samples. Prior = ', prior),
         xlim=c(0, r0_lim),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = true_r0, col = 'orange', lwd = 2)
    
    #v. Beta vs gamma
    plot(beta_mcmc, gamma_mcmc,
         xlab = 'beta', ylab = 'gamma', main = 'Beta vs Gamma',
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #Hist Beta 
    hist(beta_mcmc, freq = FALSE, breaks = 100,
         xlab = 'beta', #ylab = 'Density', 
         main = paste("Beta, True beta = ", betaX), 
         xlim=c(0, b_lim),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = betaX, col = 'blue', lwd = 2)
    
    #Hist Gamma 
    hist(gamma_mcmc, freq = FALSE, breaks = 100,
         xlab = 'gamma', #ylab = 'Density', 
         main = paste("Gamma, True gamma = ", gammaX),
         xlim=c(0, g_lim),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = gammaX, col = 'green', lwd = 2)
    
    #Final Mean Stats
    data_10_pc = 0.1*n
    a_mcmc_mean = round(mean(alpha_mcmc[n-data_10_pc:n]), 2)
    b_mcmc_mean = round(mean(beta_mcmc[n-data_10_pc:n]), 2)
    g_mcmc_mean = round(mean(gamma_mcmc[n-data_10_pc:n]), 2)
    r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
    
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
alphaX = 0.8 #1.1 #0.8 #1.1 # 0.8 #2 #0.9 #2 #2 #Without ss event, ~r0.
betaX = 0.2 #0.05 #0.2 #0.05 #0.2 #0.2 #0.05 #0.2 #0.05 #0.05
gammaX = 10
true_r0 = alphaX + betaX*gammaX
true_r0
##!!---##############################################################!!---##

#Data: Set seed 
# for (i in 1:seed_count) {
#   cat('Seed count; ', i)
#   set.seed(i)
# }
set.seed(seed_count)

#Epidemic data - Neg Bin
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')
sim_data = sim_data

#Epidemic data - Poisson 
#sim_data = simulate_ss_poisson(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
#plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')

#MCMC 
n = 30000
start_time = Sys.time()
print('Start time:')
print(start_time)
sigma = 1
sigma_b = 0.05
prior = FALSE
#mcmc_params = mcmc_super_spreading(sim_data, n, sigma, sigma_b, x0 = 1)
prior = TRUE
mcmc_params = mcmc_ss_prior(sim_data, n, sigma, sigma_b, x0 = 1)
end_time = Sys.time()
time_elap = round(end_time - start_time, 2)
print('Time elapsed:')
print(time_elap)

#Plotting 
dist_type = 'Neg Bin,'
#plot_mcmc_results_total(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count)
plot_mcmc_results_x4(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count, prior)


#Seed
seed_count = seed_count + 1
seed_count
