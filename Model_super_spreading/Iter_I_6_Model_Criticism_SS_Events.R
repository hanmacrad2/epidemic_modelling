#Model Criticism
#Super spreading Events model- MCMC Inference + Model Criticism
#Setup
library(MASS)
library(pracma)
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
  vec_data_simulated <- vector('numeric', n/thinning_factor)
  count_thin = 1
  
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
    
    #Thinning Factor
    #print(mod(i, thinning_factor))
    if(mod(i, thinning_factor) == 0){
      #cat('i =', i)
      vec_data_simulated[count_thin] = sum(simulate_branching_ss(num_days, shape_gamma, scale_gamma, alpha_vec[i], beta_vec[i], gamma_vec[i]))
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
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3,
              accept_rate4,
              vec_data_simulated))
}


#Model Criticism Function
model_criticism <- function(mcmc_params, sim_data, max_sum_val) { 
  
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



#****************
#MCMC Plots 4x4
plot_mcmc_x4_II <- function(sim_data, mcmc_params, true_r0, dist_type,
                                total_time, seed_count, prior, max_sum_val, model_crit = TRUE, joint = TRUE){
  
  #Plot Set up
  plot.new()
  
  if (joint){
    par(mfrow=c(4,4))
  } else {
    par(mfrow=c(3,4))
  }
  
  
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
  abline(h = alphaX, col = 'red', lwd = 2) #True = green
  #lines(seq_along(like_a), like_a, col = 'red')
  #lines(seq_along(prior_a), prior_a, col = 'blue')
  
  plot.ts(beta_mcmc, ylab = 'beta', ylim=c(0, b_lim),
          main = paste("MCMC SS Events, true beta = ", betaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = betaX, col = 'blue', lwd = 2) #True = green
  
  plot.ts(gamma_mcmc,  ylab = 'gamma', ylim=c(0,g_lim),
          main = paste("MCMC SS Events, true gamma = ", gammaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = gammaX, col = 'green', lwd = 2) #True = green
  
  #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
  
  #Title
  text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
       line2user(line=2, side=3), paste('MCMC SS, True R0:', true_r0, 'Prior = ', prior), xpd=NA, cex=2, font=2)
  
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
  #lines(seq_along(like_a_mean), like_a_mean, col = 'red')
  #lines(seq_along(prior_a_mean), prior_a_mean, col = 'blue')
  
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
  
  # #v. Beta vs gamma
  # plot(beta_mcmc, gamma_mcmc,
  #      xlab = 'beta', ylab = 'gamma', main = 'Beta vs Gamma',
  #      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Hist alpha 
  hist(alpha_mcmc, freq = FALSE, breaks = 100,
       xlab = 'alpha', #ylab = 'Density', 
       main = paste("alpha, True alpha = ", alphaX), 
       xlim=c(0, a_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = alphaX, col = 'red', lwd = 2)
  #Prior
  # x <- seq(from = 0, to = 20, by = 0.05)
  # exp1 = dexp(x, 1)
  # lines(seq_along(exp1), exp1, type = 'l')

  
  #Hist Beta 
  hist(beta_mcmc, freq = FALSE, breaks = 100,
       xlab = 'beta', #ylab = 'Density', 
       main = paste("Beta, True beta = ", betaX), 
       xlim=c(0, b_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = betaX, col = 'blue', lwd = 2)
  #Prior
  # x <- seq(from = 0, to = 20, by = 0.05)
  # exp1 = dexp(x, 1)
  # lines(seq_along(exp1), exp1, type = 'l')
  
  
  #Hist Gamma 
  hist(gamma_mcmc, freq = FALSE, breaks = 100,
       xlab = 'gamma', #ylab = 'Density', 
       main = paste("Gamma, True gamma = ", gammaX),
       xlim=c(0, g_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = gammaX, col = 'green', lwd = 2)
  #Prior
  # x <- seq(from = 0, to = 20, by = 0.05)
  # exp1 = 1 + dexp(x, 1)
  # lines(seq_along(exp1), exp1, type = 'l')
  
  #Final Mean Stats
  data_10_pc = 0.5*n #50%
  a_mcmc_mean = round(mean(alpha_mcmc[n-data_10_pc:n]), 2)
  b_mcmc_mean = round(mean(beta_mcmc[n-data_10_pc:n]), 2)
  g_mcmc_mean = round(mean(gamma_mcmc[n-data_10_pc:n]), 2)
  r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
  
  #Joint distrbutions
  if (model_crit){
    
    model_criticism(mcmc_params, sim_data, max_sum_val)
  }
  
  #Joint distrbutions
  if (joint){
    
    # #v. r0 vs beta
    # plot(beta_mcmc, r0_mcmc,
    #      xlab = 'beta', ylab = 'R0', main = 'Beta vs R0',
    #      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. alpha vs beta
    plot(alpha_mcmc, beta_mcmc,
         xlab = 'alpha', ylab = 'beta', main = 'alpha vs Beta',
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. alpha vs gamma
    plot(alpha_mcmc, gamma_mcmc,
         xlab = 'alpha', ylab = 'gamma', main = 'alpha vs gamma',
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. beta vs gamma
    plot(beta_mcmc, gamma_mcmc,
         xlab = 'beta', ylab = 'gamma', main = 'Beta vs gamma',
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }
  
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
    a_rte_b_g = round(mcmc_params[[8]],2),
    tot_time = total_time) 
  
  print(df_results)
  
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
mcmc_params = mcmc_ss_mod_crit(sim_data, n, sigma_a, sigma_b, sigma_g, sigma_bg, prior)
end_time = Sys.time()
time_elap = round(end_time - start_time, 2)
print('Time elapsed:')
print(time_elap)

#Apply
dist_type = 'Neg Bin,'
max_sum_val = 5000
#plot_mcmc_x4_priors(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count, prior)
plot_mcmc_x4_II(sim_data, mcmc_params, true_r0, dist_type, time_elap, seed_count, prior, max_sum_val)

#Model Criticism
par(mfrow = c(1,1))
model_criticism(mcmc_params, sim_data, max_sum_val)

par(mfrow = c(2,1))
