#Epidemic Modelling Functions
#-Simulations fucntions
#- Likelihood Functions

#*********************************************************
#*Simulation Functions

#Baseline simulation
simulate_branching = function(num_days, r0, shape_gamma, scale_gamma) {
  
  'Baseline simulation model'
  #Set up
  vec_infecteds = vector('numeric', num_days)
  vec_infecteds[1] = 2
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Total rate
    tot_rate = r0*sum(vec_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    vec_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
  }
  
  vec_infecteds
}

#*******************************************************
#Super-spreading simulation
simulate_branching_ss = function(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX) {
  'Simulate an epidemic with Superspreading events
  prop_ss = Proportion of superspreading days
  magnitude_ss = increased rate of superspreading event'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nsse_infecteds = vector('numeric', num_days)
  sse_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 2
  nsse_infecteds[1] = 2
  sse_infecteds[1] = 0
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(total_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Eg. Today is day 10. Infected on day 1. Current Infectiousness is t - day_infected 
    tot_rate = alphaX*lambda_t #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nsse_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    sse_infecteds[t] = rnbinom(1, betaX*lambda_t, 1/(1 + gammaX)) #z_t: Total infecteds due to super-spreading event - num of events x Num individuals
    
    #Total 
    total_infecteds[t] = nsse_infecteds[t] + sse_infecteds[t]
  }
  
  total_infecteds
}
# 
# #********
# # #*Implement
# num_days = 50
# #lambda params
# shape_gamma = 6
# scale_gamma = 1
# #params
# alphaX = 0.8 #Without ss event, ~r0.
# betaX = 0.1
# gammaX = 10
# true_r0 = alphaX + betaX*gammaX
# true_r0
# #Epidemic data
# sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
# plot.ts(sim_data, ylab = 'Daily Infections', main = paste('Super - Spreading Events model - Daily Infections count, true R0 = ', true_r0))
# par(mfrow = c(2,1))

#*******************************************************
#Super-spreading simulation
simulate_ss_poisson = function(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX) {
  'Simulate an epidemic with Superspreading events
  prop_ss = Proportion of superspreading days
  magnitude_ss = increased rate of superspreading event'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nsse_infecteds = vector('numeric', num_days)
  sse_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 2
  nsse_infecteds[1] = 2
  sse_infecteds[1] = 0
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum(total_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    tot_rate = alphaX*lambda_t #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nsse_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, betaX*lambda_t) #Number of super-spreading events (beta)
    sse_infecteds[t] = rpois(1, gammaX*n_t) #z_t: Total infecteds due to super-spreading event - num of events x Num individuals
    
    total_infecteds[t] = nsse_infecteds[t] + sse_infecteds[t]
  }
  
  total_infecteds
}

#*Implement
# num_days = 50
# #lambda params
# shape_gamma = 6
# scale_gamma = 1
# #params
# alphaX = 1.2 #Without ss event, ~r0. 
# betaX = 0.5
# gammaX = 10
# #Epidemic data
# sim_data2 = simulate_ss_poisson(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
# plot.ts(sim_data2, ylab = 'Daily Infections count', main = 'Daily Infections count')


#*******************************************************
#Super-spreaders simulation
simulation_super_spreaders = function(num_days, shape_gamma, scale_gamma, aX, bX, ss_mult) {
  'Simulate an epidemic with Superspreading events
  prop_ss = Proportion of superspreading days
  magnitude_ss = increased rate of superspreading event'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nss_infecteds = vector('numeric', num_days)
  ss_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 2
  nss_infecteds[1] = 2
  ss_infecteds[1] = 0 
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((total_infecteds[1:(t-1)] + ss_mult*ss_infecteds[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nss_infecteds[t] = rpois(1, aX*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss_infecteds[t] = rpois(1, bX*lambda_t)
    total_infecteds[t] = nss_infecteds[t] + ss_infecteds[t]
  }
  
  total_infecteds
}

#*Implement
num_days = 50
#lambda params
shape_gamma = 6
scale_gamma = 1
#params
aX = 0.8 #1.1 #Without ss event, ~r0.
bX = 0.1 #0.2
ss_mult = 10 #8
#Epidemic data
sim_data2 = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, aX, bX, ss_mult)
plot.ts(sim_data2, ylab = 'Daily Infections count', main = 'Super Spreaders Model - Daily Infections count')

###################################################################################
# MODELS!

###############################
# BASELINE

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

#*******************************************************************************
#SUPER-SPREADING MODEL

#Log Likelihood 
log_like_ss <- function(x, alphaX, betaX, gammaX){
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma) i,e Prob less than x2 - prob less than x1; the area in between 
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    inner_sum_xt = 0
    
    for (y_t in 0:x[t]){ #Sum for all values of y_t
      
      #Log likelihood
      inner_sum_xt = (inner_sum_xt + exp(-alphaX*lambda_t)*(1/factorial(y_t))*(alphaX*lambda_t)^y_t*
                        (gamma((x[t] - y_t) + betaX*lambda_t))/(gamma(betaX*lambda_t)*
                                                                  factorial(x[t] - y_t))*(1/(gammaX +1))^(betaX*lambda_t)*
                        (gammaX/(gammaX + 1))^(x[t] - y_t))
      
    } 
    
    logl = logl + log(inner_sum_xt) 
    
  }
  
  logl
  
}

#*****************************************************
#Log Likelihood - log-exp-sum trick 
log_like_ss_lse <- function(x, alphaX, betaX, gammaX){
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)),
                                                                                         shape = shape_gamma, scale = scale_gamma)
  logl = 0 
  
  for (t in 2:num_days) {
    
    #print(t)
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    
    if(x[t] == 0){ #y_t also equal to zero
      
      #L(x_i) for y_t, x_t = 0
      logl = logl -(alphaX*lambda_t) - 
        (betaX*lambda_t*log(gammaX +1))
      
    } else {
      
      #Terms in inner sum
      inner_sum_vec <- vector('numeric', x[t])
      
      for (y_t in 0:x[t]){ #Sum for all values of y_t up to x_t
        
        #Store inner L(x_i) term in vector position
        inner_sum_vec[y_t + 1] = (-(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) +
                                    lgamma((x[t] - y_t) + (betaX*lambda_t)) - lgamma(betaX*lambda_t) - 
                                    lfactorial(x[t] - y_t) - (betaX*lambda_t*log(gammaX +1)) + 
                                    (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1))
        
      }
      
      #Calculate max element in inner vector, for all y_t for a given t, x[t]
      lx_max = max(inner_sum_vec)
      
      #Calculate lse
      lse = lx_max + log(sum(exp(inner_sum_vec - lx_max) ))
      
      #Add to overall log likelihood 
      logl = logl + lse 
      
    }
    
  }
  
  logl
  
}

#*********************************************************************************************
#* Log Likelihood - Using package distribution functions
log_like_ss_R_funcs <- function(x, alphaX, betaX, gammaX){
  
  'Calculate the log likelihood on the simulated data using package distribution functions in R'
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  logl = 0
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Over all timepoints
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:(t-1)]*rev(prob_infect[1:(t-1)]))
    inner_sum = 0
    
    #Sum for all values of y_t
    for (y_t in 0:x[t]){ 
      
      #Log likelihood
      prob_yt = dpois(y_t, alphaX*lambda_t)
      #cat("prob_yt, dpois : ", prob_yt, "\n")
      
      prob_zt = s(x[t] - y_t, betaX*lambda_t, 1/(1 + gammaX))
      #cat("prob_zt, dnbinom : ", prob_zt, "\n")
      
      prob_xt = prob_zt*prob_yt
      inner_sum = inner_sum + prob_xt
      #cat("prob_xt, : ", prob_xt, "\n")
      
    } 
    
    #cat("log(prob_xt), : ", log(prob_xt), "\n")
    logl = logl + log(inner_sum) 
    #cat("logl, : ", logl, "\n")
    
  }
  
  logl
}

################################################################################
# Plot MCMC Total Results
################################################################################

plot_mcmc_results <- function(sim_data, mcmc_params, true_r0, 
                              dist_type, total_time, seed_count, prior, beta_gamma_flag){
  
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
  
  #Beta-gamma
  if (beta_gamma_flag){
    b_g_a_rate = mcmc_params[14]
    b_g_a_rate = round(unlist(b_g_a_rate), 2)
  } else {
    b_g_a_rate = 0
  }
  
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
  abline(h = alphaX, col = 'green', lwd = 2)
  
  plot.ts(beta_mcmc, ylab = 'beta', ylim=c(0, b_lim),
          main = paste("MCMC SS Events, true beta = ", betaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = betaX, col = 'green', lwd = 2)
  
  plot.ts(gamma_mcmc,  ylab = 'gamma', ylim=c(0,g_lim),
          main = paste("MCMC SS Events, true gamma = ", gammaX),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = gammaX, col = 'green', lwd = 2)
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
    b_g_a_rate = b_g_a_rate, 
    tot_time = total_time) 
  
  print(df_results)
  
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

###############################################################################
#FUNCTION TO PLOT 1x4 DASHBOARD OF MCMC RESULTS FOR BASE MODEL
plot_mcmc_results_r0 <- function(sim_data, mcmc_params, true_r0, time_elap, seed_count){
  
  #Plot Set up
  #par(mar=c(1,1,1,1))
  plot.new()
  par(mfrow=c(2,2))
  
  #Extract r0
  r0_mcmc = mcmc_params[1]
  r0_mcmc = unlist(r0_mcmc)
  
  #Cumulative means + param sample limits
  #r0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  r0_lim = max(true_r0, max(r0_mcmc))
  r0_lim2 = max(true_r0, r0_mean)
  
  
  #***********
  #* Plots *
  
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(seed_count, "Day Infts SS Evnts, ", "r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC Trace Plots
  plot.ts(r0_mcmc, ylab = 'alpha', ylim=c(0, a_lim),
          main = paste("MCMC SS Events, true r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = true_r0, col = 'orange', lwd = 2) #True = green
  
  #iii. Cumulative mean plots
  #r0 Mean
  plot2 = plot(seq_along(r0_mean), r0_mean,
               ylim=c(0, r0_lim),
               xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #iv. Histogram
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = paste('R0 total MCMC samples. Prior = ', prior),
       xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_r0, col = 'orange', lwd = 2)
  
  #Final Mean Stats
  data_10_pc = 0.5*n #50%
  r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
  
  #Results
  df_results <- data.frame(
    R0 = true_r0, 
    R0_mc = r0_mcmc_mean,
    accept_rate_r0 = round(mcmc_params[[2]],2),
    time_elap = round(time_elap,2)) 
  
  print(df_results)
  
}

###############################################################################
#FUNCTION TO PLOT 4x4 DASHBOARD OF MCMC RESULTS FOR SUPER SPREADING EVENTS MODEL
plot_mcmc_x4_priors <- function(sim_data, mcmc_params, true_r0, dist_type, total_time, seed_count, prior, joint = TRUE){
  
  #Plot Set up
  #par(mar=c(1,1,1,1))
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
  # text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
  #      line2user(line=2, side=3), paste('MCMC SS, True R0:', true_r0, 'Prior = ', prior), xpd=NA, cex=2, font=2)
  # 
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
  #Prior
  #x <- seq(from = 0, to = 10, by = 0.5)
  #exp1 = dexp(x, 1)
  #lines(exp1, col = 'purple')
  abline(v = alphaX, col = 'red', lwd = 2)
  
  
  #Hist Beta 
  hist(beta_mcmc, freq = FALSE, breaks = 100,
       xlab = 'beta', #ylab = 'Density', 
       main = paste("Beta, True beta = ", betaX), 
       xlim=c(0, b_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #Prior
  #lines(exp1, col = 'purple')
  abline(v = betaX, col = 'blue', lwd = 2)
  
  
  #Hist Gamma 
  hist(gamma_mcmc, freq = FALSE, breaks = 100,
       xlab = 'gamma', #ylab = 'Density', 
       main = paste("Gamma, True gamma = ", gammaX),
       xlim=c(0, g_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #Prior
  #x <- seq(from = 0, to = 20, by = 0.5)
  #exp1 = 1 + dexp(x, 1)
  #lines(exp1, col = 'purple')
  abline(v = gammaX, col = 'green', lwd = 2)
  
  
  #Final Mean Stats
  data_10_pc = 0.5*n #50%
  a_mcmc_mean = round(mean(alpha_mcmc[n-data_10_pc:n]), 2)
  b_mcmc_mean = round(mean(beta_mcmc[n-data_10_pc:n]), 2)
  g_mcmc_mean = round(mean(gamma_mcmc[n-data_10_pc:n]), 2)
  r0_mcmc_mean = round(mean(r0_mcmc[n-data_10_pc:n]), 2)
  
  #Joint distrbutions
  if (joint){
    
    #v. r0 vs beta
    plot(beta_mcmc, r0_mcmc,
         xlab = 'beta', ylab = 'R0', main = 'Beta vs R0',
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
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
