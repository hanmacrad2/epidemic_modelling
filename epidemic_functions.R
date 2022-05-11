#Epidemic Modelling Functions

#Contains:
#SIMULATION FUNCTIONS
#MODEL LIKELIHOOD FUNCTIONS
#PLOTTING MCMC FUNCTIONS

################################################################################
# SIMULATION
################################################################################

#*********************************************************
#I. BASELINE SIMULATION
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

#Apply
num_days = 50
#lambda params
shape_gamma = 6; scale_gamma = 1
r0 = 1.8
sim_data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
plot.ts(sim_data, ylab = 'Daily infection count',
        main = paste('Base Model - Daily infections count. R0 = ', r0))


#*******************************************************
#II. SUPER-SPREADING EVENTS SIMULATION
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
num_days = 50
#lambda params
shape_gamma = 6; scale_gamma = 1
#Params
alphaX = 0.8 #Without ss event, ~r0.
betaX = 0.1; gammaX = 10
true_r0 = alphaX + betaX*gammaX
true_r0
#Epidemic data
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily infection count',
        main = paste('Super-Spreading Events Model - Daily Infections. R0 = ', true_r0)) #,
                     #expression(alpha), ':', alphaX, 'b:', betaX, ':', gammaX))
# par(mfrow = c(2,1))

#*******************************************************
#II. SUPER-SPREADING EVENTS SIMULATION + POISSON
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
#III. SUPER-SRPEADING INDIVIDUALS (SUPER-SPREADERS) SIMULATION 
simulation_super_spreaders_v0 = function(num_days, shape_gamma, scale_gamma, aX, bX, cX) {
  'Simulate an epidemic with Superspreading individuals'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nss_infecteds = vector('numeric', num_days)
  ss_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 3
  nss_infecteds[1] = 2
  ss_infecteds[1] = 1 
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((total_infecteds[1:(t-1)] + cX*ss_infecteds[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nss_infecteds[t] = rpois(1, aX*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss_infecteds[t] = rpois(1, bX*lambda_t)
    total_infecteds[t] = nss_infecteds[t] + ss_infecteds[t]
  }
  
  total_infecteds
}

#*******************************************************
#III B. SUPER-SRPEADING INDIVIDUALS (SUPER-SPREADERS) SIMULATION 
simulation_super_spreaders = function(num_days, shape_gamma, scale_gamma, aX, bX, cX) {
  'Simulate an epidemic with Superspreading individuals'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nss_infecteds = vector('numeric', num_days)
  ss_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 3
  nss_infecteds[1] = 2
  ss_infecteds[1] = 1 
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((total_infecteds[1:(t-1)] + cX*ss_infecteds[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nss_infecteds[t] = rpois(1, aX*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss_infecteds[t] = rpois(1, bX*lambda_t)
    total_infecteds[t] = nss_infecteds[t] + ss_infecteds[t]
  }
  
  return(list(nss_infecteds, ss_infecteds))
}

#Return sum of n + s
get_total_x = function(data_list){
  
  'Returns sum of x = n + s in a single vector'
  n = data[[1]]; s = data[[2]]
  x = n + s
  x
} 

#*Implement
num_days = 50
#lambda params
shape_gamma = 6; scale_gamma = 1
#params
aX = 0.8 #1.1 #Without ss event, ~r0.
bX = 0.1; cX = 10 #8
r02 = aX + bX*cX
#Epidemic data
sim_data2 = simulation_super_spreaders_v0(num_days, shape_gamma, scale_gamma, aX, bX, cX)
plot.ts(sim_data2, ylab = 'Daily infection count',
        main = paste('Super Spreaders Model - Daily Infections. R0 = ', r02))

##############################################################################
# MODELS
##############################################################################

#*******************************************************************************
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
    logl = logl + y[t]*log(lambda) - lambda - lfactorial(y[t]) #Need to include normalizing constant 
    
  }
  
  logl
  
}

# BASE MCMC
mcmc_r0 <- function(data, n, sigma, burn_in, x0 = 1) {
  
  'Returns mcmc samples of R0. 
  Prior on R_0 (or alpha) = exp(1)'
  #Data
  print('sim data')
  print(data)
  
  #Set up
  r0_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  U <- runif(n)
  count_accept = 0
  
  #MCMC chain
  for(i in 2:n) {
    r0_dash <- r0_vec[i-1] + rnorm(1, sd = sigma) #, mean = 0, sd = sigma_opt)
    if(r0_dash < 0){
      r0_dash = abs(r0_dash)
    }
    
    #Alpha
    log_alpha = log_like(data, r0_dash) - log_like(data, r0_vec[i-1]) - r0_dash + r0_vec[i-1] #exponential prior
    #log_alpha = log_like(data, Y) - log_like(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) 
    #Should include: Likelihood + prior + propogsal density x2 (Previous time step & Current time step)
    
    if (is.na(log_alpha)){
      print('na value')
      sprintf("r0_dash: %i", r0_dash)
    }
    if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) { 
      r0_vec[i] <- r0_dash
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
    }
  }
  #Final stats
  accept_rate = 100*(count_accept/n)
  print(paste0("Acceptance rate = ",accept_rate))
  
  r0_vec = r0_vec[burn_in:n]
  r0_vec
  
  return(list(r0_vec, accept_rate))
}

#*******************************************************************************
# SUPER-SPREADING EVENTS MODEL
#*******************************************************************************

#Log Likelihood 
log_like_sse <- function(x, alphaX, betaX, gammaX){
  
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
    
    if((x[t] == 0) | is.na(x[t])) { #y_t also equal to zero 
      #L(x_i) for y_t, x_t = 0
      logl = logl -(alphaX*lambda_t) - 
        (betaX*lambda_t*log(gammaX +1))
      
      #print(paste0('logl 1 = ', logl))
      
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
      #print(paste0('lx_max = ', lx_max))
      
      #Calculate lse
      lse = lx_max + log(sum(exp(inner_sum_vec - lx_max) ))
      
      #print(paste0('lse = ', lse))
      
      #Add to overall log likelihood 
      logl = logl + lse 
      
    }
    
  }
  
  logl
  
}

#***************************************************************#####
#GAMMA PRIOR ON BETA
log_gamma_dist <- function(param, gamma_priors){
  
  #Params
  shapeX = gamma_priors[1]
  scaleX = gamma_priors[2]
  
  log_gamma_dist = (1/lgamma(shapeX)*shapeX*log(scaleX))*(shapeX - 1)*log(param)*(-param/shapeX)
  
  log_gamma_dist
}

##############################
#1. MCMC
MCMC_SSE <- function(data, n, sigma, model_params, gamma_prior, gamma_priors,
                            x0 = 1, prior = TRUE, alpha_transform = FALSE) {#thinning_factor, burn_in
  
  'Returns MCMC samples of SSE model parameters (alpha, beta, gamma, r0 = a + b*g) 
  w/ acceptance rates. Includes alpha transform, beta-gamma transform' 
  print('MCMC SUPERSPREADING EVENTS')
  
  'Priors
  p(alpha) = exp(1) = rate*exp(-rate*x) = 1*exp(-1*alpha) = exp(-alpha). log(exp(-alpha)) = - alpha
  p(beta) = exp(1) or p(beta) = gamma(shape, scale), for e.g gamma(3, 2)
  p(gamma) = exp(1) + 1 = 1 + exp(-gamma) = exp(gamma - 1)'
  
  #Initialise params
  alpha_vec <- vector('numeric', n); beta_vec <- vector('numeric', n)
  gamma_vec <- vector('numeric', n); r0_vec <- vector('numeric', n)
  log_like_vec <- vector('numeric', n)
  
  alpha_vec[1] <- model_params[1]; beta_vec[1] <- model_params[2] 
  gamma_vec[1] <- model_params[3]; r0_vec[1] <- model_params[4];
  log_like_vec[1] <- log_like_ss_lse(data, alpha_vec[1], beta_vec[1],  gamma_vec[1])   
  
  alpha = alpha_vec[1]; beta =  beta_vec[1]; 
  gamma = gamma_vec[1]; log_like = log_like_vec[1]
  
  #Extract params
  sigma_a = sigma[1]; sigma_b = sigma[2]
  sigma_g = sigma[3]; sigma_bg = sigma[4];
  
  #Result vectors
  count_accept1 = 0; 
  count_accept2 = 0; count_reject2 = 0;
  count_accept3 = 0; count_reject3 = 0;
  count_accept4 = 0; count_reject4 = 0;
  
  #MCMC chain
  for(i in 2:n) {
    
    #******************************************************
    #ALPHA
    alpha_dash <- alpha + rnorm(1, sd = sigma_a) 
    if(alpha_dash < 0){
      alpha_dash = abs(alpha_dash)
    }
    
    #log alpha
    logl_new = log_like_ss_lse(data, alpha_dash, beta, gamma)
    log_accept_prob = logl_new - log_like  #+ prior1 - prior
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - alpha_dash + alpha
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      alpha <- alpha_dash
      count_accept1 = count_accept1 + 1
      log_like = logl_new
    } 

    #************************************************************************
    #BETA 
    #Only if (Beta > 0){ WHY? Took it out
    beta_dash <- beta + rnorm(1, sd = sigma_b) 
    if(beta_dash < 0){
      beta_dash = abs(beta_dash)
    }
    #loglikelihood
    logl_new = log_like_ss_lse(data, alpha, beta_dash, gamma)
    log_accept_prob = logl_new - log_like #logl_prev
    
    #Priors
    if (gamma_prior){
      log_accept_prob = log_accept_prob + log_gamma_dist(beta_dash, gamma_priors) - log_gamma_dist(beta, gamma_priors) 
    } else {
      log_accept_prob = log_accept_prob - beta_dash + beta 
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      beta <- beta_dash
      log_like = logl_new
      count_accept2 = count_accept2 + 1
    } else {
      count_reject2 = count_reject2 + 1
    }
    
    #************************************************************************
    #GAMMA
    gamma_dash <- gamma + rnorm(1, sd = sigma_g) 
    if(gamma_dash < 1){
      gamma_dash = 2 - gamma_dash #Prior on gamma - gt 1
    }
    #Acceptance Probability
    logl_new = log_like_ss_lse(data, alpha, beta, gamma_dash)
    log_accept_prob = logl_new - log_like #logl_prev 
    #Priors
    if (prior){
      log_accept_prob = log_accept_prob - gamma_dash + gamma
    }
    
    #Metropolis Acceptance Step
    if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
      gamma <- gamma_dash
      log_like <- logl_new
      count_accept3 = count_accept3 + 1
    } else {
      count_reject3 = count_reject3 + 1
    }
    
    #*****************************************************
    #GAMMA-BETA
    gamma_dash <- gamma + rnorm(1, sd = sigma_bg) #Alter sigma_bg depending on acceptance rate.
    #Acc rate too big -> Make sigma bigger. Acc rate too small -> make sigma smaller
    if(gamma_dash < 1){ #If less then 1
      gamma_dash = 2 - gamma_dash #abs(gamma_dash)
    }
    #New Beta
    r0 = alpha + beta*gamma 
    beta_new = (r0 - alpha)/gamma_dash #Proposing new Gamma AND Beta. Beta_dash = f(R0 & gamma_dash)
    
    if(beta_new >= 0){ #Only accept values of beta > 0
      
      logl_new = log_like_ss_lse(data, alpha, beta_new, gamma_dash)
      log_accept_prob = logl_new - log_like 
      
      #Priors: Gamma or Exp
      if (gamma_prior){
        log_accept_prob = log_accept_prob + log_gamma_dist(beta_new, gamma_priors) - log_gamma_dist(beta, gamma_priors)
      } else {  
        log_accept_prob = log_accept_prob - beta_new + beta
      }
      
      #Metropolis Step
      if(!(is.na(log_accept_prob)) && log(runif(1)) < log_accept_prob) {
        beta <- beta_new
        gamma <- gamma_dash
        count_accept4 = count_accept4 + 1
      } else {
        count_reject4 = count_reject4 + 1
      }
    }
    
    #Populate Model parameters w/ current values
    alpha_vec[i] <- alpha; beta_vec[i] <- beta
    gamma_vec[i] <- gamma; r0_vec[i] <- alpha + beta*gamma
    log_like_vec[i] <- log_like
    
  }
  
  #Bayes Factor
  #beta_pc0 = (length(which(beta_vec == 0)))/length(beta_vec) #Check beta_mcmc
  #bayes_factor = beta_pc0/(1-beta_pc0); bayes_factor = round(bayes_factor, 6)
  
  #Final stats
  accept_rate1 = 100*count_accept1/(n-1)
  accept_rate2 = 100*count_accept2/(count_accept2 + count_reject2)
  accept_rate3 = 100*count_accept3/(count_accept3 + count_reject3)
  accept_rate4 = 100*count_accept4/(count_accept4 + count_reject4)
  
  #Return alpha, acceptance rate
  return(list(alpha_vec, beta_vec, gamma_vec, r0_vec,
              accept_rate1, accept_rate2, accept_rate3, accept_rate4,
              count_accept2, count_accept3, count_accept4))
}


################################################################################
# PLOTTING
################################################################################

#*******************************************************************************
#R0: FUNCTION TO PLOT 1x4 DASHBOARD OF MCMC RESULTS FOR BASE MODEL
plot_mcmc_results_r0 <- function(n, sim_data, mcmc_params, true_r0, time_elap, seed_count, model_type){
  
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
  
  #***********
  #* Plots *
  
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(seed_count, "Orig Infts,", model_type, ", R0 = ", true_r0), #model_type
          #main = paste(seed_count, "Infts SS Evnts, ", "r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC Trace Plots
  r0_min = min(min(r0_mcmc), true_r0)
  r0_max = max(max(r0_mcmc), true_r0)
  plot.ts(r0_mcmc, ylab = 'alpha', #ylim=c(0, a_lim),
          ylim=c(r0_min, r0_max),
          main = paste("MCMC SS Events, true r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = true_r0, col = 'orange', lwd = 2) #True = green
  
  #iii. Cumulative mean plots
  #r0 Mean
  r0_min2 = min(min(r0_mean), true_r0)
  r0_max2 = max(max(r0_mean), true_r0)
  plot2 = plot(seq_along(r0_mean), r0_mean,
               ylim=c(r0_min2, r0_max2),
               xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
               cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  print(plot2)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #iv. Histogram
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = paste('R0 total MCMC samples'), #. Prior = ', prior),
       xlim=c(r0_min, r0_max),
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

#*##############################################
#*******************************************
#*
#* GRID PLOT SUPER-SPREADING MODELS:
#* 
#* FUNCTION TO PLOT 4x4 DASHBOARD OF MCMC RESULTS FOR SUPER SPREADING MODELs
#* 
#*******************************************
#*##############################################
plot_mcmc_grid <- function(n_mcmc, sim_data, mcmc_params, true_r0, total_time,
                           seed_count, model_params,
                           model_typeX = 'SSE', prior = TRUE, joint = TRUE,
                           flag_gam_prior_on_b = FALSE, gam_priors_on_b = c(0,0), rjmcmc = FALSE,
                           data_aug = FALSE,
                           mod_par_names = c('alpha', 'beta', 'gamma')){
  #Plot
  #plot.new()
  par(mfrow=c(4,4))
  
  #Extract params
  m1_mcmc = mcmc_params[1]; m1_mcmc = unlist(m1_mcmc)
  m2_mcmc = mcmc_params[2]; m2_mcmc = unlist(m2_mcmc)
  m3_mcmc = mcmc_params[3]; m3_mcmc = unlist(m3_mcmc)
  r0_mcmc = mcmc_params[4]; r0_mcmc = unlist(r0_mcmc)
  #True vals
  m1X = model_params[1]; m2X = model_params[2]; m3X = model_params[3]; 
  
  #Cumulative means + param sample limits
  #r0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  r0_lim = max(true_r0, max(r0_mcmc))
  r0_lim2 = max(true_r0, r0_mean)
  
  #m1
  m1_mean = cumsum(m1_mcmc)/seq_along(m1_mcmc)
  a_lim =  max(m1X, max(m1_mcmc))
  a_lim2 =  max(m1X, m1_mean)
  
  #m2
  m2_mean = cumsum(m2_mcmc)/seq_along(m2_mcmc)
  b_lim = max(m2X, max(m2_mcmc))
  b_lim2 = max(m2X, m2_mean)
  
  #m3
  m3_mean = cumsum(m3_mcmc)/seq_along(m3_mcmc)
  m3_lim =  max(m3X, max(m3_mcmc))
  m3_lim2 =  max(m3X, m3_mean) 
  
  #Priors
  if (flag_gam_prior_on_b) {
    m2_prior = paste0('m3(',  gam_priors_on_b[1], ', ', gam_priors_on_b[2], ')')
  } else {
    m2_prior = 'exp(1)'
  }
  m1_prior = 'exp(1)'
  m3_prior = '1 + exp(1)'
  
  #***********
  #* Plots *
  
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = paste(seed_count, ' Day Infts, ', model_typeX, "Data, r0 = ", true_r0),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC Trace Plots 
  plot.ts(m1_mcmc, ylab = mod_par_names[1], ylim=c(0, a_lim),
          main = paste("MCMC", model_typeX, ":", mod_par_names[1], "prior:", m1_prior),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = m1X, col = 'red', lwd = 2) #True = green
  
  plot.ts(m2_mcmc, ylab = 'm2', ylim=c(0, b_lim), 
          main = paste("MCMC", model_typeX, ":", mod_par_names[2], "prior:", m2_prior),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = m2X, col = 'blue', lwd = 2) #True = green
  
  plot.ts(m3_mcmc,  ylab = 'm3', ylim=c(0,m3_lim),
          main = paste("MCMC", model_typeX, ":", mod_par_names[3], "prior:", m3_prior),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = m3X, col = 'green', lwd = 2) #True = green
  
  #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
  
  #iii. Cumulative mean plots
  #r0 Mean
  plot(seq_along(r0_mean), r0_mean,
       ylim=c(0, r0_lim),
       xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #print(plot2)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #m1 mean
  plot(seq_along(m1_mean), m1_mean,
       ylim=c(0, a_lim),
       xlab = 'Time', ylab =  mod_par_names[1],
       main = paste(mod_par_names[1], "MCMC mean, True", mod_par_names[1], "=", m1X),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #print(plot2)
  abline(h = m1X, col = 'red', lwd = 2)
  
  #m2 mean
  plot(seq_along(m2_mean), m2_mean,
       ylim=c(0, b_lim),
       xlab = 'Time', ylab = 'm2',
       main = paste(mod_par_names[2], "MCMC mean, True", mod_par_names[2], "=", m2X),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #print(plot2)
  abline(h = m2X, col = 'blue', lwd = 2)
  
  #m3 Mean
  plot(seq_along(m3_mean), m3_mean,
       xlab = 'Time', ylab = 'm3', 
       main = paste(mod_par_names[3], "MCMC mean, True", mod_par_names[3], "=", m3X),
       ylim=c(0, m3_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = m3X, col = 'green', lwd = 2)
  
  #iv. Param Histograms (Plots 9,11,12)
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = paste('R0 total MCMC samples. Prior = ', prior),
       xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_r0, col = 'orange', lwd = 2)
  
  # #v. m2 vs m3
  # plot(m2_mcmc, m3_mcmc,
  #      xlab = 'm2', ylab = 'm3', main = 'm2 vs m3',
  #      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Hist m1 
  hist(m1_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[1], #ylab = 'Density', 
       main = paste(mod_par_names[1], ", True", mod_par_names[1], "=", m1X), 
       xlim=c(0, a_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = m1X, col = 'red', lwd = 2)
  
  #Hist m2 
  hist(m2_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[2], #ylab = 'Density', 
       main = paste(mod_par_names[2], ", True", mod_par_names[2], "=", m2X), 
       xlim=c(0, b_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = m2X, col = 'blue', lwd = 2)
  
  #Hist m3 
  hist(m3_mcmc, freq = FALSE, breaks = 100,
       xlab = mod_par_names[3], #ylab = 'Density', 
       main = paste(mod_par_names[3], ", True", mod_par_names[3], "=", m3X),
       xlim=c(0, m3_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = m3X, col = 'green', lwd = 2)
  
  #Final Mean Stats
  data_10_pc = 0.5*n_mcmc #50%
  a_mcmc_mean = round(mean(m1_mcmc[n_mcmc - data_10_pc:n_mcmc]), 2) 
  b_mcmc_mean = round(mean(m2_mcmc[n_mcmc - data_10_pc:n_mcmc]), 2)
  g_mcmc_mean = round(mean(m3_mcmc[n_mcmc - data_10_pc:n_mcmc]), 2)
  r0_mcmc_mean = round(mean(r0_mcmc[n_mcmc - data_10_pc:n_mcmc]), 2)
  
  #Joint distrbutions
  if (joint){
    
    #v. r0 vs m2
    plot(m2_mcmc, r0_mcmc,
         xlab = mod_par_names[2], ylab = 'R0', main = paste(mod_par_names[2], 'vs R0'),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. m1 vs m2
    plot(m1_mcmc, m2_mcmc,
         xlab = mod_par_names[1], ylab = mod_par_names[2], main = paste(mod_par_names[1], 'vs', mod_par_names[2]),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. m1 vs m3
    plot(m1_mcmc, m3_mcmc,
         xlab = mod_par_names[1], ylab = mod_par_names[3], main = paste(mod_par_names[1], 'vs', mod_par_names[3]),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    
    #v. m2 vs m3
    plot(m2_mcmc, m3_mcmc,
         xlab = mod_par_names[2], ylab = mod_par_names[3], main = paste(mod_par_names[2], 'vs', mod_par_names[3]),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }
  
  #*****************
  #RJMCMC
  if (rjmcmc){
    
    #m1S
    par(mfrow=c(3,1))
    
    #HIST m1 
    hist(m1_mcmc, freq = FALSE, breaks = 100,
         xlab = mod_par_names[1], #ylab = 'Density', 
         main = paste("m1, True m1 = ", m1X), 
         xlim=c(0, a_lim),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = m1X, col = 'red', lwd = 2)
    
    #m1 SSE
    ind_b = which(m2_mcmc > 0)
    m1_sse = m1_mcmc[ind_b]
    
    #Hist
    hist(m1_sse, freq = FALSE, breaks = 100,
         xlab = 'm1_sse', #ylab = 'Density', 
         main = "m1_sse", 
         xlim=c(0, a_lim),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = m1X, col = 'red', lwd = 2)
    
    #m1 BASE
    ind_b = which(m2_mcmc == 0)
    
    if (length(ind_b > 2)){
      
      m1_base = m1_mcmc[ind_b]
      #Hist
      hist(m1_base, freq = FALSE, breaks = 100,
           xlab = 'm1_base', #ylab = 'Density', 
           main = "m1_base",
           xlim=c(0, a_lim),
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(v = m1X, col = 'red', lwd = 2)
    }
    
    #RESULTS DF
    df_results <- data.frame(
      rep = seed_count,
      n_mcmc = n_mcmc,
      m1 = m1X,
      a_mc = a_mcmc_mean,
      m2 = m2X,
      b_mc = b_mcmc_mean,
      m3 = m3X,
      g_mc = g_mcmc_mean,
      R0 = true_r0, 
      R0_mc = r0_mcmc_mean,
      accept_rate_m1 = round(mcmc_params[[5]],2),
      a_rte_m2 = round(mcmc_params[[6]], 2),
      n_accept_m2 = mcmc_params[[15]],
      a_rte_m3 = round(mcmc_params[[7]],2),
      n_accept_m3 = mcmc_params[[16]],
      a_rte_m2_m3 = round(mcmc_params[[8]],2),
      n_accept_2_3 = mcmc_params[[17]],
      n_accept_rj0 = mcmc_params[[11]],
      n_reject_rj0 = mcmc_params[[13]],
      a_rte_rj0 = round(mcmc_params[[9]],2),
      n_accept_rj1 = mcmc_params[[12]],
      n_reject_rj1 = mcmc_params[[14]],
      a_rte_rj1 = round(mcmc_params[[10]],2),
      m2_pc0 = mcmc_params[[18]],
      m2_pc_non_0 = 1- mcmc_params[[18]],
      bf = mcmc_params[[19]])
    #tot_time = total_time)
    
  } else if (data_aug) {
    
    df_results <- data.frame(
      rep = seed_count,
      n_mcmc = n_mcmc,
      m1 = m1X,
      m1_mc = a_mcmc_mean,
      m2 = m2X,
      m2_mc = b_mcmc_mean,
      m3 = m3X,
      m3_mc = g_mcmc_mean,
      R0 = true_r0, 
      R0_mc = r0_mcmc_mean,
      accept_rate_m1 = round(mcmc_params[[5]],2),
      a_rte_m2 = round(mcmc_params[[6]], 2),
      a_rte_m3 = round(mcmc_params[[7]],2),
      a_rte_m2_m3 = round(mcmc_params[[8]],2),
      a_rte_d_aug = round(mcmc_params[[12]],2))
      #tot_time = total_time)
    
    } else {
    df_results <- data.frame(
      rep = seed_count,
      n_mcmc = n_mcmc,
      m1 = m1X,
      m1_mc = a_mcmc_mean,
      m2 = m2X,
      m2_mc = b_mcmc_mean,
      m3 = m3X,
      m3_mc = g_mcmc_mean,
      R0 = true_r0, 
      R0_mc = r0_mcmc_mean,
      accept_rate_m1 = round(mcmc_params[[5]],2),
      a_rte_m2 = round(mcmc_params[[6]], 2),
      a_rte_m3 = round(mcmc_params[[7]],2),
      a_rte_m2_m3 = round(mcmc_params[[8]],2),
      tot_time = total_time)
  }
  
  print(df_results)
  
}
