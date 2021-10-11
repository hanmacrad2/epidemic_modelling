#Functions

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
    lambda_t = sum(nsse_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written
    tot_rate = alphaX*lambda_t #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nsse_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, betaX*lambda_t) #Number of super-spreading events (beta)
    
    if (n_t > 0){
      sse_infecteds[t] = rpois(1, gammaX*n_t) #z_t: Total infecteds due to super-spreading event - num of events x Num individuals
    }
    
    total_infecteds[t] = nsse_infecteds[t] + sse_infecteds[t]
  }
  
  total_infecteds
}

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
    lambda_t = sum((nss_infecteds[1:(t-1)] + ss_mult*ss_infecteds[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nss_infecteds[t] = rpois(1, aX*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss_infecteds[t] = rpois(1, bX*lambda_t)
    total_infecteds[t] = nss_infecteds[t] + ss_infecteds[t]
  }
  
  total_infecteds
}

#Apply
n = 50000
num_days = 15 #60 #100
shape_gamma = 6
scale_gamma = 1
aX = 2
bX = 2
ss_mult = 5
data_ss = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, aX, bX, ss_mult)

#*******************************************************************************
#Model functions  - Super-spreading Model

#Log Likelihood 
log_like_ss <- function(x, alphaX, betaX, gammaX){
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
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
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
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


#Apply
#Parameters
n = 50000
num_days = 15 #60 #100
shape_gamma = 6
scale_gamma = 1
#Priors
prior_alpha_k = 1
prior_alpha_theta = 1
#x = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, gammaX, betaX)

#Log likelihood 
#logl_1 = log_like_ss(x, alphaX, betaX, gammaX)
#print(logl_1)

#log exp sum trick
#logl_2 = log_like_ss_lse(x, alphaX, betaX, gammaX)
#print(logl_2)

#******************************************************
#Plotting
#**********************************************
#*
#Plots
plot_mcmc_super_spreading <- function(mcmc_vector1, mcmc_vector2, mcmc_vector3, alpha_true, betaX, gammaX, folder_dir_ad) {
  
  #Folder save
  pdf(paste(folder_dir_ad, "/", "ss_adpative_mc_alpha_true_", alpha_true, ".pdf", sep=""))
  
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector1, ylab = 'alpha', main = paste("MCMC Results. Alpha - Super spreading model, true alpha = ", alpha_true))
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
  
  #**************************************************
  #2. beta
  
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector2, ylab = 'beta', main = paste("MCMC Results; beta - Super spreading model, true beta = ", betaX))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  alpha_mean = cumsum(mcmc_vector2)/seq_along(mcmc_vector2)
  plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, true beta = ", betaX))
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
  plot1 = ts.plot(mcmc_vector3, ylab = 'gamma', main = paste("MCMC Results; gamma - Super spreading model, true gamma = ", gammaX))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  gamma_mean = cumsum(mcmc_vector3)/seq_along(mcmc_vector3)
  plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True gamma = ", gammaX))
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