#Description
#Metropolis algorithm with target r0

#Setup
#source("1_simulate_branching.R")
source("1_simulation.R")
par(mar=c(1,1,1,1))

#Params
num_days = 30 #100
r0_true = 3.1
shape_gamma = 6
scale_gamma = 1

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
    
    log_alpha = log_like(data, Y) - log_like(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) #log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    
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
data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
n = 20000
sigma =  0.5 #(2.38^2/dimesnion_paramter)*Posterior or sample variance ~optimal
r0_mcmc = MetropolisHastings_r0(data, n, sigma)

#Plots
#i. MCMC chain
ts.plot(r0_mcmc, ylab = 'R0', main = paste("MCMC of R0, true R0 = 3.1, sd of proposal = ", sigma))
#main = paste("Tau= ", tau, " Sigma= ", sigma),

#ii. Mean
#Plot mean
r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0', main = paste("Mean of R0 MCMC chain, True R0 = ",r0, ", sd of proposal = ", sigma))

#Histogram
hist(r0_mcmc, prob = TRUE)

#Hist
hist1 <- hist(r0_mcmc, breaks = 80)
hist1$counts <- hist1$counts/sum(hist1$counts)
plot(hist1, xlab = 'r0', ylab = 'Density', 
     main = 'Empirical density of r0 - MCMC chain')

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
  
  #Return r0, acceptance rate
  return(c(r0_vec, accept_rate))
}


#*************************************************************************************************

#Investigate range of sd on MCMC solution

MCMC_range_sd <- function(list_sd, data, n, r0_true){
  
  #Initialise
  par(mfrow = c(2, 1))
  list_accept_rate <- vector('numeric', length(list_sd))
  i = 1
  
  for (sigma in list_sd) {
    print('sigma =')
    print(sigma)
    #R0 - mcmc
    mcmc_params = MetropolisHastings_r0_vII(data, n, sigma)
    #Extract
    r0_mcmc = mcmc_params[1]
    accept_rate = mcmc_params[2]
    list_accept_rate[i] = accept_rate
    
    #Plots
    
    #i.MCMC chain
    ts.plot(r0_mcmc, ylab = 'R0', 
            main = paste("MCMC of R0, true R0 = ", r0_true, ", sd of proposal = ", sigma))
    
    #ii.Cumulative Mean
    r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
    plot(seq_along(r0_mean), r0_mean, xlab = 'Time', ylab = 'R0', main = paste("Mean of R0 MCMC chain, True R0 = ",r0, ", sd of proposal = ", sigma))
    
    #iii. Histograms
    hist1 <- hist(r0_mcmc, breaks = 80, xlab = 'r0', ylab = 'Density', 
                  main = paste("Histogram of r0 - MCMC chain. True R0 = ",r0, ", sd of proposal = ", sigma))
    #iiib. Histogram - density
    hist1$counts <- hist1$counts/sum(hist1$counts)
    plot(hist1, xlab = 'r0', ylab = 'Density', 
         main = paste("Empirical density of r0 - MCMC chain. True R0 = ",r0, ", sd of proposal = ", sigma))
    
    
  }
  
  #Create dataframe
  df_sd_results <- data.frame(
    sd = list_sd,
    acceptance_rate = list_accept_rate
    )
  
  df_sd_results
}

#Apply
list_sd = c(0.25, 0.5, 0.75) #1, 1.25, 1.5, 2, 2.5, 3
df_sd_results = MCMC_range_sd(list_sd, data, n, r0_true)


#Note
#Actual analytical (Not possible with real data - parameters unknown) alpha <- (dgamma(Y, post_gamma_shape, post_gamma_scale))/ (dgamma(r0_vec[i-1], post_gamma_shape, post_gamma_scale))


#***********************
#Bayesian Analytical plot
#Data
x = simulate_branching(num_days, r0, shape_gamma, scale_gamma)

#Posterior
get_posterior_params <- function(){
  
  #Data + params
  sum_lambda = 0
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) {
    
    lambda_t = sum(x[1:t-1]*rev(prob_infect[1:t-1]), na.rm = TRUE)
    sum_lambda = sum_lambda + lambda_t
    
  }
  
  #Posterior
  sum_x = sum(x)
  post_gamma_shape = sum_x + prior_r0_k
  #print(post_gamma_shape)
  post_gamma_scale = 1/(sum_lambda + (1/prior_r0_theta))
  
  
  return(c(post_gamma_shape, post_gamma_scale))
  
}

#Apply
post_params = get_posterior_params()
post_gamma_shape = post_params[1]
post_gamma_scale = post_params[2]

#Plot
a = seq(0.0, 20, by = 0.01)
b = rgamma(n, post_gamma_shape, post_gamma_scale)
hist2 = hist(b, breaks = 100)
hist2$counts <- hist2$counts/sum(hist2$counts)
plot(hist2, xlab = 'r0', ylab = 'Density of gamma(2, 0.5)', #shape, scale
     main = 'Density of r0 - Analytical Bayesian model')


#******************************************************************************
# Investigate affect of varying the sd in the Standard Normal proposal Distribution

#Apply MCMC
n = 20000
sigma = 1 #2.38 #optimal
data = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
r0_mcmc = MetropolisHastings_r0(data, n, sigma)

#Chain 
ts.plot(r0_mcmc, ylab = 'R0', main = 'MCMC chain of R0, proposal sd = 1')


#Check
df <- data.frame(
  sd = list_sd,
  acceptance_rate = list_accept_rate[1:3]
)
