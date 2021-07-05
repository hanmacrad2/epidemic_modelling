#Log likelihood with additional function using log_exp_sum trick
source("1b_simulation_super_spreading_events.R")

#Parameters
n = 50000
num_days = 60 #100
shape_gamma = 6
scale_gamma = 1
#Priors
prior_alpha_k = 1
prior_alpha_theta = 1


#***********************************
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
    
    lambda_t = sum(x[1:t-1]*rev(prob_infect[1:t-1]))
    inner_sum_xt = 0
    
    for (y_t in 0:x[t]){ #Sum for all values of y_t
      
      #Log likelihood
      inner_sum_xt = inner_sum_xt + exp(-alphaX*lambda_t)*(1/factorial(y_t))*(alphaX*lambda_t)^y_t*
        (gamma((x[t] - y_t) + betaX*lambda_t))/(gamma(betaX*lambda_t)*
                                                  factorial(x[t] - y_t))*(1/(gammaX +1))^(betaX*lambda_t)*
        (gammaX/(gammaX + 1))^(x[t] - y_t)
      
    } 
    print('x[t]')
    print(x[t])
    logl = logl + log(inner_sum_xt) 
    print('inner_sum_xt')
    print(inner_sum_xt)
    print('log_inner_sum_xt')
    print(log(inner_sum_xt))
    print('logl')
    print(logl)
  }
  
  logl
  
}

#Apply
num_days = 10
x = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, gammaX, betaX)
x
logl_1 = log_like_ss(x, alphaX, betaX, gammaX)
logl_1

#******************************************************************************8
#Log Likelihood - log-exp-sum trick 
log_like_ss_lse <- function(x, alphaX, betaX, gammaX){
  
  #Params
  num_days = length(x)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  logl = 0 #0.000001
  
  for (t in 2:num_days) {
    
    print(t)
    lambda_t = sum(x[1:t-1]*rev(prob_infect[1:t-1]))
    print('x[t]')
    print(x[t])
    
    if(x[t] == 0){
      
      y_t = 0
      
      #Store inner product in vector position
      logl = logl -(alphaX*lambda_t) - lgamma(betaX*lambda_t) - 
        (betaX*lambda_t*log(gammaX +1))
      
      print('logl')
      print(logl)
      
    } else {
      
      #Terms in inner sum
      inner_sum_vec <- vector('numeric', x[t])
      
      for (y_t in 1:x[t]){ #Sum for all values of y_t up to x_t
        
        #Store inner product in vector position L(x_i)
        inner_sum_vec[y_t] = -(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) 
          + lgamma((x[t] - y_t)*betaX*lambda_t) - lgamma(betaX*lambda_t) - 
                                                    lfactorial(x[t] - y_t) - (betaX*lambda_t*log(gammaX +1)) + 
          (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1)
        
        
        #inner_sum_vec[y_t] = exp(-alphaX*lambda_t)*(1/factorial(y_t))*(alphaX*lambda_t)^y_t*
          #(gamma((x[t] - y_t) + betaX*lambda_t))/(gamma(betaX*lambda_t)*
                                                    #factorial(x[t] - y_t))*(1/(gammaX +1))^(betaX*lambda_t)*
          #(gammaX/(gammaX + 1))^(x[t] - y_t)
        
      }
      
      
      #Calculate max element in inner vector, for all y_t for a given t, x[t]
      x_max = max(inner_sum_vec)
      
      #Calculate lse
      innersum2 = 0
      for (i in 1:length(inner_sum_vec)){
        innersum2 = innersum2 + exp(inner_sum_vec[i] - x_max)
      }
      lse = x_max + log(innersum2)
      
      #Add to overall log likelihood 
      logl = logl + lse 
      
      #print('lse')
      #print(innersum2)
      
      #print('log_lse')
      #print(log(innersum2))
      
      print('logl')
      print(logl)
      
      #print('x_max')
      #print(x_max)
      
    }
    
  }
  
  logl
  
}

#Apply
#num_days = 10
#x = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, gammaX, betaX)
#x
logl_1 = log_like_ss_lse(x, alphaX, betaX, gammaX)
logl_1
