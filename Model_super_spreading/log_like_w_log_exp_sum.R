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


#***********************************************************************************************************
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

#******************************************************************************************************************
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
num_days = 15
x = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, gammaX, betaX)
logl_1 = log_like_ss(x, alphaX, betaX, gammaX)
print(logl_1)

#log exp sum trick
logl_2 = log_like_ss_lse(x, alphaX, betaX, gammaX)
print(logl_2)
