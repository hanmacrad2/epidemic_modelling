#Compare likelihood outputs

setwd("~/GitHub/epidemic_modelling")
source("epidemic_functions.R")
source("helper_functions.R")

#Epidemic params
num_days = 50
shape_gamma = 6 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 1 
seed_count = 1

#***************
#BASE LIKELIHOOD

#Log Likelihood - log-exp-sum trick 
log_like_ss_lse_B0 <- function(x, alphaX, betaX, gammaX){
  
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
      
      #print(paste0('logl 1 = ', logl))
      
    } else {
      #Terms in inner sum
      inner_sum_vec <- vector('numeric', x[t])
      
      #CUT FOR LOOP; REPLACE ALL YT WITH XT (ONE ITERATION) LET YT = XT. 
      for (y_t in 0:x[t]){ #Sum for all values of y_t up to x_t
        #print(paste0('y_t  = ', y_t))
        #Store inner L(x_i) term in vector position
        # inner_sum_Lx = (-(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) +
        #                   lgamma((x[t] - y_t) + (betaX*lambda_t)) - lgamma(betaX*lambda_t) - 
        #                   lfactorial(x[t] - y_t) - (betaX*lambda_t*log(gammaX +1)) + 
        #                   (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1))
        
        inner_sum_Lx = (-(alphaX*lambda_t) - lfactorial(y_t) + y_t*log(alphaX*lambda_t) +
                          lgamma((x[t] - y_t) + (betaX*lambda_t)) #lgamma(betaX*lambda_t) - 
                        - lfactorial(x[t] - y_t)) #- (betaX*lambda_t*log(gammaX +1)) 
        #+ (x[t] - y_t)*log(gammaX) -(x[t] - y_t)*log(gammaX + 1))
        
        #print(paste0('inner_sum_Lx  = ', inner_sum_Lx))
        
        #Check inf
        if (is.infinite(inner_sum_Lx)){
          #print(paste0('inner_sum_Lx is inf:', inner_sum_Lx))
          #print(paste0('yt value = ', y_t))
        } else {
          inner_sum_vec[y_t + 1] = inner_sum_Lx
        }
        
        
      }
      
      #Calculate max element in inner vector, for all y_t for a given t, x[t]
      #print(paste0('inner_sum_vec = ', inner_sum_vec))
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


#Likelihoods
#All parameters
log_like_ss_lse(x, alphaX, betaX, gammaX)

#Parameters == 0
logl = log_like(x, alphaX)

#TO DO
#1. FIX SUM IN THE FUNCITON ABOVE - NO SUM JUST XT 
#2. CHECK REPORT - CONSTANTS OF PROPORTIONALITY
#3. COMPARISONS

#3i. COMPARISON 1: BASE
#LOG_LIKE vs log_like_ss_lse_B0

#3ii. COMPARISON 2: B = ZERO VS SMALL 
#log_like_ss_lse_B0 VS SMALL VALUES IN log_like_ss_lse




