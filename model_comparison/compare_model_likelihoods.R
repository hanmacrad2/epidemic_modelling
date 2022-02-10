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

log_like_B0 <- function(y, alphaX) {
  
  #Params
  num_days = length(y)
  shape_gamma = 6
  scale_gamma = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days - 1)),                                                                                    shape = shape_gamma, scale = scale_gamma)
  logl = 0
  
  for (t in 2:num_days) {
    
    #Data
    y_t = y[t]
    lambda_t = sum(y[1:(t - 1)] * rev(prob_infect[1:(t - 1)]))
    
    if (y_t == 0) {
      logl = logl - (alphaX * lambda_t)
      
    } else {
      #Add to log likelihood
      logl = logl +  y_t * log(alphaX * lambda_t) - (alphaX * lambda_t) -
        lfactorial(y_t) #- lfactorial(x[t] - y_t))
      
    }
    
  }
  
  logl
  
}

#SSE
#All parameters
log_like_ss_lse(x, alphaX, betaX, gammaX)

#Parameters == 0
logl = log_like(x, alphaX)

#TO DO
#1. FIX SUM IN THE FUNCITON ABOVE - NO SUM JUST XT: DONE
#2. CHECK REPORT - CONSTANTS OF PROPORTIONALITY: DONE 
#3. COMPARISONS

#3i. COMPARISON 1: BASE
#LOG_LIKE vs log_like_ss_lse_B0

#3ii. COMPARISON 2: B = ZERO VS SMALL 
#log_like_ss_lse_B0 VS SMALL VALUES IN log_like_ss_lse0

#*************************************************************
#COMPARISON I - :D
#OPTION 1
loglikeBO_1 = log_like_B0(sim_data, alphaX)
loglikeBO_1

#OPTION 2
loglikeBO_2 = log_like(sim_data, alphaX)
loglikeBO_1

#COMPARISON 2 - :D
log2 = log_like_ss_lse(sim_data, alphaX, 0.0001, 0.0001)
log2
log2b = log_like_B0(sim_data, alphaX)
log2b
