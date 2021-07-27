#Experimentation

#Log sum exp
y_t = 0
t = 2
x[t]
a = exp(-alphaX*lambda_t)*(1/factorial(y_t))*(alphaX*lambda_t)^y_t*
  (gamma((x[t] - y_t) + betaX*lambda_t))/(gamma(betaX*lambda_t)*
                                            factorial(x[t] - y_t))*(1/(gammaX +1))^(betaX*lambda_t)*
  (gammaX/(gammaX + 1))^(x[t] - y_t)



#***********************************#*******************************************
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
      logl = logl + log(exp(-alphaX*lambda_t)*(1/factorial(y_t))*(alphaX*lambda_t)^y_t*
        (gamma((x[t] - y_t) + betaX*lambda_t))/(gamma(betaX*lambda_t)*
                                                  factorial(x[t] - y_t))*(1/(gammaX +1))^(betaX*lambda_t)*
        (gammaX/(gammaX + 1))^(x[t] - y_t))
      print('logl')
      print(logl)
      
    } else {
      
      #Terms in inner sum
      inner_sum_vec <- vector('numeric', x[t])
      
      for (y_t in 1:x[t]){ #Sum for all values of y_t up to x_t
        
        #Store inner product in vector position
        inner_sum_vec[y_t] = exp(-alphaX*lambda_t)*(1/factorial(y_t))*(alphaX*lambda_t)^y_t*
          (gamma((x[t] - y_t) + betaX*lambda_t))/(gamma(betaX*lambda_t)*
                                                    factorial(x[t] - y_t))*(1/(gammaX +1))^(betaX*lambda_t)*
          (gammaX/(gammaX + 1))^(x[t] - y_t)
        
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
      
      print('lse')
      print(innersum2)
      
      print('log_lse')
      print(log(innersum2))
      
      print('logl')
      print(logl)
      
      print('x_max')
      print(x_max)
      
      
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





#**************************
#*Neyman type A distribution
seq1 = seq(0.0, 10, by = 1)
neyAdist = dCompound(seq1, parent = "pois", compound = "neymantypea", compoundDist = "neymantypea")

neyAdist = dCompound(seq1, parent = "poisson", compound = "poisson", compoundDist = "neymantypea", params = c(1,1), shape1 = 1, shape2 = 1)

neyAdist = dCompound(seq1, parent = "poisson", compoundDist = "neymantypea", params = c(1,1))

plot(seq1, neyAdist)

seq2 = seq(-1, 1, by = 0.2)
params<-c(4,5)
pgfDneymantypea(seq2,params)

#***********

for (y_t in 1:1){
  
  print(y_t)
  
}
