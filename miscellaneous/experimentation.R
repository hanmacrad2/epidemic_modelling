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


#**********************************************************************
#*Adaptive Scaling Algorithm

adaptive_scaling_metropolis_r0 <- function(data, n, sigma, alpha_star, x0 = 1, burn_in = 5000) { #burn_in = 2500
  
  'Returns mcmc samples of R0 & acceptance rate'
  
  #Set up
  r0_vec <- vector('numeric', n)
  scaling_vec <- vector('numeric', n)
  r0_vec[1] <- x0
  scaling_vec[1] <- x0
  U <- runif(n)
  count_accept = 0
  count_reject = 0
  #sd_sample = 1
  
  #MCMC chain
  for(i in 2:n) {
    
    #New Proposal
    Y <- r0_vec[i-1] + exp(scaling_vec[i-1])*rnorm(1) #sd = sigma) exp
    
    #if (is.na(Y) || is.nan(Y) || is.infinite(Y)) next
    
    if(Y < 0){
      Y = abs(Y)
    }
    # print('y')
    # print(Y)
    # 
    # #Prints
    # print('log-likelihoods')
    # loglike1 = log_like(data, Y)
    # print(loglike1)
    # loglike2 = log_like(data, r0_vec[i-1])
    # print(loglike2)
    # dg1 = dgamma(Y, shape = 1, scale = 1, log = TRUE)
    # print('priors')
    # print(dg1)
    # dg2 = dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE)
    # print(dg2)
    
    log_alpha = log_likeII(data, Y) - log_likeII(data, r0_vec[i-1]) + dgamma(Y, shape = 1, scale = 1, log = TRUE) - dgamma(r0_vec[i-1], shape = 1, scale = 1, log = TRUE) #log_prior(theta_dash) - log_prior(theta) = 1 - 1 
    # print('log_alpha')
    # print(log_alpha)
    
    #if (is.na(log_alpha)){
    #print('na log_alpha value:')
    #print(log_alpha)
    #print('Y value:')
    #print(Y)
    #}
    #if(!(is.na(log_alpha)) && log(U[i]) < log_alpha) {
    if(log(U[i]) < log_alpha) {
      r0_vec[i] <- Y
      count_accept = count_accept + 1
    } else {
      r0_vec[i] <- r0_vec[i-1]
      count_reject = count_reject + 1
    }
    log_alpha = min(0, log_alpha)
    
    #Scaling factor
    #print('log_alpha')
    #print(log_alpha)
    scaling_vec[i] = scaling_vec[i-1] + (1/i)*(exp(log_alpha) - alpha_star)
    # print('scaling_vec')
    # print(scaling_vec[i])
    
    #Adaptive MC
    #if (i == burn_in){
    #  sigma = var(r0_vec[2:i])*(2.38^2)
    #}
    
  }
  #Final stats
  total_iters = count_accept + count_reject
  accept_rate = 100*(count_accept/(count_accept+count_reject))
  num_samples = count_accept
  print("Total iterations = ")
  print(total_iters)
  print("Acceptance rate = ")
  print(accept_rate)
  print("Number samples = ")
  print(count_accept)
  
  #Burn-in
  r0_vec = r0_vec[burn_in:n]
  
  #Return r0, acceptance rate
  return(list(r0_vec, accept_rate, num_samples))
}

#Apply
alpha_star = 0.40
as_params = adaptive_scaling_metropolis_r0(data, n, sigma, alpha_star, x0 = 1, burn_in = 5000)

r0_as = as_params[1]
r0_as = unlist(r0_as)

#3d Plots
#Plots

x2 = seq(-2,1,length=5) 

for (i in x2) {
  print(i)
}

#Test
u <- seq(-5, 5, by = .1)
v <- seq(-5, 5, by = .1)
M <- expand.grid(u,v)

x <- M$Var1
y <- M$Var2

sigma <- matrix(c(1, .5, .5, 1), nrow = 2, byrow = TRUE)
z <- dmvnorm(x = M, sigma = sigma)

scatterplot3js(x, y, z, phi = 40, theta = 20,
               color=rainbow(length(z)),
               colkey = FALSE,
               cex = .3,
               main = "Bivariate Normal")

scatterplot3js(X1, X2, z, phi = 40, theta = 20,
               color=rainbow(length(z)),
               colkey = FALSE,
               cex = .3,
               main = "f(x1, x2)")

#Plots
plot
library(plotly)

# Data: volcano is provided by plotly

# Plot
persp(X1, X2, fx[3,], 
      xlab = "x1", ylab = "x2",
      main = "f(x1, x2)"
)

library(rgl)
library(plot3D)
library(threejs)
library(mvtnorm)

surf3D(x = fx[1,], y = fx[2 ,], z = fx[3 ,], type = "surface")
rglwidget(elementId = "plot3drgl")

z = fx[3 ,]
z = matrix(fx[3 ,], nrow = 1, ncol = len_x^2)
fig <- plot_ly(x = fx[1,], y = fx[2 ,], z = z) %>% add_surface()

fig


#Other
scatterplot3js(X1, X2, z, phi = 40, theta = 20,
               color=rainbow(length(z)),
               colkey = FALSE,
               cex = .3,
               axisLabels=c("x1", "x2", "f(x1, x2"),
               main = "f(x1, x2)")

#Check
h = 0 
for (i in c(1,2,3,4)) {
  h[i] = i*2
}




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


