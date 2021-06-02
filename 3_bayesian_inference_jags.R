##To do
#1. Download JAGs
#2. Inspect results

#Packages
library(R2OpenBUGS)
library(rjags)
library(coda)
library(MCMCvis)

source("1_simulation.R")
source("3_bayesian_inference_jags.R")

#Params
r0 = 3.5
num_days = 45
shape_gamma = 6
scale_gamma = 1

#Data
y = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
  

#Model
cat("model
    {

  #Likelihood
  #Params
  num_days = length(y)

  #Infectiousness
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  likel = 1

  for (t in 2:num_days) {

    lambda = r0*sum(y[1:t-1]*rev(prob_infect[1:t-1]))
    likel = likel*lambda^(y[t])*exp(-lambda)
    
  }
    
  #Priors
  r0~exp(1)

  }",file = "bayesian_model.txt")

#Initial values
inits <- NULL

#Parameters to track
params = c('r0')

#Hyperparameters
ni = 10000 # number of iterations
nb = 1000 # burn in interval
nt = 1 # thinning interval
nc = 3 # number of chains

#Compile model
jmod = jags.model(file = 'bayesian_model.txt', data = y, n.chains = nc, inits = inits, n.adapt = 1000)



#********
#Drafts
#Model
# model_bayesian = function(){
#   
#   #Priors
#   r0~exp(1)
#   
#   #Likelihood
#   #Params
#   num_days = length(y)
#   
#   #Infectiousness (Discrete gamma)
#   prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
#   likel = 1
#   
#   for (t in 2:num_days) {
#     
#     lambda = r0*sum(y[1:t-1]*rev(prob_infect[1:t-1]))
#     likel = likel*lambda^(y[t])*exp(-lambda)
#     
#   }
#   
# }

#Write Model
#model.file = "model.file"
#write.model(model_bayesian, model.file)
