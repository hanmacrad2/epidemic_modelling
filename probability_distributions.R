#Distributions
lambda = 1.8 #1.8 #10 #1.8

#1. Exponential distribution

#Samples
plot.ts(rexp(100, lambda))
#Probability Density 
plot.ts(dexp(1:100, lambda))
plot.ts(1- dexp(1:100, lambda))
#Cumulative Distribution Function
plot.ts(qexp(1:10, lambda))

#Poisson distribution

#Samples
plot.ts(rpois(100, lambda))
#Probability Density 
plot.ts(dpois(1:20, lambda))
#Cumulative Distribution Function
plot.ts(qpois(1:10, lambda))

#Negative Binomial Distribution

betaX = 0.2
gammaX = 10
plot.ts(rnbinom(100, betaX*poi_rate, 1/(1 + gammaX)))

plot.ts(dnbinom(1:20, betaX*poi_rate, 1/(1 + gammaX)))

#Gamma
#Mean = shape*scale
shape_gamma = 6 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 1 #1/rate. 

plot.ts(rgamma(100, shape = shape_gamma, scale = scale_gamma))
plot.ts(dgamma(1:20, shape = shape_gamma, scale = scale_gamma))

#Example number of goals in a match
shape_football = 25
scale_football = 1/10

plot.ts(rgamma(100, shape = shape_football, scale = scale_football))
plot.ts(dgamma(1:20, shape = shape_football, scale = scale_football))




