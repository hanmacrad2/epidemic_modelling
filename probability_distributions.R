#Distributions
lambda = 1.8 #1.8 #10 #1.8

#1. Exponential distribution

#Samples
plot.ts(rexp(100, lambda))
#Probability Density 
a = seq(0, 1.5, by = 0.005)
plot(a, dexp(a, 1))

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

#Option 2
shape_gamma = 2 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 0.5 #1/rate. 

a = seq(0, 1.5, by = 0.005)
plot(a, dgamma(a, shape = shape_gamma, scale = scale_gamma),
        ylab = 'Gamma(2, 0.5)', main = 'Gamma(2, 0.5)')

#V2
shape_gamma = 5
scale_gamma = 5
a = seq(0, 1.5, by = 0.005)
plot(a, dgamma(a, shape = shape_gamma, scale = scale_gamma),
     ylab = 'Gamma(5, 5)', main = 'Gamma(5, 5)')

plot.ts(dgamma(1:20, shape = shape_gamma, scale = scale_gamma),
        ylab = 'Gamma(2, 0.5)', main = 'Gamma(2, 0.5)')

plot.ts(rgamma(1000, shape = shape_gamma, scale = scale_gamma),
        ylab = 'Gamma(2, 0.5)', main = 'Gamma(2, 0.5)')

mean(rgamma(100, shape = shape_gamma, scale = scale_gamma))

#Example number of goals in a match
shape_football = 25
scale_football = 1/10

plot.ts(rgamma(100, shape = shape_football, scale = scale_football))
plot.ts(dgamma(1:20, shape = shape_football, scale = scale_football))


length(a)
