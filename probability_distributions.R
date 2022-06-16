#DISTRIBUTIONS
#Explain

#d - Probability density 
#p - Cumulative distribution
#q - Quantile function
#r - random sample generation from the distribution

#Rate
lambda = 1.8 #1.8 #10 #1.8

#****************************************************************
#1. EXPONENTIAL DISTRIBUTION
#****************************************************************

#1. Probability Density 
a = seq(0, 1.5, by = 0.005)
plot(a, dexp(a, 1))

plot(a, dexp(a, 25))

plot.ts(dexp(1:100, lambda))
plot.ts(1- dexp(1:100, lambda))

#2. Cumulative Distribution Function
plot.ts(qexp(1:10, lambda))

#3. Random samples
plot.ts(rexp(100, lambda))
plot.ts(rexp(100, 10))

#****************************************************************
#2. POISSON DISTRIBUTION
#****************************************************************
#Samples
plot.ts(rpois(100, lambda))
#Probability Density 
plot.ts(dpois(1:20, lambda))
#Cumulative Distribution Function
plot.ts(qpois(1:10, lambda))

#****************************************************************
#3. NEGATIVE BINOMIAL DISTRIBUTION
#****************************************************************
betaX = 0.2
gammaX = 10
plot.ts(rnbinom(100, betaX*poi_rate, 1/(1 + gammaX)))
plot.ts(dnbinom(1:20, betaX*poi_rate, 1/(1 + gammaX)))


#****************************************************************
#4. GAMMA DISTRIBUTION
#****************************************************************
#Mean = shape*scale
plot.new()
shape_gamma = 6 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 1 #1/rate. 

shape_gamma = 10 
scale_gamma = 1/100 

shape_gamma = 10 
scale_gamma = 1 

plot.ts(rgamma(100, shape = shape_gamma, scale = scale_gamma))
plot.ts(dgamma(1:20, shape = shape_gamma, scale = scale_gamma))

plot.ts(dgamma(100, shape = shape_gamma, scale = scale_gamma))

#PLOT
xseq = seq(0, 0.3, length.out = 500)
plot(xseq, dgamma(xseq, shape = shape_gamma, scale = scale_gamma), type = 'l', #,
     ylab = 'Gamma(10, 0.01)', main = 'Gamma(10, 0.01)', lwd = 2, col = 'blue')

xseq = seq(0, 35, length.out = 1000)
plot(xseq, dgamma(xseq, shape = shape_gamma, scale = scale_gamma), type = 'l', #,
     ylab = 'Gamma(10, 0.01)', main = 'Gamma(10, 0.01)', lwd = 2, col = 'red')

#Option 2
shape_gamma = 10 #Gamma params for infectiousness curve (lambda) distribution
scale_gamma = 1/100 #1/rate. 

#CENTRED ON THE MEAN :D
a = seq(0, 0.5, by = 0.00005)
plot(a, dgamma(a, shape = shape_gamma, scale = scale_gamma, log = TRUE), #,
        ylab = 'Gamma(10, 0.01)', main = 'Gamma(10, 0.01), log = TRUE')

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

#**************
#*PRIOR PLOTS
shape_gamma = 10 
scale_gamma = 1/100 

xseq = seq(0, 0.3, length.out = 500)
plot(xseq, dgamma(xseq, shape = shape_gamma, scale = scale_gamma),  #,
     ylab = 'Gamma(10, 0.01)', main = 'Gamma(10, 0.01)',
     type = 'l', lwd = 2, col = 'blue')

#c
shape_gamma = 10 
scale_gamma = 1 

xseq = seq(0, 35, length.out = 1000)
lines(xseq, dgamma(xseq, shape = shape_gamma, scale = scale_gamma), #,
     ylab = 'Gamma(10, 0.01)', main = 'Gamma(10, 0.01)',
     type = 'l', lwd = 2, col = 'green')

xseq = seq(0, 1.5, length.out = 500)
lines(xseq, dexp(xseq, 25),
     type = 'l', lwd = 2, col = 'red')
