#Simulate Branching Process
#par(mar=c(1,1,1,1))

#Parameters
num_days = 30 #100
shape_gamma = 6
scale_gamma = 1
alpha = 1 #Without ss event, ~r0. 

#Function
simulate_branching_ss = function(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX) {
  'Simulate an epidemic with Superspreading events
  prop_ss = Proportion of superspreading days
  magnitude_ss = increased rate of superspreading event'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nsse_infecteds = vector('numeric', num_days)
  sse_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 1
  nsse_infecteds[1] = 1
  sse_infecteds[1] = 0
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda = sum(nsse_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability
    tot_rate = alpha*lambda #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nsse_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, beta*lambda) #Number of super-spreading events (beta)
    
    if (n_t > 0){
      sse_infecteds[t] = rpois(1, gammaX*n_t) #z_t: Total infecteds due to super-spreading event - num of events x Num individuals
    }
    
    total_infecteds[t] = nsse_infecteds[t] + sse_infecteds[t]
  }
  
  total_infecteds
}

#Implement
betaX = 10 
gammaX = 5
start_time = Sys.time()
x = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, gammaX, betaX)
end_time = Sys.time()
time_elap = end_time - start_time
#print(time_elap)
x

#Plots
plot.ts(x, ylab = "N Daily infections", main = 'Simulation of Super-spreading events - Daily Infections')
cum_data <- cumsum(x)
plot.ts(cum_data, ylab = "Cumulated infections", main = 'Simulation of Super-spreading events - Cumulated Infections')


#*********************************************************************
#Plots
plot_variations = function(list_r0, list_scale_shape, num_days){
  #Plot
  par(mfrow=c(3,5))
  for (shape_scale in list_scale_shape){
    shape_gamma = shape_scale[[1]]
    scale_gamma = shape_scale[[2]]
    for (r0 in list_r0){
      num_daily = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
      plot.ts(num_daily,lwd = 2, ylab = "N Daily infections", col = 'purple')
    }
  }
}

#*********
#Params + Apply
num_days = 30
list_r0 = c(0.5, 1, 2, 4, 8)
#list_shape_scale = list(list(9,0.5) , list(7.5,1.0), list(6,1))
list_shape_scale = list(list(2,2), list(3,2), list(1.5, 2))
plot_variations(list_r0,  list_shape_scale, num_days)

#**********
#Inspect gamma density
seq1 = seq(0.0, 10, by = 1)
gammaX = dgamma(seq1, shape = 1.5, scale = 2)
plot(seq1, gammaX)


#***************
#*Negative Binomial - inspect
vs = c();

#Poisson of poissons
for (i in 1:10000){
  vs[i] = rpois(1, rpois(1,20))
}

print(c(mean(vs), var(vs)))

#Negative Binomial
vs2 = rnbinom(10000, 20, 0.5)
print(c(mean(vs), var(vs)))

hist(vs)
hist(vs2)

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
