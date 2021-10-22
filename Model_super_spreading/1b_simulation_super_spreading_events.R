#Simulate Branching Process
#par(mar=c(1,1,1,1))

#Parameters
num_days = 30 #100
shape_gamma = 6
scale_gamma = 1
alphaX = 1 #Without ss event, ~r0. 

#Function
simulate_ss_events = function(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX) {
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
    lambda_t = sum(nsse_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability
    tot_rate = alphaX*lambda_t #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nsse_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, betaX*lambda_t) #Number of super-spreading events (beta)
    
    if (n_t > 0){
      sse_infecteds[t] = rpois(1, gammaX*n_t) #z_t: Total infecteds due to super-spreading event - num of events x Num individuals
    }
    
    total_infecteds[t] = nsse_infecteds[t] + sse_infecteds[t]
  }
  
  total_infecteds
}

#Implement
alphaX = 3
betaX = 3
gammaX = 3
num_days = 30
start_time = Sys.time()
x = simulate_ss_events(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
end_time = Sys.time()
time_elap = end_time - start_time
#print(time_elap)
x

#Plots
#plot.ts(x, ylab = "N Daily infections",
#        col = 'orange', lwd = 2)
#cum_data <- cumsum(x)
#plot.ts(cum_data, ylab = "Cumulated infections")

#Return nsse & sse
get_list_ss_events = function(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX) {
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
    lambda_t = sum(nsse_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability
    tot_rate = alphaX*lambda_t #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nsse_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    
    #Super-spreaders
    n_t = rpois(1, betaX*lambda_t) #Number of super-spreading events (beta)
    
    if (n_t > 0){
      sse_infecteds[t] = rpois(1, gammaX*n_t) #z_t: Total infecteds due to super-spreading event - num of events x Num individuals
    }
    
    total_infecteds[t] = nsse_infecteds[t] + sse_infecteds[t]
  }
  
  total_infecteds
  
  return (list(total_infecteds, nsse_infecteds, sse_infecteds))
}

#Implement
alphaX = 3
betaX = 3
gammaX = 3
num_days = 30
start_time = Sys.time()
ss_count_params = get_list_ss_events(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
end_time = Sys.time()
time_elap = end_time - start_time

#Extract elements
x = ss_count_params[1]
x = unlist(x)

y = ss_count_params[2]
y = unlist(y)

z = ss_count_params[3]
z= unlist(z)

#Plots
par(mfrow=c(1,1))

plot.ts(x, ylab = "Daily infections Count", main = 'Daily Infections Count',
        col = 'blue', lwd = 2.5)
lines(y, col = 'green', lwd = 2.5)
lines(z, col = 'red', lwd = 2.5)
#legend
legend('topleft', legend=c("Total", "Non Super-Spreading Events", "Super-Spreading Events"),
       #bty = "n", # Removes the legend box,
       col=c("blue", "green", "red"), lty=1:3, pch = c(3,3))


#Cumulative data
cum_data_x <- cumsum(x)
cum_data_y <- cumsum(y)
cum_data_z <- cumsum(z)

plot.ts(cum_data_x, ylab = "Cumulative Infections", main = 'Cumulative Count of Infections',
        col = 'blue', lwd = 3)
lines(cum_data_y, col = 'green', lwd = 2.5)
lines(cum_data_z, col = 'red', lwd = 2.5)
legend('topleft', legend=c("Total", "Non Super-Spreading Events", "Super-Spreading Events"),
       #bty = "n", # Removes the legend box,
       col=c("blue", "green", "red"), lty=1:3, pch = c(3,3))