#*******************************************************
#Super-spreading simulation

simulation_super_spreaders = function(num_days, shape_gamma, scale_gamma, aX, bX, ss_mult) {
  'Simulate an epidemic with Superspreading events
  prop_ss = Proportion of superspreading days
  magnitude_ss = increased rate of superspreading event'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nss_infecteds = vector('numeric', num_days)
  ss_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 2
  nss_infecteds[1] = 2
  ss_infecteds[1] = 0
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((nss_infecteds[1:(t-1)] + ss_mult*ss_infecteds[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nss_infecteds[t] = rpois(1, aX*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss_infecteds[t] = rpois(1, bX*lambda_t)
    total_infecteds[t] = nss_infecteds[t] + ss_infecteds[t]
  }
  
  total_infecteds
}

#Apply
n = 50000
num_days = 15 #60 #100
shape_gamma = 6
scale_gamma = 1
aX = 2
bX = 2
ss_mult = 5
data_ss = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, aX, bX, ss_mult)


#Plot
plot.ts(data_ss)

#********************************************************************************
#Get individual 

get_list_super_spreaders = function(num_days, shape_gamma, scale_gamma, aX, bX, ss_mult) {
  'Simulate an epidemic with Superspreading events
  prop_ss = Proportion of superspreading days
  magnitude_ss = increased rate of superspreading event'
  
  #Set up
  total_infecteds = vector('numeric', num_days)
  nss_infecteds = vector('numeric', num_days)
  ss_infecteds = vector('numeric', num_days)
  total_infecteds[1] = 2
  nss_infecteds[1] = 2
  ss_infecteds[1] = 0
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Regular infecteds (tot_rate = lambda) fix notation
    lambda_t = sum((nss_infecteds[1:(t-1)] + ss_mult*ss_infecteds[1:(t-1)])*rev(prob_infect[1:(t-1)])) #?Why is it the reversed probability - given the way prob_infect is written. Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    nss_infecteds[t] = rpois(1, aX*lambda_t) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
    ss_infecteds[t] = rpois(1, bX*lambda_t)
    total_infecteds[t] = nss_infecteds[t] + ss_infecteds[t]
  }
  
  return (list(total_infecteds, nss_infecteds, ss_infecteds))
}


#Implement
num_days = 30
shape_gamma = 6
scale_gamma = 1
aX = 2
bX = 5 #2
ss_mult = 3 #5
start_time = Sys.time()
ss_count_params = get_list_super_spreaders(num_days, shape_gamma, scale_gamma, aX, bX, ss_mult)
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

plot.ts(x, ylab = "Daily infections Count", main = 'Daily Infections Count - Super-Spreaders model',
        col = 'blue', lwd = 2.5, xlim=c(2,29))
lines(y, col = 'green', lwd = 2.5)
lines(z, col = 'red', lwd = 2.5)
#legend
legend('topleft', legend=c("Total", "Non Super-Spreaders", "Super-Spreaders"),
       #bty = "n", # Removes the legend box,
       col=c("blue", "green", "red"), lty=1:3, pch = c(3,3))


#Cumulative data
cum_data_x <- cumsum(x)
cum_data_y <- cumsum(y)
cum_data_z <- cumsum(z)

plot.ts(cum_data_x, ylab = "Cumulative Infections", main = 'Cumulative Count of Infections',
        col = 'blue', lwd = 3, xlim=c(2,29))
lines(cum_data_y, col = 'green', lwd = 2.5)
lines(cum_data_z, col = 'red', lwd = 2.5)
legend('topleft', legend=c("Total", "Non Super-Spreading Events", "Super-Spreading Events"),
       #bty = "n", # Removes the legend box,
       col=c("blue", "green", "red"), lty=1:3, pch = c(3,3))