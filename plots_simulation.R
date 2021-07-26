#Simulation - plots for presentation

#Simulate Branching Process
par(mar=c(1,1,1,1))

#Parameters
num_days = 30 #100
r0 = 3.1
shape_gamma = 6
scale_gamma = 1


#Function
simulate_branching = function(num_days, r0, shape_gamma, scale_gamma) {
  #Set up
  vec_infecteds = vector('numeric', num_days)
  vec_infecteds[1] = 2
  
  #Infectiousness (Discrete gamma) - I.e 'Infectiousness Pressure' - Sum of all people
  #Explanation: Gamma is a continuous function so integrate over the density at that point in time (today - previous day)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  #Days of Infection Spreading
  for (t in 2:num_days) {
    
    #Total rate
    tot_rate = r0*sum(vec_infecteds[1:(t-1)]*rev(prob_infect[1:(t-1)])) #Product of infecteds & their probablilty of infection along the gamma dist at that point in time
    vec_infecteds[t] = rpois(1, tot_rate) #Assuming number of cases each day follows a poisson distribution. Causes jumps in data 
  }
  
  vec_infecteds
}

#Implement
start_time = Sys.time()
x2 = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
end_time = Sys.time()
time_elap = end_time - start_time
#print(time_elap)
x2

#Plots
plot.ts(x2, ylab = "N Daily infections", main = 'Number of Daily Infections, R0 = 3.1')

#Cumulative data
cum_data <- cumsum(x2)
plot.ts(cum_data, ylab = "Cumulative infections", main = 'Cumulative Infections, R0 = 3.1',
        col = 'orange', lwd = 2)


#*********************************************************************
#Plots
#Plot range of R0

plot_range_r0 = function(list_r0, num_days, shape_gamma, scale_gamma, colorsX){
  #Plot
  for (i in seq_along(list_r0)) {
    r0X = list_r0[i]
    print(r0X)
    num_daily = simulate_branching(num_days, r0X, shape_gamma, scale_gamma)
    #Plot first val
    if (i == 1) {
      plot.ts(num_daily, lwd = 2, xlab = 'Time', ylab = "N Daily infections",
              main = expression(paste('No. Daily Infections, ' ~ Gamma, "(6,1)")),
              col = colorsX[i])
    } else {
      lines(num_daily, lwd = 2, col = colorsX[i])
    }
    
  }
  legend("topleft", inset=.03, legend = list_r0,
         col= colorsX, lwd = c(2,2,2,2,2), cex=0.8)
}

#Apply
num_days = 30
list_r0 = rev(c(1.0, 1.5, 2.0, 2.5, 3.0))
colorsX = c('orange', 'red', 'green', 'blue', 'cyan')
plot_range_r0(list_r0, num_days, shape_gamma, scale_gamma, colorsX)

#*********************************
#Plot range of gamma distributions

#Plot range of R0
plot_range_gamma = function(r0X, num_days, list_gamma_params, colorsX){
  i = 1
  #Plot
  for (shape_scale in list_gamma_params) {
    #Extract params
    #shape_scale = list_gamma_params[[i]]
    #print(shape_scale)
    shape_gamma = shape_scale[[1]]
    scale_gamma = shape_scale[[2]]
    print(shape_gamma)
    print(scale_gamma)
    
    num_daily = simulate_branching(num_days, r0X, shape_gamma, scale_gamma)
    
    #Plot first val
    if (i == 1) {
      plot.ts(num_daily, lwd = 2, xlab = 'Time', ylab = "N Daily infections",
              main = paste('No. Daily Infections, R0 = 2.0'),# r0X),
              col = colorsX[i])
    } else {
      lines(num_daily, lwd = 2, col = colorsX[i])
    }
    i = i + 1 
  }
  
  legend("topleft", inset=.03,
         legend = c(expression(paste(Gamma, "(2,2)")),
                    expression(paste(Gamma, "(3,2)")),
                    expression(paste(Gamma, "(1.5, 2)")),
                    expression(paste(Gamma, "(4,1)")),
                    expression(paste(Gamma, "(6,1)"))),
         col= colorsX, lwd = c(2,2,2,2,2), cex=0.8)
  
}
#list_gamma_params[1][[1]]
#Apply
num_days = 30
#list_r0 = rev(c(1.0, 1.5, 2.0, 2.5, 3.0))
colorsX = c('orange', 'red', 'green', 'blue', 'cyan')

list_gamma_params = list(list(2,2), list(3,2), list(1.5, 2), list(4, 1), list(6, 1))
r0X = 2.0
plot_range_gamma(r0X, num_days, list_gamma_params, colorsX)

#*************************************
#Plot variation of r0 on multiple plots
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
seq1 = seq(0.0, 9.5, by = 0.01)
gammaX = dgamma(seq1, shape = 1.5, scale = 2)
plot(seq1, gammaX, type =  'l', xlab = "Time", ylab = '', 
     main = expression(paste('Time-varying infectiousness of each individual,' ~ Gamma, "(1.5, 2)")),
     col = 'green', lwd = 2.5)
legend()


