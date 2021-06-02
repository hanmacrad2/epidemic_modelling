#Simulate Branching Process
par(mar=c(1,1,1,1))

#Parameters
num_days = 100
r0 = 3.1
shape_gamma = 6
scale_gamma = 1


#Function
simulate_branching = function(num_days, r0, shape_gamma, scale_gamma) {
  #Set up
  vec_infecteds = vector('numeric', num_days)
  vec_infecteds[1] = 1
  
  #Infectiousness (Discrete gamma)
  prob_infect = pgamma(c(1:num_days), shape = shape_gamma, scale = scale_gamma) - pgamma(c(0:(num_days-1)), shape = shape_gamma, scale = scale_gamma)
  
  for (t in 2:num_days) {
    
    #Total rate
    tot_rate = r0*sum(vec_infecteds[1:t-1]*rev(prob_infect[1:t-1]))
    vec_infecteds[t] = rpois(1, tot_rate)
  }
  
  vec_infecteds
}

#Implement
start_time = Sys.time()
x = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
end_time = Sys.time()
time_elap = end_time - start_time
#print(time_elap)
x

#Plots
plot.ts(x, ylab = "N Daily infections")
cum_data <- cumsum(x)
plot.ts(cum_data)


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
gammaX = dgamma(a, shape = 1.5, scale = 2)
plot(seq1, gammaX)