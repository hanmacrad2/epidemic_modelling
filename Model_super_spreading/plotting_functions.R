#Plotting Functions

#**********************************************************************************************
#Plotting
plot_points_comparison <- function(vec1, vec2){
  
  'Plot points on a diagonal line in ggplot for comparison, same x-axis'
  
  
  #Dataframe for plot
  group = c(rep.int(1, length(vec1)), rep.int(2, length(vec1)))
  results_total = c(vec1, vec2)
  vec_x_axis = c(seq_along(vec1), seq_along(vec1))
  df_results  = data.frame(group, vec_x_axis, results_total)
  
  #Plot
  ggplot(df_results, aes(x = vec_x_axis, y = results_total, color = group) ) + #3.5 == 3.41, 3 = 3.65,
    geom_point(shape=21, fill = group, size=3 ) +
    theme_bw() +
    xlab("Iteration") + ylab("Original (black), R packages (red)") +
    ggtitle("Log(acceptance_prob); Original (black), R packages(red)")
  
}

#Difference between lot points comparison
jk = 1:10
#plot_points_comparison(jk, jk^2)

#Difference between points
plot_diff_points_comparison <- function(vec1, vec2, titleX){
  
  'Plot points on a diagonal line in ggplot for comparison, same x-axis'
  plot(seq_along(vec1), abs(vec1 - vec2), type="b", pch=19, col="blue", main = titleX, xlab="index", ylab="Absolute difference")
  
}

#Plot points comparison
#jk = 1:10
#plot_diff_points_comparison(exp(jk), jk^3, 'Abs difference btwn log(acceptance_probabilites) of Original ver & R packages ver; ')

plot_points_comparison_straight <- function(vec1, vec2){
  
  'Plot points on ggplot for comparison, same x-axis'
  
  #Setup data frame
  x = seq_along(vec1)
  df <- data.frame(x, vec1, vec2)
  
  ggplot(df, aes(x)) +                    # basic graphical object
    geom_point(aes(y = vec1), colour="red") +  # first layer
    geom_point(aes(y = vec2), colour="green") +  # second layer
    theme_bw() +
    xlab("R0") + ylab("Original (green), R packages (red)") + 
    ggtitle("Log(acceptance_prob); Original (green), R packages(red)") 
  
}

#Plot points comparison
#jk = 1:10
#plot_points_comparison(jk, jk^2)


#Plots
plot_mcmc_super_spreading <- function(sim_data, mcmc_vector1, mcmc_vector2, mcmc_vector3, alpha_true, betaX, gammaX, file_name, folder_dir_ad) {
  
  #Folder save
  ifelse(!dir.exists(file.path(folder_dir_ad)), dir.create(file.path(folder_dir_ad)), FALSE)
  
  pdf(paste(folder_dir_ad, "/", file_name, alpha_true, ".pdf", sep=""))
  
  #i. Epidemic data
  plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')
  
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector1, ylab = 'alpha', main = paste("MCMC Results. Alpha - Super spreading model, true alpha = ", alpha_true))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  alpha_mean = cumsum(mcmc_vector1)/seq_along(mcmc_vector1)
  plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alpha_true))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector1, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'alpha', ylab = 'Density', 
               main = 'Empirical density of alpha - MCMC chain')
  print(hist3)
  
  #**************************************************
  #2. beta
  
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector2, ylab = 'beta', main = paste("MCMC Results; beta - Super spreading model, true beta = ", betaX))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  beta_mean = cumsum(mcmc_vector2)/seq_along(mcmc_vector2)
  plot2 = plot(seq_along(beta_mean), beta_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, true beta = ", betaX))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector2, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'alpha', ylab = 'Density', 
               main = 'Empirical density of beta - MCMC chain')
  print(hist3)
  
  #**************************************************
  #*gamma
  
  #i. MCMC chain
  plot1 = ts.plot(mcmc_vector3, ylab = 'gamma', main = paste("MCMC Results; gamma - Super spreading model, true gamma = ", gammaX))
  print(plot1)
  
  #ii. Mean
  #Plot mean
  gamma_mean = cumsum(mcmc_vector3)/seq_along(mcmc_vector3)
  plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True gamma = ", gammaX))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector3, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'gamma', ylab = 'Density', 
               main = 'Empirical density of gamma - MCMC chain')
  print(hist3)
  
  dev.off()
  
}

#Plots
plot_mcmc_super_spreading_to_screen <- function(sim_data, mcmc_vector1, mcmc_vector2, mcmc_vector3, alpha_true, betaX, gammaX, file_name, folder_dir_ad) {
  
  #Folder save
  #ifelse(!dir.exists(file.path(folder_dir_ad)), dir.create(file.path(folder_dir_ad)), FALSE)
  #pdf(paste(folder_dir_ad, "/", file_name, alpha_true, ".pdf", sep=""))
  
  #i. Epidemic data
  plot.ts(sim_data, ylab = 'Daily Infections count', main = 'Daily Infections count')
  
  #i. MCMC chain
  plot.ts(mcmc_vector1, ylab = 'alpha', main = paste("MCMC Results. Alpha - Super spreading model, true alpha = ", alpha_true))
  
  #ii. Mean
  #Plot mean
  alpha_mean = cumsum(mcmc_vector1)/seq_along(mcmc_vector1)
  plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alpha_true))
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector1, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'alpha', ylab = 'Density', 
               main = 'Empirical density of alpha - MCMC chain')
  print(hist3)
  
  #**************************************************
  #2. beta
  
  #i. MCMC chain
  plot.ts(mcmc_vector2, ylab = 'beta', main = paste("MCMC Results; beta - Super spreading model, true beta = ", betaX))
  #print(plot1)
  
  #ii. Mean
  #Plot mean
  beta_mean = cumsum(mcmc_vector2)/seq_along(mcmc_vector2)
  plot2 = plot(seq_along(beta_mean), beta_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, true beta = ", betaX))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector2, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'alpha', ylab = 'Density', 
               main = 'Empirical density of beta - MCMC chain')
  print(hist3)
  
  #**************************************************
  #*gamma
  
  #i. MCMC chain
  plot.ts(mcmc_vector3, ylab = 'gamma', main = paste("MCMC Results; gamma - Super spreading model, true gamma = ", gammaX))
  #print(plot1)
  
  #ii. Mean
  #Plot mean
  gamma_mean = cumsum(mcmc_vector3)/seq_along(mcmc_vector3)
  plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True gamma = ", gammaX))
  print(plot2)
  
  #Histogram
  #hist1 = hist(mcmc_vector, prob = TRUE)
  #print(hist1)
  
  #Hist
  hist2 <- hist(mcmc_vector3, breaks = 80)
  hist2$counts <- hist2$counts/sum(hist2$counts)
  hist3 = plot(hist2, xlab = 'gamma', ylab = 'Density', 
               main = 'Empirical density of gamma - MCMC chain')
  print(hist3)
  
  dev.off()
  
}

#Plotting + Ggplot
#Setup data frame
ggplot_fn_temp <- function() {
  
  x = seq_along(jk)
  df <- data.frame(x, jk^2, jk^2.5)
  ggplot(df, aes(x)) +                    # basic graphical object
    geom_point(aes(y = jk^2), colour="red") +  # first layer
    geom_point(aes(y = jk^2.5), colour="green") +  # second layer
    theme_bw() +
    xlab("R0") + ylab("True R0 (black), MCMC sample (red)") + 
    ggtitle("True R0 vs Mean of MCMC sample")
  
}

#Plot main title in grid plot

line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}


#Brainstorm

#********************
#Dataframe
df <- data.frame(
  ab = c(1,5,7,8,9,5,11,3,1,4),
  bc = c(2,3,4,2,6,8,9,10,2,8)
)
df

#Save
write.csv(df, file = paste('df_p_vals_', iter, '.csv'), iter, row.names = FALSE)
iter = 3
write.csv(df, file = paste("df_p_vals_.csv", iter), row.names = FALSE)

#Add row
df_results[nrow(df_results) + 1,] = list(4, 7)
df_results

#Apply to a row
apply(df, 2, min)

#Appply function to a row

#Get p values from summarys stats
check_p_val <- function(sim_data, vec2) {
  
  true_sum_inf = sum(sim_data)
  print(true_sum_inf)
  
  #P value
  lt = length(which(vec2 < true_sum_inf))
  print(lt)
  gt = length(which(vec2 > true_sum_inf))
  print(gt)
  min_val = min(lt, gt)
  pvalue = min_val/length(vec2)
  
  print((pvalue))
  
}

apply(df, 2, FUN = function(vec2) check_p_val(sim_data, vec2))

apply(df_results,2,min)

#*********************************************
#Apply function #2

#Get p values from summarys stats
compare_row_vals <- function(column) {
  
  #Final val
  last_el = column[length(column)]
  cat('last element = ', last_el)
  #P value
  lt = length(which(column < last_el))
  gt = length(which(column > last_el))
  min_val = min(lt, gt)
  pvalue = min_val/length(column)
  
  #print(pvalue)
  pvalue
  
}

apply(df, 2, FUN = function(vec) compare_row_vals(vec))


#*************
#*Model Criticism
#Model Criticism Function
plot_model_criticism <- function(mcmc_params, sim_data, max_sum_val) { 
  
  #Plot Model Criticism
  vec_mod_crit = mcmc_params[9]
  vec_mod_crit = unlist(vec_mod_crit)
  true_sum_inf = sum(sim_data)
  
  #P value
  lt = length(which(vec_mod_crit < true_sum_inf))
  print(lt)
  gt = length(which(vec_mod_crit > true_sum_inf))
  print(gt)
  min_val = min(lt, gt)
  pvalue = min_val/length(vec_mod_crit)
  
  #Check
  if (lt < gt){
    flag = 'lt (<)'
  } else if (gt < lt){
    flag = 'gt (>)'
  }
  
  #Histogram
  hist(vec_mod_crit[vec_mod_crit < max_sum_val], breaks = 100, #freq = FALSE, 
       #xlim = c(xmin, xmax),
       xlab = paste('Sum of Infecteds <', max_sum_val), ylab = 'Density',
       main = paste('Model criticism, true R0 = ', true_r0, '.',
                    'P value', flag, '=', pvalue),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_sum_inf, col = 'red', lwd = 2)
  
}



#****

#Get p values - comparing  summary stat columns to true value 
get_p_value <- function(column) {
  
  #Final val
  last_el = column[length(column)] #True value 
  cat('last element = ', last_el)
  #P value
  lt = length(which(column < last_el))
  gt = length(which(column > last_el))
  min_val = min(lt, gt)
  pvalue = min_val/length(column)
  pvalue = pvalue/2
  
  #Return p value 
  pvalue
  
}

apply(df, 2, FUN = function(vec) get_p_value(vec))


#Plot Epidemics
par(mfrow = c(4,3))

# #*Implement
num_days = 50
#lambda params
shape_gamma = 6
scale_gamma = 1
#params
alphaX = 0.8 #Without ss event, ~r0.
betaX = 0.1
gammaX = 10
true_r0 = alphaX + betaX*gammaX
true_r0
#Epidemic data
sim_data = simulate_branching_ss(num_days, shape_gamma, scale_gamma, alphaX, betaX, gammaX)
plot.ts(sim_data, ylab = 'Daily Infections', col = 'blue', lwd = 2,
        main = paste('Super Spreading Events Infs, R0 = ', true_r0))



#*Implement
num_days = 50
#lambda params
shape_gamma = 6
scale_gamma = 1
#params
aX = 0.8 #1.1 #Without ss event, ~r0.
bX = 0.1 #0.2
ss_mult = 10 #8
#Epidemic data
sim_data2 = simulation_super_spreaders(num_days, shape_gamma, scale_gamma, aX, bX, ss_mult)
plot.ts(sim_data2, ylab = 'Daily Infections count', col = 'red',
        main = 'Super Spreaders Infections')

#Orig
r0 = 1.8
sim_data0 = simulate_branching(num_days, r0, shape_gamma, scale_gamma)
plot.ts(sim_data0, ylab = 'Daily Infections', lwd = 2,
        main = paste('Regular Infections, R0 = ', r0))
