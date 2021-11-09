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