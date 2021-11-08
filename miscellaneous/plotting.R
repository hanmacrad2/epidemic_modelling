#Plotting brainstorming

#Plots
#alpha
plot.ts(alpha_mcmc, ylab = 'alpha', main = paste("MCMC Super spreading model, simulated alpha = ", alphaX))
#alpha mean
alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
plot2 = plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alphaX))
print(plot2)

#beta
plot.ts(beta_mcmc, ylab = 'beta', main = paste("MCMC Super spreading model, simulated beta = ", betaX))
#beta mean
beta_mean = cumsum(beta_mcmc)/seq_along(beta_mcmc)
plot2 = plot(seq_along(beta_mean), beta_mean, xlab = 'Time', ylab = 'beta', main = paste("Mean of beta MCMC chain, True beta = ",betaX))
print(plot2)

#gamma
plot.ts(gamma_mcmc,  ylab = 'gamma', main = paste("MCMC Super spreading model, simulated gamma = ", gammaX))
#gamma Mean
gamma_mean = cumsum(gamma_mcmc)/seq_along(gamma_mcmc)
plot2 = plot(seq_along(gamma_mean), gamma_mean, xlab = 'Time', ylab = 'gamma', main = paste("Mean of gamma MCMC chain, True gamma = ",gammaX))
print(plot2)

#r0
plot.ts(r0_total_mcmc,  ylab = 'r0', main = paste("R0 total - MCMC Super spreading model, true total r0 = ", true_tot_r0))
#r0 mean
r0_tot_mean = cumsum(r0_total_mcmc)/seq_along(r0_total_mcmc)
plot2 = plot(seq_along(r0_tot_mean), r0_tot_mean, xlab = 'Time', ylab = 'r0 total', main = paste("Mean of R0 total MCMC chain, True R0 total = ", true_tot_r0))
print(plot2)

#r0 hist
#hist(r0_total_mcmc, prob = TRUE, breaks = 80,main = paste("Histogram of R0_total MCMC samples, True R0 total = ", true_tot_r0))
#Hist - density
hist2 <- hist(r0_total_mcmc, breaks = 80)
hist2$counts <- hist2$counts/sum(hist2$counts)
hist3 = plot(hist2, xlab = 'R0 total', ylab = 'Density', 
             main = 'Empirical density of R0 total - MCMC samples')
print(hist3)

#*******************
#PLOT ALL ON ONE PLOT
par(mfrow=c(2,3))

#Sim data
flag_dist = 'Negative Binomial distribution,'
plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count', main = paste("Daily Infections count, ", flag_dist, "r0 = ", true_tot_r0),
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

#Plots
#alpha
plot.ts(alpha_mcmc, xlab = 'Time', ylab = 'alpha', main = paste("MCMC Super spreading model, simulated alpha = ", alphaX),
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
#alpha mean
alpha_mean = cumsum(alpha_mcmc)/seq_along(alpha_mcmc)
plot(seq_along(alpha_mean), alpha_mean, xlab = 'Time', ylab = 'alpha', main = paste("Mean of alpha MCMC chain, True alpha = ",alphaX),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

#r0
plot.ts(r0_total_mcmc,  xlab = 'Time', ylab = 'r0', main = paste("R0 total - MCMC Super spreading model, true total r0 = ", true_tot_r0),
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
#r0 mean
r0_tot_mean = cumsum(r0_total_mcmc)/seq_along(r0_total_mcmc)
plot2 = plot(seq_along(r0_tot_mean), r0_tot_mean, xlab = 'Time', ylab = 'r0 total', main = paste("Mean of R0 total MCMC chain, True R0 total = ", true_tot_r0),
             cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
print(plot2)

#r0 hist
#hist(r0_total_mcmc, prob = TRUE, breaks = 80,main = paste("Histogram of R0_total MCMC samples, True R0 total = ", true_tot_r0))
#Hist - density

hist(r0_total_mcmc, freq = FALSE, xlab = 'R0 total', ylab = 'Density', 
     main = 'Empirical density of R0 total - MCMC samples',
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

#hist2 <- hist(r0_total_mcmc, breaks = 80)
#hist2$counts <- hist2$counts/sum(hist2$counts)
#hist3 = plot(hist2, xlab = 'R0 total', ylab = 'Density', 
#             main = 'Empirical density of R0 total - MCMC samples')
#print(hist3)
#hist()


#*************************************************************************
#* Densities ****

#Exponential
x <- seq(from = 0, to = 10, by = 1)
exp1 = rexp(n, 1)
plot.ts(exp1)
exp1 = dexp(x, 1)
plot.ts(exp1)

#Prior/Posterior/Likelihood
x <- seq(from = -90, to = 90, by = 1)
data <- dnorm(x, mean = 30, sd = 10)
prior <- dnorm(x, mean = 10, sd = 5)
posterior <- 0.5 * dnorm(x, mean = 10, sd = 5) + 0.5 * dnorm(x, mean = 30, sd = 10)

plot(x, prior, type = "l", col = "red")
lines(x, posterior, type = "l", col = "green")
lines(x, data , type = "l", col = "blue")