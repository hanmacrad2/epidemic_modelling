#*##############################################
#*******************************************
#*
#* GRID PLOT SUPER-SPREADING MODELS:
#* 
#* FUNCTION TO PLOT 4x4 DASHBOARD OF MCMC RESULTS FOR SUPER SPREADING MODELs
#* 
#******************************************
PLOT_MCMC_GRID <- function(sim_data, mcmc_output,
                           mcmc_inputs = list(n_mcmc = n_mcmc,
                                              model_params = model_params,
                                              mod_par_names = c('alpha', 'beta', 'gamma'),
                                              sigma = sigma,
                                              model_typeX = 'SSE',
                                              total_time = 0, seed_count = 1,
                                              x0 = 1),
                           priors_list = list(a_prior = c(1, 0), b_prior = c(10, 1/100), b_prior_exp = c(1,0),
                                              c_prior = c(10, 1), c_prior_exp = c(0.1,0)),
                           FLAGS_LIST = list(DATA_AUG = TRUE, BC_TRANSFORM = TRUE,
                                             PRIOR = TRUE, JOINT = TRUE,
                                             B_PRIOR_GAMMA = TRUE, C_PRIOR_GAMMA = TRUE,
                                             RJMCMC = FALSE)) { 
  #PLOT
  #plot.new()
  par(mfrow=c(4,4))
  
  #Extract params
  m1_mcmc = mcmc_output[1]; m1_mcmc = unlist(m1_mcmc)
  m2_mcmc = mcmc_output[2]; m2_mcmc = unlist(m2_mcmc)
  m3_mcmc = mcmc_output[3]; m3_mcmc = unlist(m3_mcmc)
  r0_mcmc = mcmc_output[4]; r0_mcmc = unlist(r0_mcmc)
 
  #Cumulative means + param sample limits
  #r0
  r0_mean = cumsum(r0_mcmc)/seq_along(r0_mcmc)
  r0_lim = max(true_r0, max(r0_mcmc))
  r0_lim2 = max(true_r0, r0_mean)
  
  #m1
  m1_mean = cumsum(m1_mcmc)/seq_along(m1_mcmc)
  print(mcmc_inputs$model_params$m1[[1]])
  print(max(m1_mcmc))
  a_lim =  max(mcmc_inputs$model_params$m1[[1]], max(m1_mcmc))
  a_lim2 =  max(mcmc_inputs$model_params$m1[[1]], m1_mean)
  
  #m2
  m2_mean = cumsum(m2_mcmc)/seq_along(m2_mcmc)
  b_lim = max(mcmc_inputs$model_params$m2[[1]], max(m2_mcmc))
  b_lim2 = max(mcmc_inputs$model_params$m2[[1]], m2_mean)
  
  #m3
  m3_mean = cumsum(m3_mcmc)/seq_along(m3_mcmc)
  m3_lim =  max(mcmc_inputs$model_params$m3[[1]], max(m3_mcmc))
  m3_lim2 =  max(mcmc_inputs$model_params$m3[[1]], m3_mean) 
  
  #Priors
  if (FLAGS_LIST$B_PRIOR_GAMMA) {
    m2_prior = paste0('m3(', priors_list$b_prior[1], ', ', priors_list$b_prior[2], ')')
  } else {
    m2_prior = paste0('exp(', priors_list$b_prior[1], ')')
  }
  
  if (FLAGS_LIST$C_PRIOR_GAMMA) {
    m3_prior = paste0('1 + Ga(',   priors_list$c_prior[1], ', ',  priors_list$c_prior[2], ')')
  } else {
    m3_prior = paste0('1 + exp(',   priors_list$c_prior_exp, ')')
  }
  
  m1_prior =  paste0('exp(', priors_list$a_prior[1], ')')
  
  #***********
  #* Plots *
  
  #i.Infections
  if(!FLAGS_LIST$DATA_AUG) inf_tite = paste0(mcmc_inputs$seed_count, ', ', mcmc_inputs$model_typeX, " Data, r0 = ", true_r0) # 'Day Infts, '
  else inf_tite = paste0(mcmc_inputs$seed_count, ', ', mcmc_inputs$model_typeX, " Data, r0 = ", true_r0, ", + Data Aug")
  
  #i.Infections
  plot.ts(sim_data, xlab = 'Time', ylab = 'Daily Infections count',
          main = inf_tite, 
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #ii. MCMC Trace Plots 
  plot.ts(m1_mcmc, ylab = mcmc_inputs$mod_par_names[1], ylim=c(0, a_lim),
          main = paste("MCMC", mcmc_inputs$model_typeX, ":", mcmc_inputs$mod_par_names[1], "prior:", m1_prior),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = mcmc_inputs$model_params$m1[[1]], col = 'red', lwd = 2) #True = green
  
  plot.ts(m2_mcmc, ylab = 'm2', ylim=c(0, b_lim), 
          main = paste("MCMC", mcmc_inputs$model_typeX, ":", mcmc_inputs$mod_par_names[2], "prior:", m2_prior),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = mcmc_inputs$model_params$m2[[1]], col = 'blue', lwd = 2) #True = green
  
  plot.ts(m3_mcmc,  ylab = 'm3', ylim=c(0,m3_lim),
          main = paste("MCMC", mcmc_inputs$model_typeX, ":", mcmc_inputs$mod_par_names[3], "prior:", m3_prior),
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = mcmc_inputs$model_params$m3[[1]], col = 'green', lwd = 2) #True = green
  
  #plot.ts(r0_mcmc,  ylab = 'r0', main = paste("MCMC SS Events, true r0 = ", r0_true))
  
  #iii. Cumulative mean plots
  #r0 Mean
  plot(seq_along(r0_mean), r0_mean,
       ylim=c(0, r0_lim),
       xlab = 'Time', ylab = 'R0', main = paste("R0 MCMC Mean, True R0 = ", true_r0),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #print(plot2)
  abline(h = true_r0, col = 'orange', lwd = 2)
  
  #m1 mean
  plot(seq_along(m1_mean), m1_mean,
       ylim=c(0, a_lim),
       xlab = 'Time', ylab =  mcmc_inputs$mod_par_names[1],
       main = paste(mcmc_inputs$mod_par_names[1], "MCMC mean, True", mcmc_inputs$mod_par_names[1], "=", mcmc_inputs$model_params$m1[[1]]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #print(plot2)
  abline(h = mcmc_inputs$model_params$m1[[1]], col = 'red', lwd = 2)
  
  #m2 mean
  plot(seq_along(m2_mean), m2_mean,
       ylim=c(0, b_lim),
       xlab = 'Time', ylab = 'm2',
       main = paste(mcmc_inputs$mod_par_names[2], "MCMC mean, True", mcmc_inputs$mod_par_names[2], "=", mcmc_inputs$model_params$m2[[1]]),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  #print(plot2)
  abline(h = mcmc_inputs$model_params$m2[[1]], col = 'blue', lwd = 2)
  
  #m3 Mean
  plot(seq_along(m3_mean), m3_mean,
       xlab = 'Time', ylab = 'm3', 
       main = paste(mcmc_inputs$mod_par_names[3], "MCMC mean, True", mcmc_inputs$mod_par_names[3], "=", mcmc_inputs$model_params$m3[[1]]),
       ylim=c(0, m3_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h = mcmc_inputs$model_params$m3[[1]], col = 'green', lwd = 2)
  
  #iv. Param Histograms (Plots 9,11,12)
  hist(r0_mcmc, freq = FALSE, breaks = 100,
       xlab = 'R0 total', #ylab = 'Density', 
       main = paste('R0 total MCMC samples'),
       xlim=c(0, r0_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = true_r0, col = 'orange', lwd = 2)
  
  # #v. m2 vs m3
  # plot(m2_mcmc, m3_mcmc,
  #      xlab = 'm2', ylab = 'm3', main = 'm2 vs m3',
  #      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  
  #Hist m1 
  hist(m1_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_inputs$mod_par_names[1], #ylab = 'Density', 
       main = paste(mcmc_inputs$mod_par_names[1], ", True", mcmc_inputs$mod_par_names[1], "=", mcmc_inputs$model_params$m1[[1]]), 
       xlim=c(0, a_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = mcmc_inputs$model_params$m1[[1]], col = 'red', lwd = 2)
  
  #Hist m2 
  hist(m2_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_inputs$mod_par_names[2], #ylab = 'Density', 
       main = paste(mcmc_inputs$mod_par_names[2], ", True", mcmc_inputs$mod_par_names[2], "=", mcmc_inputs$model_params$m2[[1]]), 
       xlim=c(0, b_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = mcmc_inputs$model_params$m2[[1]], col = 'blue', lwd = 2)
  
  #Hist m3 
  hist(m3_mcmc, freq = FALSE, breaks = 100,
       xlab = mcmc_inputs$mod_par_names[3], #ylab = 'Density', 
       main = paste(mcmc_inputs$mod_par_names[3], ", True", mcmc_inputs$mod_par_names[3], "=", mcmc_inputs$model_params$m3[[1]]),
       xlim=c(0, m3_lim),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(v = mcmc_inputs$model_params$m3[[1]], col = 'green', lwd = 2)
  
  #Final Mean Stats
  data_10_pc = 0.5*mcmc_inputs$n_mcmc #50%
  m1_mean_tail = round(mean(m1_mcmc[mcmc_inputs$n_mcmc - data_10_pc:mcmc_inputs$n_mcmc]), 2) 
  m2_mean_tail = round(mean(m2_mcmc[mcmc_inputs$n_mcmc - data_10_pc:mcmc_inputs$n_mcmc]), 2)
  m3_mean_tail = round(mean(m3_mcmc[mcmc_inputs$n_mcmc - data_10_pc:mcmc_inputs$n_mcmc]), 2)
  m4_mean_tail = round(mean(r0_mcmc[mcmc_inputs$n_mcmc - data_10_pc:mcmc_inputs$n_mcmc]), 2)
  
  #FLAGS_LIST$JOINT distrbutions
  if (FLAGS_LIST$JOINT){
    
    #v. r0 vs m2
    plot(m2_mcmc, r0_mcmc,
         xlab = mcmc_inputs$mod_par_names[2], ylab = 'R0', main = paste(mcmc_inputs$mod_par_names[2], 'vs R0'),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
         cex.main = 0.5)
    
    #v. m1 vs m2
    plot(m1_mcmc, m2_mcmc,
         xlab = mcmc_inputs$mod_par_names[1], ylab = mcmc_inputs$mod_par_names[2], main = paste(mcmc_inputs$mod_par_names[1], 'vs', mcmc_inputs$mod_par_names[2]),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
         cex = 0.5)
    
    #v. m1 vs m3
    plot(m1_mcmc, m3_mcmc,
         xlab = mcmc_inputs$mod_par_names[1], ylab = mcmc_inputs$mod_par_names[3], main = paste(mcmc_inputs$mod_par_names[1], 'vs', mcmc_inputs$mod_par_names[3]),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
         cex = 0.5)
    
    #v. m2 vs m3
    plot(m2_mcmc, m3_mcmc,
         xlab = mcmc_inputs$mod_par_names[2], ylab = mcmc_inputs$mod_par_names[3], main = paste(mcmc_inputs$mod_par_names[2], 'vs', mcmc_inputs$mod_par_names[3]),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
         cex = 0.5)
  }
  
  #*****************
  #RJMCMC
  if (FLAGS_LIST$RJMCMC){
    
    #m1S
    par(mfrow=c(3,1))
    
    #HIST m1 
    hist(m1_mcmc, freq = FALSE, breaks = 100,
         xlab = mcmc_inputs$mod_par_names[1], #ylab = 'Density', 
         main = paste("m1, True m1 = ", mcmc_inputs$model_params$m1[[1]]), 
         xlim=c(0, a_lim),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = mcmc_inputs$model_params$m1[[1]], col = 'red', lwd = 2)
    
    #m1 SSE
    ind_b = which(m2_mcmc > 0)
    m1_sse = m1_mcmc[ind_b]
    
    #Hist
    hist(m1_sse, freq = FALSE, breaks = 100,
         xlab = 'm1_sse', #ylab = 'Density', 
         main = "m1_sse", 
         xlim=c(0, a_lim),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    abline(v = mcmc_inputs$model_params$m1[[1]], col = 'red', lwd = 2)
    
    #m1 BASE
    ind_b = which(m2_mcmc == 0)
    
    if (length(ind_b > 2)){
      
      m1_base = m1_mcmc[ind_b]
      #Hist
      hist(m1_base, freq = FALSE, breaks = 100,
           xlab = 'm1_base', #ylab = 'Density', 
           main = "m1_base",
           xlim=c(0, a_lim),
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
      abline(v = mcmc_inputs$model_params$m1[[1]], col = 'red', lwd = 2)
    }
    
    #RESULTS DF
    df_results <- data.frame(
      rep = mcmc_inputs$seed_count,
      n_mcmc = mcmc_inputs$n_mcmc,
      m1 = mcmc_inputs$model_params$m1[[1]],
      a_mc = a_mcmc_mean,
      m2 = mcmc_inputs$model_params$m2[[1]],
      b_mc = b_mcmc_mean,
      m3 = mcmc_inputs$model_params$m3[[1]],
      c_mc = c_mcmc_mean,
      R0 = true_r0, 
      R0_mc = r0_mcmc_mean,
      accept_rate_m1 = round(mcmc_output[[5]],2),
      a_rte_m2 = round(mcmc_output[[6]], 2),
      n_accept_m2 = mcmc_output[[15]],
      a_rte_m3 = round(mcmc_output[[7]],2),
      n_accept_m3 = mcmc_output[[16]],
      a_rte_m2_m3 = round(mcmc_output[[8]],2),
      n_accept_2_3 = mcmc_output[[17]],
      n_accept_rj0 = mcmc_output[[11]],
      n_reject_rj0 = mcmc_output[[13]],
      a_rte_rj0 = round(mcmc_output[[9]],2),
      n_accept_rj1 = mcmc_output[[12]],
      n_reject_rj1 = mcmc_output[[14]],
      a_rte_rj1 = round(mcmc_output[[10]],2),
      m2_pc0 = mcmc_output[[18]],
      m2_pc_non_0 = 1- mcmc_output[[18]],
      bf = mcmc_output[[19]])
    #tot_time = mcmc_inputs$total_time)
    
  } else if (FLAGS_LIST$DATA_AUG) {
    
    df_results <- data.frame(
      rep = mcmc_inputs$seed_count,
      n_mcmc = mcmc_inputs$n_mcmc,
      m1 = mcmc_inputs$model_params$m1[[1]],
      m1_mc = m1_mean_tail,
      m2 = mcmc_inputs$model_params$m2[[1]],
      m2_mc = m2_mean_tail,
      m3 = mcmc_inputs$model_params$m3[[1]],
      m3_mc = m3_mean_tail,
      R0 = true_r0, 
      R0_mc = m4_mean_tail,
      accept_rate_m1 = round(mcmc_output$list_accept_rates$accept_rate1, 2),
      a_rte_m2 = round(mcmc_output$list_accept_rates$accept_rate2, 2),
      a_rte_m3 = round(mcmc_output$list_accept_rates$accept_rate3, 2),
      a_rte_m2_m3 = round(mcmc_output$list_accept_rates$accept_rate4, 2),
      a_rte_d_aug = round(mcmc_output$list_accept_rates$accept_rate5, 2),
      a_es = effectiveSize(as.mcmc(m1_mcmc))[[1]],
      b_es = effectiveSize(as.mcmc(m2_mcmc))[[1]],
      c_es = effectiveSize(as.mcmc(m3_mcmc))[[1]],
      d_es = effectiveSize(as.mcmc(r0_mcmc))[[1]],
      time_elap = format(mcmc_output$time_elap, format = "%H:%M:%S")[1])
    
  } else {
    df_results <- data.frame(
      rep = mcmc_inputs$seed_count,
      n_mcmc = mcmc_inputs$n_mcmc,
      m1 = mcmc_inputs$model_params$m1[[1]],
      m1_mc = a_mcmc_mean,
      m2 = mcmc_inputs$model_params$m2[[1]],
      m2_mc = b_mcmc_mean,
      m3 = mcmc_inputs$model_params$m3[[1]],
      m3_mc = c_mcmc_mean,
      R0 = true_r0, 
      R0_mc = r0_mcmc_mean,
      accept_rate_m1 = round(mcmc_output$list_accept_rates$accept_rate1, 2),
      a_rte_m2 = round(mcmc_output$list_accept_rates$accept_rate2, 2),
      a_rte_m3 = round(mcmc_output$list_accept_rates$accept_rate3, 2),
      a_rte_m2_m3 = round(mcmc_output$list_accept_rates$accept_rate4, 2),
      a_es = effectiveSize(as.mcmc(a)),
      b_es = effectiveSize(as.mcmc(b)),
      c_es = effectiveSize(as.mcmc(c)),
      d_es = effectiveSize(as.mcmc(d)),
      tot_time = mcmc_inputs$total_time)
  }
  
  print(df_results)
  
}
