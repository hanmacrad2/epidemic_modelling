#Results - Dispaly results of MCMC

#Final MCMC mean vs R0 mean

#Results
vec_r0 = c(0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5.0, 8.0, 10.0)
vec_final_mean = c(0.58, 0.52, 1.01, 0.82, 1.1, 1.32, 0.925, 2.66, 2.89, 3.39, 2.6, 3.65, 3.41, 4.31, 3.725, 4.66, 7.5, 9.06) 

#Dataframe for plot
group = c(rep.int(1, length(vec_r0)), rep.int(2, length(vec_r0)))
r0_results = c(vec_r0, vec_final_mean)
vec_r02 = c(vec_r0, vec_r0)

df_results  = data.frame(group, vec_r02, r0_results)

ggplot(df_results, aes(x = vec_r02, y = r0_results, color = group) ) + #3.5 == 3.41, 3 = 3.65, 
  geom_point(shape=21, fill = group, size=3 ) +
  theme_bw() +
  xlab("R0") + ylab("True R0 (black), MCMC sample (red)") + 
  ggtitle("True R0 vs Mean of MCMC sample") 
 # scale_color_identity(name = 'Data', label = c('True R0', 'MCMC'),
                       #guide = guide_legend())
  #labs(colour = group)
  #geom_smooth(method = "lm", se = FALSE)


#Quantile plot
r_quants = c(c(3.174966, 3.726409)) #3.5

#*********************************
#2 Acceptance Rate

ggplot(df_acc_rate, aes(x=vec_r0, y=vec_rate)) +
  geom_line( color="grey") +
  ylim(0, 100) +
  theme_bw() + 
  xlab("R0") + ylab("Acceptance Rate %") + 
  geom_point(shape=21, color="black", fill="black", size=6) +
  ggtitle("Acceptance rate % for varying R0") 