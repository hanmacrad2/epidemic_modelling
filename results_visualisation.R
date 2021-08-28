#Results - Dispaly results of MCMC
library(ggplot2)
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


#************************************************
#Quantile plot
vec_final_mean = c(0.99, 0.61, 0.50, 0.60, 1.09, 1.01, 1.75, 2.00, 2.41, 2.37, 2.69, 2.92, 3.23, 4.31, 3.725, 4.66, 7.5, 9.06)
vec_q1 = c(0.51, 0.12, 0.06, 0.12, 0.52, 0.46, 1.45, 1.67, 1.96, 2.03, 2.34, 2.67, 3.01, 3.5, 3.5, 4.35, 7.0, 8.5)
vec_q2 = c(1.62, 1.48, 1.38, 1.44, 1.86, 1.77, 2.06, 2.35, 2.92, 2.75, 3.07, 3.18, 3.46, 4.5,  3.9, 4.85, 8.0, 9.5)

#Dataframe for plot
group = c(rep.int(1, length(vec_r0)), rep.int(2, length(vec_r0)))
r0_results = c(vec_r0, vec_final_mean)
vec_r02 = c(vec_r0, vec_r0)
q1 = c(vec_final_mean, vec_q1) #0, length(vec_r0)) 
q2 = c(vec_final_mean, vec_q2) #0, length(vec_r0))
col_group = c(rep('green', length(vec_r0)), rep('black', length(vec_r0)))

col_group = c(rep(A = "#333BFF", length(vec_r0)), rep( B = "#9633FF", length(vec_r0)))

df_results  = data.frame(group, vec_r02, r0_results, q1, q2, col_group)
df_results

r_quants = c(c(3.174966, 3.726409)) #3.5

#Quantile df_ad_results_formI
ggplot(df_results, aes(x = vec_r02, y = r0_results, color = group) ) + #3.5 == 3.41, 3 = 3.65, 
  geom_point(shape=21, size=3 ) +
  theme_bw() +
  geom_errorbar(aes(ymin=q1, ymax= q2), width=.2,
                position=position_dodge(0.05)) +
  xlab("R0") + ylab("True R0 vs Mean of MCMC sample)") + 
  ggtitle("True R0 (black), MCMC sample (red)") #+
  #scale_fill_manual(values=col_group)

group.colors <- c(A = "#333BFF", B = "#CC6600") #, C ="#9633FF", D = "#E2FF33", E = "#E3DB71")


#Attempt 2
ggplot(df_results, aes(x = vec_r02, y = r0_results, color = group) ) + #3.5 == 3.41, 3 = 3.65, 
  geom_point(aes(fill=factor(group)), shape=21, size=4 ) +
  theme_bw() +
  geom_errorbar(aes(ymin=q1, ymax= q2), width=.2,
                position=position_dodge(0.05)) +
  xlab("R0") + ylab("True R0 vs Mean of MCMC sample") + 
  ggtitle("True R0 (red), MCMC sample + 95% quantile (blue)") #+
#scale_fill_manual(values=col_group)



ggplot(df_results, aes(x = vec_r02, y = r0_results, color = color=factor(group)) ) + #3.5 == 3.41, 3 = 3.65, 
  geom_errorbar(aes(ymin=q1, ymax= q2, color = factor(group)), width=.2, position=position_dodge(0.05)) +
  scale_color_manual("group", breaks=c(1,2),values=c("#0072B2", "#009E73"))+ ##E69F00" "#D55E00"
  geom_point(aes(fill=factor(group)),size=5, shape=21)+
  scale_fill_manual("group",breaks=c(1,2),values=c("#0072B2", "#009E73"))+
  theme_bw() +
  xlab("R0") + ylab("True R0 vs Mean of MCMC sample)") + 
  ggtitle("True R0 (black), MCMC sample (red)") #+
#scale_fill_manual(values=col_group)

#
#scale_color_manual("group", breaks=c(1,2),values=c("#0072B2", "#009E73"))+ "#E69F00", "#D55E00"
#  geom_point(aes(fill=factor(group)),size=3, shape=21)+
#  scale_fill_manual("group",breaks=c(1,2),values=c("#0072B2", "#009E73"))+ #, "#E69F00", "#D55E00"
  
#Example
    ggplot(data, aes(a, mean)) +
    geom_point()+
    geom_errorbar(aes(ymax=mean+CI,ymin=mean-CI, color=factor(stress)), width=0.3)+
    scale_color_manual("Stress", breaks=c(1,2,3,4),values=c("#0072B2", "#009E73", "#E69F00", "#D55E00"))+
    geom_point(aes(fill=factor(stress)),size=8, shape=21)+
    scale_fill_manual("Stress",breaks=c(1,2,3,4),values=c("#0072B2", "#009E73", "#E69F00", "#D55E00"))+
    scale_x_continuous("Level A",breaks=c(10,20))+
    ylab(expression("Level B"))+
    theme_bw(17)




#Colour
#Create a custom color scale
library(RColorBrewer)
myColors <- brewer.pal(2,"Set1")
names(myColors) <- levels(dat$grp)
colScale <- scale_colour_manual(name = "grp",values = myColors)

#*********************************
#2 Acceptance Rate

ggplot(df_acc_rate, aes(x=vec_r0, y=vec_rate)) +
  geom_line( color="grey") +
  ylim(0, 100) +
  theme_bw() + 
  xlab("R0") + ylab("Acceptance Rate %") + 
  geom_point(shape=21, color="black", fill="black", size=6) +
  ggtitle("Acceptance rate % for varying R0") 

