## This file is aimed at fitting the data collected by Moral_asso_behavioral_exp7 Using Ex-Gaussian Model
#  Author     Date           History of change
# ============================================
#  Hcp     06/07/2016        Save rawData2 as the processed data
#  hcp     06/07/2016        begin with saved csv data
#  hcp     19/07/2016        only the code for final analysis
#  hcp     20/02/2017        used for the final data
#  hcp     08/04/2017        save to pdf
#  hcp     15/07/2018        update to latest format and data
#  hcp     Jan 09, 2019      delete extra comments and prepared to open

########## preparation #################################################
# define directories
curDir <- dirname(rstudioapi::getSourceEditorContext()$path) #Get the directory ofcurrent script
setwd(curDir)

source('Initial_exgaussian.r')  # run the initilization code
library(retimes)     # load the package for ex-gaussian

########## load and clear data ########################################
df.exG <- utils::read.csv("MS_categ_exG.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE) 

# render the numeric column
cols.num <- c("RT","ACC")
df.exG[cols.num] <- sapply(df.exG[cols.num], as.numeric)

# check the class for each column
sapply(df.exG,class)

rm(cols.num) # remove temporary variables

########## begin analysis#############################################
nSub <- length(unique(df.exG$Subject))             # number of subjects, save for later for loop
nTask <- length(unique(df.exG$Task))               # number of tasks, save for later for loop
nIdentity <- length(unique(df.exG$Identity))       # number of identies, save for later for loop
nMorality <- length(unique(df.exG$Morality))       # number of moral levels, save for later for loop

rowNum <- nSub*nTask*nIdentity*nMorality           # predefine the number of rows for saving data

# define the colnames and rows for the results dataframe
trialNum <- rep(-1,rowNum)
subIdtemp <- rep(-1,rowNum)
taskTemp <- as.character(rep(1,rowNum))
identityTemp <- as.character(rep(1,rowNum))
moralitytemp <- as.character(rep(1,rowNum))
Pmu <- rep(-1,rowNum)
Psigma <- rep (-1,rowNum)
Ptau <- rep(-1,rowNum)
AIC <- rep(-1,rowNum)
BIC <- rep(-1,rowNum)

# create an empty data frame for results parameters
results.params <- data.frame(trialNum=trialNum,
                             subId=subIdtemp,
                             Task=taskTemp,
                             Identity=identityTemp,
                             Morality=moralitytemp,
                             Pmu =Pmu,
                             Psigma=Psigma,
                             Ptau=Ptau,
                             AIC=AIC,
                             BIC=BIC,
                             stringsAsFactors=FALSE)
rm(trialNum,subIdtemp,taskTemp,identityTemp,moralitytemp,Pmu,Psigma,Ptau,AIC,BIC)

# Key part for ex-Gaussian analysis: select each condition for each particiant, and fit the data with ex-Gaussian fucntion
# then, save three parameters' estimation in 'results.params'
index <- 1     
for (ii in 1:nSub){
        for (jj in 1:nTask){
                for (kk in 1:nIdentity) {
                        for (ll in 1: nMorality){
                                dat <- df.exG[df.exG$Subject == RTsubId[ii] &
                                                             df.exG$Task == RTtask[jj] & 
                                                             df.exG$Identity == RTidentity[kk] &
                                                             df.exG$Morality == RTMorality[ll],]
                                if (nrow(dat) < 10) {
                                        results.params$trialNum[index] <- nrow(dat)
                                        results.params$subId[index] <- RTsubId[ii]
                                        results.params$Task[index] <-  RTtask[jj]
                                        results.params$Identity[index] <- RTidentity[kk]
                                        results.params$Morality[index] <- RTMorality[ll]
                                }else{
                                        exGassparams <- timefit(dat$RT,iter = 1000)
                                        results.params$trialNum[index] <- nrow(dat)
                                        results.params$subId[index] <- RTsubId[ii]
                                        results.params$Task[index] <-  RTtask[jj]
                                        results.params$Identity[index] <- RTidentity[kk]
                                        results.params$Morality[index] <- RTMorality[ll]
                                        results.params$Pmu[index] <- exGassparams@par[1]
                                        results.params$Psigma[index] <- exGassparams@par[2]
                                        results.params$Ptau[index] <- exGassparams@par[3]
                                        results.params$AIC[index] <- exGassparams@AIC
                                        results.params$BIC[index] <- exGassparams@BIC
                                }
                                index <- index + 1
                        }
                        
                }
                
        }
}

# remove invalid data in results of params
results.params_valid <- results.params[results.params$trialNum > 36,]
valid_condition <- count(results.params_valid,'subId') # check if there some invalid data 
params1 <- results.params_valid[-which(results.params_valid$subId %in% valid_condition[valid_condition$freq < 8,]$subId),]

# save to csv
write.csv(params1,'ExGaussian_params_2.csv',row.names = FALSE)

# transfer to wide-format
# mu
params_wide_mu <- reshape2::dcast(params1,subId ~ Task + Identity + Morality, value.var = "Pmu")
# sigma
params_wide_sigma <- reshape2::dcast(params1,subId ~ Task + Identity + Morality, value.var = "Psigma")
# tau
params_wide_tau <- reshape2::dcast(params1,subId ~ Task + Identity + Morality, value.var = "Ptau")
# save the data
write.csv(params_wide_mu,'ExGaussian_params_wide_mu.csv',row.names = FALSE)
write.csv(params_wide_sigma,'ExGaussian_params_wide_sigma.csv',row.names = FALSE)
write.csv(params_wide_tau,'ExGaussian_params_wide_tau.csv',row.names = FALSE)

# analysis of params
Pmu_anova    <- ez::ezANOVA(params1,dv = Pmu, wid = subId, within=.(Task,Morality,Identity), type=3)
Psigma_anova <- ez::ezANOVA(params1,dv = Psigma, wid = subId, within=.(Task,Morality,Identity), type=3)
Ptau_anova   <- ez::ezANOVA(params1,dv = Ptau, wid = subId, within=.(Task,Morality,Identity), type=3)

# calculate the summary data for each parameter
Pmu.sum1 <- summarySEwithin(params1,measurevar = 'Pmu',withinvars = c('Task','Identity','Morality'),idvar = 'subId')
Pmu.sum1$Task <- factor(Pmu.sum1$Task, levels=c("Val","Id"))
Pmu.sum1$Identity <- factor(Pmu.sum1$Identity,levels = c('Self','Other'))
Pmu.sum1$Morality <- factor(Pmu.sum1$Morality,levels = c('Good','Bad'))

# analyze the main effect of identity
Pmu.sum2 <- summarySEwithin(params1,measurevar = 'Pmu',withinvars = c('Identity'),idvar = 'subId')

Psigma.sum1 <- summarySEwithin(params1,measurevar = 'Psigma',withinvars = c('Task','Identity','Morality'),idvar = 'subId')
Psigma.sum1$Task <- factor(Psigma.sum1$Task, levels=c("Val","Id"))
Psigma.sum1$Identity <- factor(Psigma.sum1$Identity,levels = c('Self','Other'))
Psigma.sum1$Morality <- factor(Psigma.sum1$Morality,levels = c('Good','Bad'))

# analyzing the interaction between morality and identity (combined different tasks)

params_sigma_noTask <- summarySEwithin(params1,measurevar = 'Psigma',withinvars = c('subId','Identity','Morality'),idvar = 'subId')
params_wide_sigma_noTask <- dcast(params_sigma_noTask,subId ~ Identity + Morality, value.var = "Psigma")
write.csv(params_wide_sigma_noTask,'ExGaussian_params_wide_sigma_noTask.csv',row.names = FALSE)

Psigma_anova_noTask <- reshape2::ezANOVA(params_sigma_noTask,dv = Psigma, wid = subId, within=.(Morality,Identity), type=3)
Psigma_anova_noTask.sum <- summarySEwithin(params1,measurevar = 'Psigma',withinvars = c('Identity','Morality'),idvar = 'subId')
Psigma_anova_noTask.sum$Identity <- factor(Psigma_anova_noTask.sum$Identity,levels = c('Self','Other'))
Psigma_anova_noTask.sum$Morality <- factor(Psigma_anova_noTask.sum$Morality,levels = c('Good','Bad'))

# t tests
Psigma.t.mrl_imm_self <- t.test(params_wide_sigma_noTask$self_good,params_wide_sigma_noTask$self_bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.mrl_imm_self[3]),method = 'bonferroni',4)
params_wide_sigma_noTask$mrl_imm_self <- params_wide_sigma_noTask$self_good - params_wide_sigma_noTask$self_bad
Psigma.t.mrl_imm_self.ci <- bootES(params_wide_sigma_noTask$mrl_imm_self,R = 20000, effect.type = "cohens.d")

Psigma.t.mrl_imm_other <- t.test(params_wide_sigma_noTask$Other_Good,params_wide_sigma_noTask$Other_Bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.mrl_imm_other[3]),method = 'bonferroni',4)
params_wide_sigma_noTask$mrl_imm_other <- params_wide_sigma_noTask$Other_Good - params_wide_sigma_noTask$Other_Bad
Psigma.t.mrl_imm_other.ci <- bootES(params_wide_sigma_noTask$mrl_imm_other,R = 20000, effect.type = "cohens.d")

Psigma.t.slf_oth_mrl <- t.test(params_wide_sigma_noTask$Self_Good,params_wide_sigma_noTask$Other_Good,paired = TRUE)
p.adjust(as.numeric(Psigma.t.slf_oth_mrl[3]),method = 'bonferroni',4)
params_wide_sigma_noTask$mrl_self_other <- params_wide_sigma_noTask$Self_Good - params_wide_sigma_noTask$Other_Good
Psigma.t.mrl_self_other.ci <- bootES(params_wide_sigma_noTask$mrl_self_other,R = 20000, effect.type = "cohens.d")

Psigma.t.slf_oth_imm <- t.test(params_wide_sigma_noTask$Self_Bad,params_wide_sigma_noTask$Other_Bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.slf_oth_imm[3]),method = 'bonferroni',4)
params_wide_sigma_noTask$slf_oth_imm <- params_wide_sigma_noTask$Self_Bad - params_wide_sigma_noTask$Other_Bad
Psigma.t.slf_oth_imm.ci <- bootES(params_wide_sigma_noTask$slf_oth_imm,R = 20000, effect.type = "cohens.d")

Ptau.sum1 <- summarySEwithin(params1,measurevar = 'Ptau',withinvars = c('Task','Identity','Morality'),idvar = 'subId')
Ptau.sum1$Task <- factor(Ptau.sum1$Task, levels=c("Val","Id"))
Ptau.sum1$Identity <- factor(Ptau.sum1$Identity,levels = c('Self','Other'))
Ptau.sum1$Morality <- factor(Ptau.sum1$Morality,levels = c('Good','Bad'))

# library(ggplot2)
p_mu1 <- ggplot(data = Pmu.sum1, aes(y=Pmu,x=Morality,group=Task,shape = Task,fill = Task)) +
        geom_bar(position = position_dodge(),stat = "identity",colour = "black", size=.3) +         # Thinner lines
        #geom_errorbar(aes(ymin = RT-ci,ymax=RT+ci),colour="black",width=.1,position = position_dodge()) +
        #geom_line(position = pd) + geom_point(shape=21,fill="white",position = pd) +
        #geom_point(position = pd) +
        geom_errorbar(aes(ymin = Pmu-se, ymax = Pmu + se),
                      size = 1.2,
                      width = .2,
                      position=position_dodge(.9)) +
        xlab("Moral valence") +
        ylab(" Paramter mu") + 
        ggtitle("mu in each condition") +
        coord_cartesian(ylim=c(0.4,0.5)) +
        scale_y_continuous(breaks=seq(0.4,0.5,0.05),expand = c(0, 0)) +
        #scale_fill_grey (start=0.2, end=0.8) +   # using grey scale, start from darker, end to lighter. 
        apatheme + 
        theme(axis.text = element_text (size = 20,color = "black")) + 
        theme(axis.title = element_text (size = 20)) + 
        theme(plot.title = element_text(size = 20)) +
        theme(legend.text = element_text(size =20)) +
        theme(axis.title.y = element_text(margin=margin(0,20,0,0))) +  # increase the space between title and y axis
        theme(axis.title.x = element_text(margin=margin(20,0,0,0))) +   # increase the sapce betwen title and x axis
        scale_fill_manual(values=c("grey20", "grey80"),labels=c("Val-based    ","Id-based"))+
        facet_grid(.~ Identity) +
        theme(strip.text = element_text(size = 20, colour = "black")) + # set the label of facet
        theme(axis.line.x = element_line(color="black", size = 1),
                axis.line.y = element_line(color="black", size = 1))

ggsave('Test_exGass_mu.pdf',p_mu1, width=6, height=6,family = "Times")  # save the plot

p_sigma1 <- ggplot(data = Psigma.sum1, aes(y=Psigma,x=Morality,group=Task,shape = Task,fill = Task)) +
        geom_bar(position = position_dodge(),stat = "identity",colour = "black", size=.3) +         # Thinner lines
        #geom_errorbar(aes(ymin = RT-ci,ymax=RT+ci),colour="black",width=.1,position = position_dodge()) +
        #geom_line(position = pd) + geom_point(shape=21,fill="white",position = pd) +
        #geom_point(position = pd) +
        geom_errorbar(aes(ymin = Psigma-se, ymax = Psigma + se),
                      size = 1.2,
                      width = .2,
                      position=position_dodge(.9)) +
        xlab("Moral Valence") +
        ylab(" paramter sigma") + 
        ggtitle("sigma in each task") +
        #scale_y_continuous("sigma",expand = c(0, 0)) +
        coord_cartesian(ylim=c(0.04,0.09)) +
        scale_y_continuous(breaks=seq(0.04,0.09,0.01),expand = c(0, 0)) +
        #scale_fill_grey (start=0.2, end=0.8) +   # using grey scale, start from darker, end to lighter. 
        apatheme + 
        theme(axis.text = element_text (size = 20,color = "black")) + 
        theme(axis.title = element_text (size = 20)) + 
        theme(plot.title = element_text(size = 20)) +
        theme(legend.text = element_text(size =20)) +
        theme(axis.title.y = element_text(margin=margin(0,20,0,0))) +  # increase the space between title and y axis
        theme(axis.title.x = element_text(margin=margin(20,0,0,0))) +   # increase the sapce betwen title and x axis
        scale_fill_manual(values=c("grey20","grey80"),labels=c("Val-based    ","Id-based"))+
        facet_grid(.~ Identity) +
        theme(strip.text = element_text(size = 20, colour = "black")) + # set the label of facet
        theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))
ggsave('Test_exGass_sigma1.pdf',p_sigma1, width=6, height=6,family = "Times")  # save the plot

p_sigma_noTask <- ggplot(data = Psigma_anova_noTask.sum, aes(y=Psigma,x=Morality,group=Identity,shape = Identity,fill = Identity)) +
        geom_bar(position = position_dodge(),stat = "identity",colour = "black", size=.3) +         # Thinner lines
        #geom_errorbar(aes(ymin = RT-ci,ymax=RT+ci),colour="black",width=.1,position = position_dodge()) +
        #geom_line(position = pd) + geom_point(shape=21,fill="white",position = pd) +
        #geom_point(position = pd) +
        geom_errorbar(aes(ymin = Psigma-se, ymax = Psigma + se),
                      size = 1.2,
                      width = .2,
                      position=position_dodge(.9)) +
        xlab("Moral Valence") +
        ylab(" Paramter sigma") + 
        #ggtitle("sigma in each task (combined data from both tasks)") +
        #scale_y_continuous("sigma",expand = c(0, 0)) +
        coord_cartesian(ylim=c(0.04,0.09)) +
        scale_y_continuous(breaks=seq(0.04,0.09,0.01),expand = c(0, 0)) +
        #scale_fill_grey (start=0.2, end=0.8) +   # using grey scale, start from darker, end to lighter. 
        apatheme + 
        theme(axis.text = element_text (size = 20,color = "black")) + 
        theme(axis.title = element_text (size = 20)) + 
        theme(plot.title = element_text(size = 20)) +
        theme(legend.text = element_text(size =20)) +
        theme(axis.title.y = element_text(margin=margin(0,20,0,0))) +  # increase the space between title and y axis
        theme(axis.title.x = element_text(margin=margin(20,0,0,0))) +   # increase the sapce betwen title and x axis
        scale_fill_manual(values=c("grey20","grey80"),labels=c("Self     ","Other"))+
        theme(axis.line.x = element_line(color="black", size = 1),
              axis.line.y = element_line(color="black", size = 1))
ggsave('Test_exGass_sigma_noTask.pdf',p_sigma_noTask, width=6, height=6,family = "Times")  # save the plot

p_tau1 <- ggplot(data = Ptau.sum1, aes(y=Ptau,x=Morality,group=Task,shape = Task,fill = Task)) +
        geom_bar(position = position_dodge(),stat = "identity",colour = "black", size=.3) +         # Thinner lines
        #geom_errorbar(aes(ymin = RT-ci,ymax=RT+ci),colour="black",width=.1,position = position_dodge()) +
        #geom_line(position = pd) + geom_point(shape=21,fill="white",position = pd) +
        #geom_point(position = pd) +
        geom_errorbar(aes(ymin = Ptau-se, ymax = Ptau + se),
                      size = 1.2,
                      width = .2,
                      position=position_dodge(.9)) +
        xlab("Moral valence") +
        ylab(" paramter tau") + 
        ggtitle("tau in each condition") +
        scale_y_continuous("tau",expand = c(0, 0)) +
        coord_cartesian(ylim=c(0.06,0.11)) +
        scale_y_continuous(breaks=seq(0.06,0.11,0.01),expand = c(0, 0)) +
        #scale_fill_grey (start=0.2, end=0.8) +   # using grey scale, start from darker, end to lighter.
        facet_grid(.~ Identity)+
        apatheme+ 
        theme(axis.text = element_text (size = 20)) + 
        theme(axis.title = element_text (size = 20)) + 
        theme(plot.title = element_text(size = 20)) +
        theme(legend.text = element_text(size =20)) +
        theme(axis.title.y = element_text(margin=margin(0,20,0,0))) +  # increase the space between title and y axis
        theme(axis.title.x = element_text(margin=margin(20,0,0,0))) +   # increase the sapce betwen title and x axis
        scale_fill_manual(values=c("grey20","grey80"),labels=c("Val-based    ","Id-based"))+
        theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))
ggsave('Test_exGass_tau.pdf',p_tau1, width=6, height=6,family = "Times")  # save the plot


# plot ex-Gaussian ditriubtion
datT.MoralSelf <- df.exG[df.exG$Morality == 'Good' & df.exG$Identity == 'Self',]
datT.ImmoralSelf <- df.exG[df.exG$Morality == 'Bad' & df.exG$Identity == 'Self',]
datT.MoralOther <- df.exG[df.exG$Morality == 'Good' & df.exG$Identity == 'Other',]
datT.ImmoralOther <- df.exG[df.exG$Morality == 'Bad' & df.exG$Identity == 'Other',]

timefit(datT.MoralSelf$RT,iter = 2000, plot = TRUE)
timefit(datT.ImmoralSelf$RT,iter = 2000, plot = TRUE)
timefit(datT.MoralOther$RT,iter = 2000, plot = TRUE)
timefit(datT.ImmoralOther$RT,iter = 2000, plot = TRUE)

# fit the mean data
df.exG_sum <- summarySEwithin(df.exG,measurevar = 'RT',withinvars = c('Subject','Identity','Morality'),idvar = 'Subject')
datT.MoralSelf_sum <- df.exG_sum[df.exG_sum$Morality == 'moral' & df.exG_sum$Identity == 'self',]
datT.ImmoralSelf_sum <- df.exG_sum[df.exG_sum$Morality == 'immoral' & df.exG_sum$Identity == 'self',]
datT.MoralOther_sum <- df.exG_sum[df.exG_sum$Morality == 'moral' & df.exG_sum$Identity == 'other',]
datT.ImmoralOther_sum <- df.exG_sum[df.exG_sum$Morality == 'immoral' & df.exG_sum$Identity == 'other',]

timefit(datT.MoralSelf_sum$RT,iter = 2000, plot = TRUE)
timefit(datT.ImmoralSelf_sum$RT,iter = 2000, plot = TRUE)
timefit(datT.MoralOther_sum$RT,iter = 2000, plot = TRUE)
timefit(datT.ImmoralOther_sum$RT,iter = 2000, plot = TRUE)

# plot ex-Gaussian distribution by condition and individual and save their corrsponding plot
par(mfrow = c(1,1))
sub7006 <- df.exG[df.exG$Subject == 7006,]
sub7006.Id.good.self <- subset(sub7006,Task == 'identity' & Morality == 'good' & Identity == 'self')
fit.sub7006.id.good.self <- timefit(sub7006.Id.good.self$RT,iter = 2000,plot = T)

x = sub7006.Id.good.self$RT
D <- density(x)
lim <- list(x = range(D$x), y = range(D$y))
lim$y[2] <- lim$y[2] + lim$y[2]/2.5
hist(x, breaks = seq(0,1,by=0.05),xlim = c(0.2,1), ylim = c(0,10), main = "", xlab = "", 
     ylab = "", lty = 1,col = "blue",lwd = 2,prob=TRUE)
par(new = TRUE)
curve(dexgauss(x, mu = fit.sub7006.id.good.self@par[1], sigma = fit.sub7006.id.good.self@par[2], 
               tau = fit.sub7006.id.good.self@par[3]), xlim = c(0.2,1), ylim = c(0,10), 
      main = "Distribution", xlab = "Observed data", ylab = "Density", 
      lty = 1,col = "red",lwd = 2)
legend('topright',c('observed','fitted'),lty = c(1,1),lwd = 2,col = c("blue","red"),cex=2)

# plot every condition, just for checking
pdf("Individual_plot_2.pdf",height = 8, width = 8, family = "Times")
index <- 1 

par(mfrow = c(6, 8),mar = c(2,1,2,1))
for (ii in 7:12){
        for (jj in 1:nTask){
                for (kk in 1:nIdentity) {
                        for (ll in 1: nMorality){
                                dat <- df.exG[df.exG$Subject == RTsubId[ii] &
                                                            df.exG$Task == RTtask[jj] & 
                                                            df.exG$Identity == RTidentity[kk] &
                                                            df.exG$Morality == RTMorality[ll],]
                                paramEst <- params1[params1$subId == RTsubId[ii] &
                                                            params1$Task == RTtask[jj] &
                                                            params1$Identity == RTidentity[kk] &
                                                            params1$Morality == RTMorality[ll],]
                                pName <- paste(RTsubId[ii],RTtask[jj],RTidentity[kk],RTMorality[ll],sep = "_")
                                
                                D <- density(dat$RT)
                                #lim <- list(x = range(D$x), y = range(D$y))
                                #lim$y[2] <- lim$y[2] + lim$y[2]/2.5
                                plot(D, xlim = c(0.2,1), ylim = c(0,10), main = "", xlab = "", 
                                     ylab = "", lty = 1,col = "blue",lwd = 1.5,cex.axis = 0.8,
                                     xaxs="i", yaxs="i")
                                par(new = TRUE)
                                curve(dexgauss(x, mu = paramEst$Pmu, sigma = paramEst$Psigma, 
                                      tau = paramEst$Ptau), xlim = c(0.2,1), ylim = c(0,10), 
                                      main = pName, ylab = "Density", 
                                      lty = 1,col = "red",lwd = 1.5,cex.main = 0.8,cex.axis = 0.8,
                                      xaxs="i", yaxs="i")
                                legend('topright',c('observed','fitted'),lty = c(1,1),lwd = 1.5,col = c("blue","red"),cex = 0.6)
                                
                                index <- index + 1
                        }
                        
                }
                
        }

}
dev.off()

# save for matlab analysis, each condition of each participants will be a separate file
index <- 1 
for (ii in 1:27){
        for (jj in 1:nTask){
                for (kk in 1:nIdentity) {
                        for (ll in 1: nMorality){
                                dat <- df.exG[df.exG$Subject == RTsubId[ii] &
                                                            df.exG$Task == RTtask[jj] & 
                                                            df.exG$Identity == RTidentity[kk] &
                                                            df.exG$Morality == RTMorality[ll],]
                                fName <- paste(RTsubId[ii],RTtask[jj],RTidentity[kk],RTMorality[ll],sep = "_")
                                fName <- paste('./matlab/',fName,'.csv',sep = '')
                                write.table(dat$RT,fName,row.names = FALSE,col.names = F,sep = ',')
                                index <- index + 1
                        }
                        
                }
                
        }
        
}
