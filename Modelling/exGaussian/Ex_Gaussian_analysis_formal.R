## This file is aimed at fitting the data collected by Moral_asso_behavioral_exp7 Using Ex-Gaussian Model
#  Author     Date           History of change
# ============================================
#  Hcp     06/07/2016        Save rawData2 as the processed data
#  hcp     06/07/2016        begin with saved csv data
#  hcp     19/07/2016        only the code for final analysis
#  hcp     20/02/2017        used for the final data
#  hcp     08/04/2017        save to pdf

# run summarySE, summarySEwithin, normDatawithin, and multiplot before running this script.
Sys.setenv(LANG = "en") # set the feedback language to English
options(scipen = 999)   # force R to output in decimal instead of scientifc notion
rm(list = setdiff(ls(), lsf.str()))  # remove all data but keep functions


############## load and clear data ##########################
# read data from 6 participants
df.T.hddm <- read.csv("data_Categ_hddm.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE) 
df.T.hddm$moral[df.T.hddm$moral=='moral'] <- 'good' 
df.T.hddm$moral[df.T.hddm$moral=='immoral'] <- 'bad'
df.T_sum <- summarySEwithin(df.T.hddm,measurevar = 'rt',withinvars = c('subj_idx','crit','moral','id'),idvar = 'subj_idx')


# render the numeric column
cols.num <- c("rt","response")
df.T.hddm[cols.num] <- sapply(df.T.hddm[cols.num], as.numeric)

# check the class for each column
sapply(df.T.hddm,class)

rm(cols.num) # remove temporary variables

#########################################################################################

##################### begin analysis
library(plyr)     # using plyr package
library(ez)       # using ez package for ANOVA
library(ggplot2)  # using ggplot2 for plot
library(reshape2)
library(retimes)


##  extract the trials with correct response
datT.validRT <- df.T.hddm[df.T.hddm$response == 1 & df.T.hddm$rt > 0.2,]  # RT > 200ms for accurate trials

colnames(datT.validRT) <- c('SubjectID','Task','Morality','Identity','RT','ACC')

## analyzing the distribution parameters for correct trials
RTsubId <- unique(datT.validRT$SubjectID)
RTtask <- unique(datT.validRT$Task)
RTidentity <- unique(datT.validRT$Identity)
RTMorality <- unique(datT.validRT$Morality)

nSub <- length(RTsubId)
nTask <- length(RTtask)
nIdentity <- length(RTidentity)
nMorality <- length(RTMorality)

rowNum <- nSub*nTask*nIdentity*nMorality

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

index <- 1     
for (ii in 1:nSub){
        for (jj in 1:nTask){
                for (kk in 1:nIdentity) {
                        for (ll in 1: nMorality){
                                dat <- datT.validRT[datT.validRT$SubjectID == RTsubId[ii] &
                                                             datT.validRT$Task == RTtask[jj] & 
                                                             datT.validRT$Identity == RTidentity[kk] &
                                                             datT.validRT$Morality == RTMorality[ll],]
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

#rm(RTidentity,RTMorality,RTsubId,RTtask,index,exGassparams) # remove the temporary variables

# remove invalid data in results of params
results.params_valid <-results.params[results.params$trialNum > 36,]
results.params_valid <- results.params_valid[-which(results.params_valid$subId %in% results.params[results.params$trialNum < 36,]$subId),]

# add one factor to params data
# results.params$condition[results.params$Identity == "self" & results.params$Morality == "moral"] <- "moralself"
# results.params$condition[results.params$Identity == "self" & results.params$Morality == "immoral"] <- "immoralself"
# results.params$condition[results.params$Identity == "other" & results.params$Morality == "moral"] <- "moralother"
# results.params$condition[results.params$Identity == "other" & results.params$Morality == "immoral"] <- "immoralother"
# results.params$condition <- factor(results.params$condition, levels=c("moralself", "immoralself","moralother","immoralother")) # specify the order of x-axis
# results.params$Task <- factor(results.params$Task, levels=c("identity", "morality")) # specify the order of x-axis
#results.params$Task <- factor(results.params$Task, levels=c("self", "moral","importance")) # specify the order of x-axis

# params1  <- results.params_valid[results.params_valid$Task != "importance",] # separate dataset for firs two tasks
# params2 <- results.params[results.params$Task == "importance",] # dataset for impor

params1 <- results.params

# save to csv
write.csv(params1,'ExGaussian_params_2.csv',row.names = FALSE)
#params1 <- read.csv('ExGaussian_params.csv',header = TRUE, sep = ',', stringsAsFactors=FALSE)
# analysis of params
Pmu_anova <- ezANOVA(params1,dv = Pmu, wid = subId, within=.(Task,Morality,Identity), type=3)
Psigma_anova <- ezANOVA(params1,dv = Psigma, wid = subId, within=.(Task,Morality,Identity), type=3)
Ptau_anova <- ezANOVA(params1,dv = Ptau, wid = subId, within=.(Task,Morality,Identity), type=3)

# calculate the summary data for each parameter
Pmu.sum1 <- summarySEwithin(params1,measurevar = 'Pmu',withinvars = c('Task','Identity','Morality'),idvar = 'subId')
Pmu.sum1$Task <- factor(Pmu.sum1$Task, levels=c("morality","identity"))
Pmu.sum1$Identity <- factor(Pmu.sum1$Identity,levels = c('self','other'))
Pmu.sum1$Morality <- factor(Pmu.sum1$Morality,levels = c('good','bad'))

# analyze the main effect of identity
Pmu.sum2 <- summarySEwithin(params1,measurevar = 'Pmu',withinvars = c('Identity'),idvar = 'subId')

Psigma.sum1 <- summarySEwithin(params1,measurevar = 'Psigma',withinvars = c('Task','Identity','Morality'),idvar = 'subId')
Psigma.sum1$Task <- factor(Psigma.sum1$Task, levels=c("morality","identity"))
Psigma.sum1$Identity <- factor(Psigma.sum1$Identity,levels = c('self','other'))
Psigma.sum1$Morality <- factor(Psigma.sum1$Morality,levels = c('good','bad'))

# analyzing the interaction between morality and identity
Psigma.sum2 <- summarySEwithin(params1,measurevar = 'Psigma',withinvars = c('subId','Identity','Morality'),idvar = 'subId')
Psigma_anova2 <- ezANOVA(Psigma.sum2,dv = Psigma, wid = subId, within=.(Morality,Identity), type=3)
Psigma.sum2_w <- dcast(Psigma.sum2, subId ~ Identity + Morality ,value.var = "Psigma")
Psigma.sum3 <- summarySEwithin(params1,measurevar = 'Psigma',withinvars = c('Identity','Morality'),idvar = 'subId')

# t tests
Psigma.t.mrl_imm_self <- t.test(Psigma.sum2_w$self_good,Psigma.sum2_w$self_bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.mrl_imm_self[3]),method = 'bonferroni',4)

Psigma.t.mrl_imm_other <- t.test(Psigma.sum2_w$other_good,Psigma.sum2_w$other_bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.mrl_imm_other[3]),method = 'bonferroni',4)

Psigma.t.slf_oth_mrl <- t.test(Psigma.sum2_w$self_good,Psigma.sum2_w$other_good,paired = TRUE)
p.adjust(as.numeric(Psigma.t.slf_oth_mrl[3]),method = 'bonferroni',4)

Psigma.t.slf_oth_imm <- t.test(Psigma.sum2_w$self_bad,Psigma.sum2_w$other_bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.slf_oth_imm[3]),method = 'bonferroni',4)

Ptau.sum1 <- summarySEwithin(params1,measurevar = 'Ptau',withinvars = c('Task','Identity','Morality'),idvar = 'subId')
Ptau.sum1$Task <- factor(Ptau.sum1$Task, levels=c("morality","identity"))
Ptau.sum1$Identity <- factor(Ptau.sum1$Identity,levels = c('self','other'))
Ptau.sum1$Morality <- factor(Ptau.sum1$Morality,levels = c('good','bad'))

library(ggplot2)
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
        facet_grid(.~ Identity) +
        apatheme + 
        theme(axis.text = element_text (size = 20)) + 
        theme(axis.title = element_text (size = 20)) + 
        theme(plot.title = element_text(size = 20)) +
        theme(legend.text = element_text(size =20)) +
        theme(axis.title.y = element_text(margin=margin(0,20,0,0))) +  # increase the space between title and y axis
        theme(axis.title.x = element_text(margin=margin(20,0,0,0))) +   # increase the sapce betwen title and x axis
        scale_fill_manual(values=c("grey20", "grey80"),labels=c("Morality     ","Identity"))+
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
        scale_y_continuous("sigma",expand = c(0, 0)) +
        coord_cartesian(ylim=c(0.04,0.09)) +
        scale_y_continuous(breaks=seq(0.04,0.09,0.01),expand = c(0, 0)) +
        #scale_fill_grey (start=0.2, end=0.8) +   # using grey scale, start from darker, end to lighter. 
        facet_grid(.~ Identity)+
        apatheme + 
        theme(axis.text = element_text (size = 20)) + 
        theme(axis.title = element_text (size = 20)) + 
        theme(plot.title = element_text(size = 20)) +
        theme(legend.text = element_text(size =20)) +
        theme(axis.title.y = element_text(margin=margin(0,20,0,0))) +  # increase the space between title and y axis
        theme(axis.title.x = element_text(margin=margin(20,0,0,0))) +   # increase the sapce betwen title and x axis
        scale_fill_manual(values=c("grey20","grey80"),labels=c("Morality     ","Identity"))+
        theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))
ggsave('Test_exGass_sigma.pdf',p_sigma1, width=6, height=6,family = "Times")  # save the plot


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
        scale_fill_manual(values=c("grey20","grey80"),labels=c("Morality     ","Identity"))+
        theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))
ggsave('Test_exGass_tau.pdf',p_tau1, width=6, height=6,family = "Times")  # save the plot


# plot ex-Gaussian ditriubtion
datT.MoralSelf <- datT.validRT[datT.validRT$Morality == 'moral' & datT.validRT$Identity == 'self',]
datT.ImmoralSelf <- datT.validRT[datT.validRT$Morality == 'immoral' & datT.validRT$Identity == 'self',]
datT.MoralOther <- datT.validRT[datT.validRT$Morality == 'moral' & datT.validRT$Identity == 'other',]
datT.ImmoralOther <- datT.validRT[datT.validRT$Morality == 'immoral' & datT.validRT$Identity == 'other',]

timefit(datT.MoralSelf$RT,iter = 2000, plot = TRUE)
timefit(datT.ImmoralSelf$RT,iter = 2000, plot = TRUE)
timefit(datT.MoralOther$RT,iter = 2000, plot = TRUE)
timefit(datT.ImmoralOther$RT,iter = 2000, plot = TRUE)

# fit the mean data
datT.validRT_sum <- summarySEwithin(datT.validRT,measurevar = 'RT',withinvars = c('SubjectID','Identity','Morality'),idvar = 'SubjectID')
datT.MoralSelf_sum <- datT.validRT_sum[datT.validRT_sum$Morality == 'moral' & datT.validRT_sum$Identity == 'self',]
datT.ImmoralSelf_sum <- datT.validRT_sum[datT.validRT_sum$Morality == 'immoral' & datT.validRT_sum$Identity == 'self',]
datT.MoralOther_sum <- datT.validRT_sum[datT.validRT_sum$Morality == 'moral' & datT.validRT_sum$Identity == 'other',]
datT.ImmoralOther_sum <- datT.validRT_sum[datT.validRT_sum$Morality == 'immoral' & datT.validRT_sum$Identity == 'other',]

timefit(datT.MoralSelf_sum$RT,iter = 2000, plot = TRUE)
timefit(datT.ImmoralSelf_sum$RT,iter = 2000, plot = TRUE)
timefit(datT.MoralOther_sum$RT,iter = 2000, plot = TRUE)
timefit(datT.ImmoralOther_sum$RT,iter = 2000, plot = TRUE)

# plot ex-Gaussian distribution by condition and individual and save their corrsponding plot
par(mfrow = c(1,1))
sub7006 <- datT.validRT[datT.validRT$SubjectID == 7006,]
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

# plot
pdf("Individual_plot_2.pdf",height = 8, width = 8, family = "Times")
index <- 1 

par(mfrow = c(6, 8),mar = c(2,1,2,1))
for (ii in 7:12){
        for (jj in 1:nTask){
                for (kk in 1:nIdentity) {
                        for (ll in 1: nMorality){
                                dat <- datT.validRT[datT.validRT$SubjectID == RTsubId[ii] &
                                                            datT.validRT$Task == RTtask[jj] & 
                                                            datT.validRT$Identity == RTidentity[kk] &
                                                            datT.validRT$Morality == RTMorality[ll],]
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
                                dat <- datT.validRT[datT.validRT$SubjectID == RTsubId[ii] &
                                                            datT.validRT$Task == RTtask[jj] & 
                                                            datT.validRT$Identity == RTidentity[kk] &
                                                            datT.validRT$Morality == RTMorality[ll],]
                                fName <- paste(RTsubId[ii],RTtask[jj],RTidentity[kk],RTMorality[ll],sep = "_")
                                fName <- paste('./matlab/',fName,'.csv',sep = '')
                                write.table(dat$RT,fName,row.names = FALSE,col.names = F,sep = ',')
                                index <- index + 1
                        }
                        
                }
                
        }
        
}
