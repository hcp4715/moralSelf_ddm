## This file is aimed at fitting the data collected by Moral_asso_behavioral_exp7 Using Ex-Gaussian Model
#  Author     Date           History of change
# ============================================
#  Hcp     06/07/2016        Save rawData2 as the processed data
#  hcp     06/07/2016        begin with saved csv data
#  hcp     19/07/2016        only the code for final analysis
#  hcp     20/02/2017        used for the final data
#  hcp     08/04/2017        save to pdf

#### This part of the code is preparation ####
# defube the the system parameters 
Sys.setenv(LANG = "en")              # set the feedback language to English
options(scipen = 999)                # force R to output in decimal instead of scientifc notion
windowsFonts(Times=windowsFont("TT Times New Roman"))  # explicit mapping to "times"

rm(list = setdiff(ls(), lsf.str()))  # remove all data but keep functions

# define the d prime function
dprime <- function(hit,fa) {
        qnorm(hit) - qnorm(fa)
}
# define fuctions: summarySE, summarySEwithin, normDatawithin, and multiplot.
## code for calculate the summary with sE, adopted from cook book for R
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
        library(plyr)
        
        # New version of length which can handle NA's : if na.rm == T, don't count the
        length2 <- function(x, na.rm=FALSE){
                if(na.rm) sum(!is.na(x))
                else      length(x)
        }
        
        # this does the summary. For each group's data frame, return a vector with
        # N, mean, and sd
        datac <- ddply(data,groupvars, .drop=.drop,
                       .fun = function(xx,col){
                               c(N    = length2(xx[[col]],na.rm=na.rm),
                                 mean = mean(xx[[col]],na.rm=na.rm),
                                 sd   = sd  (xx[[col]],na.rm=na.rm)
                               )
                       },
                       measurevar
        )
        # Rename the "mean" column
        
        datac <- rename(datac,c("mean" = measurevar))
        
        datac$se <- datac$sd /sqrt(datac$N)   # calculate standard error of the mean
        
        # Confidence interval mltiplier for standard error
        # calculate t-statistic for confidence interval:
        # e.g., if conf.interval is .95, use .975 (above/below), and use df.L1 = N-1
        ciMult <- qt(conf.interval/2 + .5, datac$N-1)
        datac$ci <- datac$se * ciMult
        
        return (datac)
}

## code for calculate the summary with sE for within subject data, adopted from cook book for R
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
        
        # Ensure that the betweenvars and withinvars are factors
        factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                             FUN=is.factor, FUN.VALUE=logical(1))
        
        if (!all(factorvars)) {
                nonfactorvars <- names(factorvars)[!factorvars]
                message("Automatically converting the following non-factors to factors: ",
                        paste(nonfactorvars, collapse = ", "))
                data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
        }
        
        # Get the means from the un-normed data
        datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                           na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
        
        # Drop all the unused columns (these will be calculated with normed data)
        datac$sd <- NULL
        datac$se <- NULL
        datac$ci <- NULL
        
        # Norm each subject's data
        ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
        
        # This is the name of the new column
        measurevar_n <- paste(measurevar, "_norm", sep="")
        
        # Collapse the normed data - now we can treat between and within vars the same
        ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                            na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
        
        # Apply correction from Morey (2008) to the standard error and confidence interval
        #  Get the product of the number of conditions of within-S variables
        nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                        FUN.VALUE=numeric(1)))
        correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
        
        # Apply the correction factor
        ndatac$sd <- ndatac$sd * correctionFactor
        ndatac$se <- ndatac$se * correctionFactor
        ndatac$ci <- ndatac$ci * correctionFactor
        
        # Combine the un-normed means with the normed results
        merge(datac, ndatac)
}

### code for normalizing the SE
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
        library(plyr)
        
        # Measure var on left, idvar + between vars on right of formula.
        data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                               .fun = function(xx, col, na.rm) {
                                       c(subjMean = mean(xx[,col], na.rm=na.rm))
                               },
                               measurevar,
                               na.rm
        )
        
        # Put the subject means with original data
        data <- merge(data, data.subjMean)
        
        # Get the normalized data in a new column
        measureNormedVar <- paste(measurevar, "_norm", sep="")
        data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
                mean(data[,measurevar], na.rm=na.rm)
        
        # Remove this subject mean column
        data$subjMean <- NULL
        
        return(data)
}

#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        library(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}

#### load necessary libraries ####
# define a function to install the necessary function is not installed before
pkgTest <- function(x)
{
        if (!require(x,character.only = TRUE))
        {
                install.packages(x,dep = TRUE)
                if(!require(x,character.only = TRUE)) stop("Package not found")
        }
}

pkgNeeded <- (c("plyr","ez","ggplot2", "reshape2","retimes","bootES"))

lapply(pkgNeeded,pkgTest) # apply pkgTest function to pkgNeeded to load and install library
rm(pkgNeeded)     # remove variable pkgNeeded

#library(plyr)     # using plyr package
#library(ez)       # using ez package for ANOVA
#library(ggplot2)  # using ggplot2 for plot
#library(reshape2)
#library(retimes)

# using APA style plot
# Save some time and stor APA format-related code in an object so you can easily
apatheme=theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              text=element_text(family='Times'),
              legend.title=element_blank(),
              legend.position='top',
              axis.line.x = element_line(color='black'),
              axis.line.y = element_line(color='black'))


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

##########begin analysis#############################################

##  extract the trials with correct response
datT.validRT <- df.T.hddm[df.T.hddm$response == 1 & df.T.hddm$rt > 0.2,]  # RT > 200ms for accurate trials

# rename the colnames
colnames(datT.validRT) <- c('SubjectID','Task','Morality','Identity','RT','ACC')

## analyzing the distribution parameters for correct trials

RTsubId <- unique(datT.validRT$SubjectID)  # all unique subject's ID
RTtask <- unique(datT.validRT$Task)        # all unique tasks
RTidentity <- unique(datT.validRT$Identity)# all unique Identities, key independent varialbe
RTMorality <- unique(datT.validRT$Morality)# all unique Moral valence level, key independent variable

nSub <- length(RTsubId)          # number of subjects, save for later for loop
nTask <- length(RTtask)          # number of tasks, save for later for loop
nIdentity <- length(RTidentity)  # number of identies, save for later for loop
nMorality <- length(RTMorality)  # number of moral levels, save for later for loop

rowNum <- nSub*nTask*nIdentity*nMorality # predefine the number of rows for saving data

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
results.params_valid <- results.params[results.params$trialNum > 36,]
valid_condition <- count(results.params_valid,'subId') # check if there some invalid data 
# results.params_final <- results.params_valid[-which(results.params_valid$subId %in% valid_condition[valid_condition$freq < 8,]$subId),]

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

params1 <- results.params_valid

# save to csv
write.csv(params1,'ExGaussian_params_2.csv',row.names = FALSE)
#params1 <- read.csv('ExGaussian_params.csv',header = TRUE, sep = ',', stringsAsFactors=FALSE)

# transfer to wide-format
# mu
params_wide_mu <- dcast(params1,subId ~ Task + Identity + Morality, value.var = "Pmu")
# sigma
params_wide_sigma <- dcast(params1,subId ~ Task + Identity + Morality, value.var = "Psigma")
# tau
params_wide_tau <- dcast(params1,subId ~ Task + Identity + Morality, value.var = "Ptau")
# save the data
write.csv(params_wide_mu,'ExGaussian_params_wide_mu.csv',row.names = FALSE)
write.csv(params_wide_sigma,'ExGaussian_params_wide_sigma.csv',row.names = FALSE)
write.csv(params_wide_tau,'ExGaussian_params_wide_tau.csv',row.names = FALSE)

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

# analyzing the interaction between morality and identity (combined different tasks)

params_sigma_noTask <- summarySEwithin(params1,measurevar = 'Psigma',withinvars = c('subId','Identity','Morality'),idvar = 'subId')
params_wide_sigma_noTask <- dcast(params_sigma_noTask,subId ~ Identity + Morality, value.var = "Psigma")
write.csv(params_wide_sigma_noTask,'ExGaussian_params_wide_sigma_noTask.csv',row.names = FALSE)

Psigma_anova_noTask <- ezANOVA(params_sigma_noTask,dv = Psigma, wid = subId, within=.(Morality,Identity), type=3)
Psigma_anova_noTask.sum <- summarySEwithin(params1,measurevar = 'Psigma',withinvars = c('Identity','Morality'),idvar = 'subId')
Psigma_anova_noTask.sum$Identity <- factor(Psigma_anova_noTask.sum$Identity,levels = c('self','other'))
Psigma_anova_noTask.sum$Morality <- factor(Psigma_anova_noTask.sum$Morality,levels = c('good','bad'))

# t tests
Psigma.t.mrl_imm_self <- t.test(params_wide_sigma_noTask$self_good,params_wide_sigma_noTask$self_bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.mrl_imm_self[3]),method = 'bonferroni',4)
params_wide_sigma_noTask$mrl_imm_self <- params_wide_sigma_noTask$self_good - params_wide_sigma_noTask$self_bad
Psigma.t.mrl_imm_self.ci <- bootES(params_wide_sigma_noTask$mrl_imm_self,R = 20000, effect.type = "cohens.d")

Psigma.t.mrl_imm_other <- t.test(params_wide_sigma_noTask$other_good,params_wide_sigma_noTask$other_bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.mrl_imm_other[3]),method = 'bonferroni',4)
params_wide_sigma_noTask$mrl_imm_other <- params_wide_sigma_noTask$other_good - params_wide_sigma_noTask$other_bad
Psigma.t.mrl_imm_other.ci <- bootES(params_wide_sigma_noTask$mrl_imm_other,R = 20000, effect.type = "cohens.d")

Psigma.t.slf_oth_mrl <- t.test(params_wide_sigma_noTask$self_good,params_wide_sigma_noTask$other_good,paired = TRUE)
p.adjust(as.numeric(Psigma.t.slf_oth_mrl[3]),method = 'bonferroni',4)
params_wide_sigma_noTask$mrl_self_other <- params_wide_sigma_noTask$self_good - params_wide_sigma_noTask$other_good
Psigma.t.mrl_self_other.ci <- bootES(params_wide_sigma_noTask$mrl_self_other,R = 20000, effect.type = "cohens.d")

Psigma.t.slf_oth_imm <- t.test(params_wide_sigma_noTask$self_bad,params_wide_sigma_noTask$other_bad,paired = TRUE)
p.adjust(as.numeric(Psigma.t.slf_oth_imm[3]),method = 'bonferroni',4)
params_wide_sigma_noTask$slf_oth_imm <- params_wide_sigma_noTask$self_bad - params_wide_sigma_noTask$other_bad
Psigma.t.slf_oth_imm.ci <- bootES(params_wide_sigma_noTask$slf_oth_imm,R = 20000, effect.type = "cohens.d")

Ptau.sum1 <- summarySEwithin(params1,measurevar = 'Ptau',withinvars = c('Task','Identity','Morality'),idvar = 'subId')
Ptau.sum1$Task <- factor(Ptau.sum1$Task, levels=c("Morality","Identity"))
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
        scale_fill_manual(values=c("grey20", "grey80"),labels=c("Morality     ","Identity"))+
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
        scale_fill_manual(values=c("grey20","grey80"),labels=c("Morality     ","Identity"))+
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

# plot every condition, just for checking
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
