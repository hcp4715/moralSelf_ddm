### This script is used to preprocessing the data of the moral self experiment (DDM), experiment 2. 
#
### Hu, C-P., Lan, Y., Macrae, N., & Sui. J. (2019) 
### Script author: Chuan-Peng Hu
### Email = hcp4715@gmail.com       twitter= @hcp4715
#
### Input: 
###      MS_rep_matchTask_raw.csv -- raw data of matching task;
###      MS_rep_categTask_raw.csv -- raw data of categorization task;
#
### Output:
###     df.M.hddm_m_stim.csv.csv  -- clean data of matching trials from matchign task, for HDDM analysis;
###     df.M.hddm_nm_stim.csv.csv -- clean data of nonmatching trials from matchign task, for HDDM analysis;
###     df.C.hddm_val_stim.csv    -- clean data of valence-based categorization trials, for HDDM analysis;
###     df.C.hddm_id_stim.csv     -- clean data of identity-based categorization trials, for HDDM analysis;
# 
###     MS_categ_behav_wide.csv        -- cleaned summary results (wide-format) of categorization task for 
###                                             statistical analysis in JASP
###     MS_categ_behav_noTask_wide.csv -- cleaned summary results (wide-format) of categorization task (collapsed
###                                             data from different task) for statistical analysis in JASP
###     MS_categ__rt_acc_long.csv      -- cleaned summary results of RT and accuracy (long-format) of 
###                                             categorization task, for JASP analysis
###     MS_categ__rt_acc_noTask_long.csv -- cleaned summary reuslt of RT and accuracy (long-format) of categorization 
###                                             task (collapsed data from different task) for JASP analysis
###     MS_cross_taskeffect_wide.csv    -- cleaned summary results for cross task analysis

# ---------- Table of Contents ----------------------------------------------------------
# ---------- 1. Initializing and prepare ------------------------------------------------
# ---------- 2. Loading data and clean the data -----------------------------------------
# ---------- 3. Matching task: prepare the Accuracy, d-prime and RT ---------------------
# ---------- 4. Categorization task: prepare the Accuracy and RT ------------------------
# ---------- 5. Plots -------------------------------------------------------------------


# ---------------------------------------------------------------------------------------
# ---------- 1. Initializing and prepare               ----------------------------------
# ---------------------------------------------------------------------------------------

curDir  <- dirname(rstudioapi::getSourceEditorContext()$path)   # get the directory for preprocessing
setwd(curDir)
source('Initial_ms_rep.r')  # initializing (clear global environment; load packages and functions)
curDir  <- dirname(rstudioapi::getSourceEditorContext()$path)   # get the directory for preprocessing
rootDir <- gsub('.{9}$', '', curDir)                            # get the parental folder
traDir  <- paste(rootDir,'2_trad_analysis',sep = '')        # folder for traditional analsysi
ddmDir  <- paste(rootDir,'3_hddm',sep = '')                        # folder for DDM analysis

# ---------------------------------------------------------------------------------------
# ---------- 2. Loading data and clean the data        ----------------------------------
# ---------------------------------------------------------------------------------------

df.M <- read.csv('MS_rep_matchingTask_raw.csv', header=TRUE, sep=",",stringsAsFactors=F)
df.C <- read.csv('MS_rep_categTask_raw.csv', header=TRUE, sep=",",stringsAsFactors=F)

###############################
###  Excluding particpants ####
###############################

### Rule 1: wrong trials numbers because of procedure errors
excldSub1_M <- df.M %>%
   dplyr::mutate(ACC = ifelse(ACC == 1, 1, 0))  %>%  # no response as wrong
   dplyr::group_by(Subject, Match, Identity,Morality) %>%
   dplyr::summarise(N = length(ACC)) %>%  # count the trial # for each condition of each subject
   dplyr::ungroup() %>%
   dplyr::filter(N != 75) %>%             # filter the rows that trial Number is not 75
   dplyr::distinct(Subject) %>%           # find the unique subject ID
   dplyr::pull(Subject)                   # pull the subj ID as vector


excldSub1_C <- df.C %>%
   dplyr::mutate(ACC = ifelse(ACC == 1, 1, 0))  %>%  # no response as wrong
   dplyr::group_by(Subject, Task, Identity,Morality) %>%
   dplyr::summarise(N = length(ACC)) %>%  # count the trial # for each condition of each subject
   dplyr::ungroup() %>%
   dplyr::filter(N != 90) %>%             # filter the rows that trial Number is not 90
   dplyr::distinct(Subject) %>%           # find the unique subject ID
   dplyr::pull(Subject)                   # pull the subj ID as vector


### Rule 2:  overall accuracy < 0.5
excldSub2_M <- df.M %>%
   dplyr::mutate(ACC = ifelse(ACC == 1, 1, 0))  %>%  # no response as wrong
   dplyr::group_by(Subject) %>%
   dplyr::summarise(N = length(ACC),
                    countN = sum(ACC),
                    ACC = sum(ACC)/length(ACC)) %>%  # count the trial # for each condition of each subject
   dplyr::ungroup() %>%
   dplyr::filter(ACC < .5) %>%             # filter the subjects with over all ACC < 0.5
   dplyr::distinct(Subject) %>%             # find the unique subject ID
   dplyr::pull(Subject)                     # pull the subj ID as vector

excldSub2_C <- df.C %>%
   dplyr::mutate(ACC = ifelse(ACC == 1, 1, 0))  %>%  # no response as wrong
   dplyr::group_by(Subject) %>%
   dplyr::summarise(N = length(ACC),
                    countN = sum(ACC),
                    ACC = sum(ACC)/length(ACC)) %>%  # count the trial # for each condition of each subject
   dplyr::ungroup() %>%
   dplyr::filter(ACC < .5) %>%             # filter the subjects with over all ACC < 0.5
   dplyr::distinct(Subject) %>%             # find the unique subject ID
   dplyr::pull(Subject)                     # pull the subj ID as vector

### Rule 3:  one condition with zero ACC
excldSub3_M <- df.M %>%
   dplyr::mutate(ACC = ifelse(ACC == 1, 1, 0))  %>%  # no response as wrong
   dplyr::group_by(Subject, Match, Identity,Morality) %>%
   dplyr::summarise(N = length(ACC),
                    countN = sum(ACC),
                    ACC = sum(ACC)/length(ACC)) %>%  # count the trial # for each condition of each subject
   dplyr::ungroup() %>%
   dplyr::filter(ACC == 0) %>%             # filter the subjects with over all ACC < 0.5
   dplyr::distinct(Subject) %>%             # find the unique subject ID
   dplyr::pull(Subject)                     # pull the subj ID as vector

excldSub3_C <- df.C %>%
   dplyr::mutate(ACC = ifelse(ACC == 1, 1, 0))  %>%  # no response as wrong
   dplyr::group_by(Subject, Task, Identity,Morality) %>%
   dplyr::summarise(N = length(ACC),
                    countN = sum(ACC),
                    ACC = sum(ACC)/length(ACC)) %>%  # count the trial # for each condition of each subject
   dplyr::ungroup() %>%
   dplyr::filter(ACC == 0) %>%             # filter the subjects with over all ACC < 0.5
   dplyr::distinct(Subject) %>%             # find the unique subject ID
   dplyr::pull(Subject)                     # pull the subj ID as vector
   

# all participants excluded
excldSub_M   <- c(excldSub1_M, excldSub2_M,excldSub3_M) # 7302, 7303
excldSub_C   <- c(excldSub1_C, excldSub2_C,excldSub3_C) # 7302, 7303, 7338

# select valid data for further analysis
df.M1.V <- df.M %>%
   dplyr::mutate(ACC = ifelse(ACC == 1, 1, 0))  %>%  # no response as wrong
   dplyr::filter(Subject,!Subject %in% excldSub_M)   # exclude the invalid subjects

df.C1.V <- df.C %>%
   dplyr::mutate(ACC = ifelse(ACC == 1, 1, 0))  %>%  # no response as wrong
   dplyr::filter(Subject,!Subject %in% excldSub_C)   # exclude the invalid subjects

# check the number of participants are correct
length(unique(df.M1.V$Subject)) + length(excldSub_M) == length(unique(df.M$Subject))
length(unique(df.C1.V$Subject)) + length(excldSub_C) == length(unique(df.C$Subject))


# excluded correct trials with < 200ms RT
ratio.excld.trials.M <- nrow(df.M1.V[df.M1.V$RT*1000 <= 200 & df.M1.V$ACC == 1,])/nrow(df.M1.V)  # ratio of invalid trials
ratio.excld.trials.C <- nrow(df.C1.V[df.C1.V$RT*1000 <= 200 & df.C1.V$ACC == 1,])/nrow(df.C1.V)

df.M1.V <- df.M1.V %>% dplyr::filter(!(RT <= 0.2 & ACC==1)) # filter invalid trials
df.C1.V <- df.C1.V %>% dplyr::filter(!(RT <= 0.2 & ACC==1)) # filter invalid trials

## Basic information of the data ####
df.M1.T.basic <- df.M %>%
   dplyr::select(Subject, Age, Sex) %>%
   dplyr::distinct(Subject, Age, Sex) %>%
   dplyr::summarise(subj_N = length(Subject),
                    female_N = sum(Sex == 'female'),
                    male_N = sum(Sex == 'male'),
                    Age_mean = round(mean(Age),2),
                    Age_sd   = round(sd(Age),2))

# valide data for matching task
df.M1.V.basic <- df.M1.V %>%
   dplyr::select(Subject, Age, Sex) %>%
   dplyr::distinct(Subject, Age, Sex) %>%
   dplyr::summarise(subj_N = length(Subject),
                    female_N = sum(Sex == 'female'),
                    male_N = sum(Sex == 'male'),
                    Age_mean = round(mean(Age),2),
                    Age_sd   = round(sd(Age),2))

# valid data for categorization task
df.C1.V.basic <- df.C1.V %>%
   dplyr::select(Subject, Age, Sex) %>%
   dplyr::distinct(Subject, Age, Sex) %>%
   dplyr::summarise(subj_N = length(Subject),
                    female_N = sum(Sex == 'female'),
                    male_N = sum(Sex == 'male'),
                    Age_mean = round(mean(Age),2),
                    Age_sd   = round(sd(Age),2))

###########################################################
########   prep the data file for hddm analysis     #######
###########################################################
# exclusion trials, criterion: trials without response or wrong key
# Note: we didn't remove the correct response less than 200 ms
##################################################
### preparing the data for Stimuli-code hddm  ####
##################################################

df.M.hddm_m <- df.M %>%                                 # start from the valid data, RT in seconds.
   dplyr::filter(Subject,!Subject %in% excldSub_M) %>%
   dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
   dplyr::filter(Match == 'match') %>%
   dplyr::select(Subject,Match, Morality,Identity,RT,ACC) %>%  # select column
   dplyr::rename(subj_idx = Subject, match = Match, val = Morality, id = Identity, rt = RT, response = ACC) # rename column

df.M.hddm_nm <- df.M %>%
   dplyr::filter(Subject,!Subject %in% excldSub_M) %>%
   dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
   dplyr::filter(Match == 'mismatch') %>%
   dplyr::select(Subject,Morality,Identity,RT,ACC) %>%  # select column
   dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, response = ACC) # rename column

df.M.hddm_stim <- df.M %>%
   dplyr::filter(Subject,!Subject %in% excldSub_M) %>%
   dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
   dplyr::mutate(response = ifelse((Match == "match" & ACC ==1) | (Match == "mismatch" & ACC ==0), 1, 0)) %>%
   dplyr::mutate(stim = ifelse(Match == "match", 1, 0)) %>%
   dplyr::select(Subject,Match,Morality,Identity,stim,response,RT) %>%  # select column
   dplyr::rename(subj_idx = Subject, match = Match, val = Morality, id = Identity, rt = RT) # rename column


##################################################
### preparing the data for Stimuli-code hddm  ####
##################################################
# for valence based, stim: moral = 1, immoral = 0; response: moralKey = 1; immoralKey = 0
df.C.hddm_val_stim <- df.C %>%
   dplyr::filter(Subject,!Subject %in% excldSub_C) %>%
   dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
   dplyr::filter(Task == 'Val') %>%                     # select the valence-based task
   dplyr::mutate(response = ifelse((Morality == "Good" & ACC ==1) | (Morality == "Bad" & ACC ==0), 1, 0)) %>%
   dplyr::select(Subject,Morality,Identity,Shape, response, RT) %>%  # select column
   dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, stim = Shape) %>% # rename column
   dplyr::mutate(stim = ifelse(stim == "moralSelf" | stim == "moralOther" ,1, 0))

# for Id based, stim: self = 1, other = 0; response: selfKey = 1; otherKey = 0
df.C.hddm_id_stim <- df.C%>%
   dplyr::filter(Subject,!Subject %in% excldSub_C) %>%
   dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
   dplyr::filter(Task == 'Id') %>%                      # select the valence-based task
   dplyr::mutate(response = ifelse((Identity =="Self" & ACC==1) | (Identity =="Other" & ACC==0), 1, 0)) %>%
   dplyr::select(Subject,Morality,Identity,Shape, response, RT) %>%  # select column
   dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, stim = Shape) %>% # rename column
   dplyr::mutate(stim = ifelse(stim == "moralSelf" | stim == "immoralSelf",1, 0))

setwd(ddmDir)
write.csv(df.M.hddm_m,'MS_match_hddm.csv',row.names = F)
write.csv(df.M.hddm_nm,'MS_mismatch_hddm.csv',row.names = F)
write.csv(df.M.hddm_stim,'MS_match_hddm_stim.csv',row.names = F)
write.csv(df.C.hddm_val_stim,'MS_categ_val_hddm_stim.csv',row.names = F)
write.csv(df.C.hddm_id_stim,'MS_categ_id_hddm_stim.csv',row.names = F)
setwd(curDir)
rm('df.M.hddm_m','df.M.hddm_nm','df.M.hddm_stim','df.C.hddm_val_stim','df.C.hddm_id_stim')


# ---------------------------------------------------------------------------------------
# -------- 3. Matching task: prepare the Accuracy, d-prime and RT -----------------------
# ---------------------------------------------------------------------------------------

####################################
###############   ACC    ###########
####################################
df.M1.V.acc  <- df.M1.V %>%
      dplyr::group_by(Subject, Match, Morality, Identity)  %>% 
      dplyr::summarise(N = length(ACC),
                       countN = sum(ACC),
                       ACC = sum(ACC)/length(ACC))

# transform to wide-format
df.M1.V.acc_w <- reshape2::dcast(df.M1.V.acc, Subject ~ Match + Morality + Identity,value.var = "ACC")

# rename the column number
colnames(df.M1.V.acc_w)[2:9] <- paste("ACC", colnames(df.M1.V.acc_w[,2:9]), sep = "_")

####################################
############   d prime   ###########
####################################
# calculate the number of hit,CR,miss or FA 
df.M1.V$sdt <- NA
for (i in 1:nrow(df.M1.V)){
      if (df.M1.V$ACC[i] == 1 & df.M1.V$Match[i] == "match"){
            df.M1.V$sdt[i] <- "hit"
      } else if (df.M1.V$ACC[i] == 1 & df.M1.V$Match[i] == "mismatch" ){
            df.M1.V$sdt[i] <- "CR"
      } else if (df.M1.V$ACC[i] == 0 & df.M1.V$Match[i] == "match"){
            df.M1.V$sdt[i] <- "miss"
      } else if (df.M1.V$ACC[i] == 0 & df.M1.V$Match[i] == "mismatch" ){
            df.M1.V$sdt[i] <- "FA"
      }
}

# calculate the number of each for each condition
df.M1.V.SDT <- df.M1.V %>%
   dplyr::group_by(Subject,Age, Sex, Morality, Identity,sdt) %>%
   dplyr::summarise(N = length(sdt)) %>%
   dplyr::filter(!is.na(sdt))           # no NAs

# long format to wide
df.M1.V.SDT_w <- reshape2::dcast(df.M1.V.SDT, Subject + Age + Sex + Morality + Identity ~ sdt,value.var = "N")
df.M1.V.SDT_w <- df.M1.V.SDT_w %>%
   dplyr::mutate(miss = ifelse(is.na(miss),0,miss), # if not miss trial, to 0
                 FA   = ifelse(is.na(FA),0,FA),     # if not FA trial, to 0
                 hitR = hit/(hit + miss),           # calculate the hit rate
                 faR  = FA/(FA+CR))                 # calculate the FA rate


# standardized way to deal with the extreme values
for (i in 1:nrow(df.M1.V.SDT_w)){
      if (df.M1.V.SDT_w$hitR[i] == 1){
            df.M1.V.SDT_w$hitR[i] <- 1 - 1/(2*(df.M1.V.SDT_w$hit[i] + df.M1.V.SDT_w$miss[i]))
      }
}

for (i in 1:nrow(df.M1.V.SDT_w)){
      if (df.M1.V.SDT_w$faR[i] == 0){
            df.M1.V.SDT_w$faR[i] <- 1/(2*(df.M1.V.SDT_w$FA[i] + df.M1.V.SDT_w$CR[i]))
      }
}

# calculate the d prime for each condition
df.M1.V.SDT_w$dprime <- mapply(dprime,df.M1.V.SDT_w$hitR,df.M1.V.SDT_w$faR)

# transfor from long format to wide format
df.M1.V.SDT_ww <- reshape2::dcast(df.M1.V.SDT_w, Subject + Age + Sex ~ Morality + Identity ,value.var = "dprime")
# rename the column number
colnames(df.M1.V.SDT_ww)[4:7] <- paste("d", colnames(df.M1.V.SDT_ww[,4:7]), sep = "_")

df.M1.V.SDT_l <- df.M1.V.SDT_w[,c(1:5,12)]

####################################
############      RT     ###########
####################################

df.M1.V.RT <- df.M1.V[df.M1.V$ACC == 1,]

# get the summary results
df.M1.V.RT.subj <- summarySEwithin(df.M1.V.RT,measurevar = 'RT', withinvar = c('Subject','Match','Morality','Identity'),idvar = 'Subject', na.rm = TRUE)
df.M1.V.RT.subj_w <- reshape2::dcast(df.M1.V.RT.subj, Subject ~ Match + Morality + Identity ,value.var = "RT") 

# rename the columns of RT data
colnames(df.M1.V.RT.subj_w)[2:9] <- paste("RT", colnames(df.M1.V.RT.subj_w[,2:9]), sep = "_")

## saving data ####
# merge the dprime and RT data and save (wide-format)
df.M1.V.sum_w <- merge(df.M1.V.acc_w,  df.M1.V.SDT_ww,by = "Subject")
df.M1.V.sum_w <- merge(df.M1.V.sum_w,df.M1.V.RT.subj_w,by = 'Subject')

# merge the RT and ACC data (long-format)
df.M1.V.sum_rt_acc_l <- merge(df.M1.V.acc,df.M1.V.RT.subj,by = c("Subject","Match","Morality",'Identity'))
df.M1.V.sum_rt_acc_l <- df.M1.V.sum_rt_acc_l[order(df.M1.V.sum_rt_acc_l$Subject),]

df.M1.V.sum_rt_acc_l <- df.M1.V.sum_rt_acc_l[,c("Subject","Match","Morality",'Identity',"N.x","countN","ACC","RT")]
colnames(df.M1.V.sum_rt_acc_l) <- c("Subject","Match","Morality",'Identity',"Ntrials","corrTrials","ACC","RT")

# order the columns
df.M1.V.sum_w <- df.M1.V.sum_w[,c(colnames(df.M1.V.sum_w)[c(1,10:11,2:9,12:23)])]

### plot the density of d prime and rt
# Create a scatterplot with density margin plots 
# Note, this part of the code is just for exploration, not presented in the manuscript

## prepare the data
df.M1.V.SDT_l$Cond <- paste(df.M1.V.SDT_l$Morality,df.M1.V.SDT_l$Identity,sep = '_')
df.M1.V.RT_m <- df.M1.V.sum_rt_acc_l %>%
   dplyr::filter(Match == 'match') %>%
   dplyr::mutate(Cond = paste(Morality, Identity, sep = '_')) %>%
   dplyr::mutate(RT = RT*1000)

df.M1.V_dens_p <- merge(df.M1.V.SDT_l, df.M1.V.RT_m)

# The plotHolder() function from C-3PR creates a blank plot template that will hold the figures
plotHolder <- function(useArial = F,afmPATH="~/Dropbox"){
   require(ggplot2)
   ggplot() +
      geom_blank(aes(1,1)) +
      theme(line = element_blank(),
            text  = element_blank(),
            title = element_blank(),
            plot.background = element_blank(),
            #           panel.grid.major = element_blank(),
            #           panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()
            #           axis.title.x = element_blank(),
            #           axis.title.y = element_blank(),
            #           axis.text.x = element_blank(),
            #           axis.text.y = element_blank(),
            #           axis.ticks = element_blank()
      )
}

gg.theme <- function(type=c("clean","noax")[1],useArial = F, afmPATH="~/Dropbox"){
   require(ggplot2)
   if(useArial){
      set.Arial(afmPATH)
      bf_font="Arial"
   } else {bf_font="Helvetica"}
   
   switch(type,
          clean = theme_bw(base_size = 16, base_family=bf_font) +
             theme(axis.text.x     = element_text(size = 14),
                   axis.title.y    = element_text(vjust = +1.5),
                   panel.grid.major  = element_blank(),
                   panel.grid.minor  = element_blank(),
                   legend.background = element_blank(),
                   legend.key = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.line  = element_line(colour = "black")),
          
          noax = theme(line = element_blank(),
                       text  = element_blank(),
                       title = element_blank(),
                       plot.background = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank())
   )
}


blankPlot <- plotHolder()

# X margin density plot (note: gg.theme() from C-3PR can be used directly in a ggplot2() call)

xDense <- ggplot(df.M1.V_dens_p, aes(x=dprime, fill=Cond)) + 
   geom_density(aes(y= ..count..),trim=F,alpha=.5) + 
   xlab("") + ylab("") + xlim(0,5) +
   gg.theme("noax") + 
   theme(legend.position = "none",plot.margin = unit(c(0,0,0,4), "lines"))

## Uncomment to save subplot
# ggsave("RPP_F3_xDense.png",plot=xDense)

# Y margin density plot (note: gg.theme() from C-3PR can be used directly in a ggplot2() call)

yDense <- ggplot(df.M1.V_dens_p, aes(x=RT, fill=Cond)) + 
   geom_density(aes(y= ..count..),trim=F,alpha=.5) + 
   xlab("") + ylab("") + xlim(400,900) + 
   coord_flip() + 
   gg.theme("noax") + 
   theme(legend.position = "none", plot.margin = unit(c(0,0,3,0), "lines")) 

## Uncomment to save subplot
# ggsave("RPP_F3_yDense.png",plot=yDense)

# The main scatterplot (note: gg.theme() from C-3PR can be used directly in a ggplot2() call)
scatterP<-
   ggplot(df.M1.V_dens_p,aes(x=dprime,y=RT)) +  
   #geom_hline(aes(yintercept=0),linetype=2) +
   #geom_abline(intercept=0,slope=1,color="Grey60")+
   geom_point(aes(size= 0.8,fill=Cond),color="Grey30",shape=21,alpha=.8) + 
   #geom_rug(aes(color=Cond),size=1,sides="b",alpha=.6) + 
   #geom_rug(aes(color=Cond),,size=1,sides="l",alpha=.6) + 
   scale_x_continuous(name="d prime",limits=c(0,5),breaks=c(0,1,2,3,4,5)) + 
   scale_y_continuous(name="RT",limits=c(400,900),breaks=c(400,500,600,700,800,900)) + 
   ggtitle("") + xlab("") + ylab("") + 
   #scale_size_continuous(name="Replication Power",range=c(2,9)) + 
   scale_color_discrete(name='Conditions') +
   #scale_fill_discrete(name="p-value") +
   gg.theme("clean") + 
   theme(legend.position=c(.9,.6), plot.margin = unit(c(-2,-1.5,2,2), "lines")) 

## Uncomment to save subplot
# ggsave("RPP_F3_scatter.png",plot=scatterP)

# Yet another way to organise plots: grid.arrange() from the gridExtra package.
gridExtra::grid.arrange(xDense, blankPlot, scatterP, yDense, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

# SAVE combined plots as PDF
pdf("Figure_RT_drpime_density_explore.pdf",pagecentre=T, width=15,height=12 ,paper = "special")
gridExtra::grid.arrange(xDense, blankPlot, scatterP, yDense, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()


####################################
#####  effects in matching  ########
####################################
# calculate the effect of self-ref and valence
df.M1.v.sum_eff_w <- df.M1.V.sum_w[,1:3]
df.M1.v.sum_eff_w$d_goodslf_goodoth <- df.M1.V.sum_w$d_Good_Self - df.M1.V.sum_w$d_Good_Other
df.M1.v.sum_eff_w$d_goodslf_badslf  <- df.M1.V.sum_w$d_Good_Self - df.M1.V.sum_w$d_Bad_Self
df.M1.v.sum_eff_w$d_goodoth_badoth  <- df.M1.V.sum_w$d_Good_Other - df.M1.V.sum_w$d_Bad_Other
df.M1.v.sum_eff_w$d_badslf_badoth   <- df.M1.V.sum_w$d_Bad_Self - df.M1.V.sum_w$d_Bad_Other

df.M1.v.sum_eff_w$RT_goodslf_goodoth <- df.M1.V.sum_w$RT_match_Good_Other -  df.M1.V.sum_w$RT_match_Good_Self
df.M1.v.sum_eff_w$RT_goodslf_badslf  <- df.M1.V.sum_w$RT_match_Bad_Self -  df.M1.V.sum_w$RT_match_Good_Self
df.M1.v.sum_eff_w$RT_goodoth_badoth  <- df.M1.V.sum_w$RT_match_Bad_Other -  df.M1.V.sum_w$RT_match_Good_Other
df.M1.v.sum_eff_w$RT_badslf_badoth   <- df.M1.V.sum_w$RT_match_Bad_Self -  df.M1.V.sum_w$RT_match_Bad_Other

# write files
setwd(traDir)
write.csv(df.M1.V.sum_w,'MS_match_behav_wide.csv',row.names = F)
write.csv(df.M1.V.SDT_l,'MS_match__dprime_long.csv',row.names = F)
write.csv(df.M1.V.sum_rt_acc_l,'MS_match__rt_acc_long.csv',row.names = F)
setwd(curDir)


# ---------------------------------------------------------------------------------------
# ---------- 4. Categorization task: prepare the Accuracy and RT ------------------------
# ---------------------------------------------------------------------------------------

####################################
###############   ACC    ###########
####################################

df.C1.V.acc <- plyr::ddply(df.C1.V,.(Subject,Age, Sex, Task,Morality,Identity), summarise,
                           N = length(ACC),
                           countN = sum(ACC),
                           ACC = sum(ACC)/length(ACC))

# wide-format
df.C1.V.acc_w  <- reshape2::dcast(df.C1.V.acc, Subject ~ Task + Morality + Identity ,value.var = "ACC")
# rename the column number
colnames(df.C1.V.acc_w)[2:9] <- paste("ACC", colnames(df.C1.V.acc_w[,2:9]), sep = "_")

# combing data from diff task for analyzing the interaction btw val and id
df.C1.V.acc_noTask  <-  plyr::ddply(df.C1.V,.(Subject, Morality, Identity), summarise,
                                    N = length(ACC),
                                    countN = sum(ACC),
                                    ACC = sum(ACC)/length(ACC))

df.C1.V.acc_noTask_w <- reshape2::dcast(df.C1.V.acc_noTask, Subject ~ Morality + Identity,value.var = "ACC")
# rename the column number
colnames(df.C1.V.acc_noTask_w)[2:5] <- paste("ACC", colnames(df.C1.V.acc_noTask_w[,2:5]), sep = "_")

# calculate the mean differences(?)
df.C1.V.acc_noTask_w$good_bad_slf <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Bad_Self
df.C1.V.acc_noTask_w$good_bad_oth <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Bad_Other
df.C1.V.acc_noTask_w$slf_oth_good <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Good_Other
df.C1.V.acc_noTask_w$slf_oth_bad  <- df.C1.V.acc_noTask_w$ACC_Bad_Self - df.C1.V.acc_noTask_w$ACC_Bad_Other

####################################
###############   RT     ###########
####################################

df.C1.V.RT <- df.C1.V[df.C1.V$ACC == 1,]  # exclued inaccurate data
df.C1.V.RT.subj <- summarySEwithin(df.C1.V.RT,measurevar = 'RT', withinvar = c('Subject','Task','Morality','Identity'), idvar = 'Subject',na.rm = TRUE)
df.C1.V.RT.subj_w <- reshape2::dcast(df.C1.V.RT.subj, Subject ~ Task + Morality + Identity ,value.var = "RT") 

# rename the columns of RT data
colnames(df.C1.V.RT.subj_w)[2:9] <- paste("RT", colnames(df.C1.V.RT.subj_w[,2:9]), sep = "_")

# combining data form different task for analyszing interaction of val and id
df.C1.V.RT.subj_noTask <- summarySEwithin(df.C1.V.RT,measurevar = 'RT', withinvar = c('Subject','Morality','Identity'), idvar = 'Subject',na.rm = TRUE)
df.C1.V.RT.subj_noTask_w <- reshape2::dcast(df.C1.V.RT.subj_noTask, Subject ~ Morality + Identity ,value.var = "RT") 

# rename the columns of RT data
colnames(df.C1.V.RT.subj_noTask_w)[2:5] <- paste("RT", colnames(df.C1.V.RT.subj_noTask_w[,2:5]), sep = "_")

## saving data ####
# merge the accuracy and RT data and save
df.C1.V.sum_w <- merge(df.C1.V.acc_w,  df.C1.V.RT.subj_w,by = "Subject")
df.C1.V.sum_noTask_w <- merge(df.C1.V.acc_noTask_w,  df.C1.V.RT.subj_noTask_w,by = "Subject")

####################################
## effects in categorization  ######
####################################
# calculate the effect of self-ref and valence
df.C1.v.sum_eff_w <- data.frame(df.C1.V.sum_w[,c('Subject')])
colnames(df.C1.v.sum_eff_w) <- 'Subject'
df.C1.v.sum_eff_w$Val_RT_goodslf_goodoth <- df.C1.V.sum_w$RT_Val_Good_Other - df.C1.V.sum_w$RT_Val_Good_Self
df.C1.v.sum_eff_w$Val_RT_goodslf_badslf  <- df.C1.V.sum_w$RT_Val_Bad_Self - df.C1.V.sum_w$RT_Val_Good_Self
df.C1.v.sum_eff_w$Val_RT_goodoth_badoth  <- df.C1.V.sum_w$RT_Val_Bad_Other - df.C1.V.sum_w$RT_Val_Good_Other
df.C1.v.sum_eff_w$Val_RT_badslf_badoth   <- df.C1.V.sum_w$RT_Val_Bad_Self - df.C1.V.sum_w$RT_Val_Bad_Other

df.C1.v.sum_eff_w$Id_RT_goodslf_goodoth  <- df.C1.V.sum_w$RT_Id_Good_Other - df.C1.V.sum_w$RT_Id_Good_Self
df.C1.v.sum_eff_w$Id_RT_goodslf_badslf   <- df.C1.V.sum_w$RT_Id_Bad_Self - df.C1.V.sum_w$RT_Id_Good_Self
df.C1.v.sum_eff_w$Id_RT_goodoth_badoth   <- df.C1.V.sum_w$RT_Id_Bad_Other - df.C1.V.sum_w$RT_Id_Good_Other
df.C1.v.sum_eff_w$Id_RT_badslf_badoth    <- df.C1.V.sum_w$RT_Id_Bad_Self - df.C1.V.sum_w$RT_Id_Bad_Other

df.C1.v.sum_eff_w$Val_ACC_goodslf_goodoth <- df.C1.V.sum_w$ACC_Val_Good_Self - df.C1.V.sum_w$ACC_Val_Good_Other
df.C1.v.sum_eff_w$Val_ACC_goodslf_badslf  <- df.C1.V.sum_w$ACC_Val_Good_Self - df.C1.V.sum_w$ACC_Val_Bad_Self
df.C1.v.sum_eff_w$Val_ACC_goodoth_badoth  <- df.C1.V.sum_w$ACC_Val_Good_Other - df.C1.V.sum_w$ACC_Val_Bad_Other
df.C1.v.sum_eff_w$Val_ACC_badslf_badoth   <- df.C1.V.sum_w$ACC_Val_Bad_Self - df.C1.V.sum_w$ACC_Val_Bad_Other

df.C1.v.sum_eff_w$Id_ACC_goodslf_goodoth  <- df.C1.V.sum_w$ACC_Id_Good_Self - df.C1.V.sum_w$ACC_Id_Good_Other
df.C1.v.sum_eff_w$Id_ACC_goodslf_badslf   <- df.C1.V.sum_w$ACC_Id_Good_Self - df.C1.V.sum_w$ACC_Id_Bad_Self
df.C1.v.sum_eff_w$Id_ACC_goodoth_badoth   <- df.C1.V.sum_w$ACC_Id_Good_Other - df.C1.V.sum_w$ACC_Id_Bad_Other
df.C1.v.sum_eff_w$Id_ACC_badslf_badoth    <- df.C1.V.sum_w$ACC_Id_Bad_Self - df.C1.V.sum_w$ACC_Id_Bad_Other

# merge the effect file
df.v.sum_eff_all_w <- merge(df.M1.v.sum_eff_w,df.C1.v.sum_eff_w,by="Subject")

# merge the RT and ACC data (long-format) ####
df.C1.V.sum_rt_acc_l <- merge(df.C1.V.acc,df.C1.V.RT.subj,by = c("Subject","Task","Morality",'Identity'))
df.C1.V.sum_rt_acc_l <- df.C1.V.sum_rt_acc_l[order(df.C1.V.sum_rt_acc_l$Subject),]

df.C1.V.sum_rt_acc_l <- df.C1.V.sum_rt_acc_l[,c("Subject","Task","Morality",'Identity',"N.x","countN","ACC","RT")]
colnames(df.C1.V.sum_rt_acc_l) <- c("Subject","Task","Morality",'Identity',"Ntrials","corrTrials","ACC","RT")

# merge the RT and ACC data without task (long-format) ####
df.C1.V.sum_rt_acc_noTask_l <- merge(df.C1.V.acc_noTask,df.C1.V.RT.subj_noTask,by = c("Subject","Morality",'Identity'))
df.C1.V.sum_rt_acc_noTask_l <- df.C1.V.sum_rt_acc_noTask_l[order(df.C1.V.sum_rt_acc_noTask_l$Subject),]

df.C1.V.sum_rt_acc_noTask_l <- df.C1.V.sum_rt_acc_noTask_l[,c("Subject","Morality",'Identity',"N.x","countN","ACC","RT")]
colnames(df.C1.V.sum_rt_acc_noTask_l) <- c("Subject","Morality",'Identity',"Ntrials","corrTrials","ACC","RT")

# write files to an upper-level folder
setwd(traDir)
write.csv(df.C1.V.sum_w,'MS_categ_behav_wide.csv',row.names = F)
write.csv(df.C1.V.sum_noTask_w,'MS_categ_behav_noTask_wide.csv',row.names = F)
write.csv(df.C1.V.sum_rt_acc_l,'MS_categ__rt_acc_long.csv',row.names = F)
write.csv(df.C1.V.sum_rt_acc_noTask_l,'MS_categ__rt_acc_noTask_long.csv',row.names = F)
write.csv(df.v.sum_eff_all_w,'MS_cross_taskeffect_wide.csv',row.names = F)
setwd(curDir)


# ---------------------------------------------------------------------------------------
# ------------ 5. Plots  ----------------------------------------------------------------
# ---------------------------------------------------------------------------------------
# plots here are made by pre-defined functions in initial.r:
#     "Mplots"  -- plot for matching task
#     "CAplots" -- plot for categorization task
## Matching task, plots are save to 'saveDir'

# first convert time to ms
df.M1.V.sum_rt_acc_l$RT <- df.M1.V.sum_rt_acc_l$RT*1000
df.C1.V.sum_rt_acc_l$RT <- df.C1.V.sum_rt_acc_l$RT*1000
df.C1.V.sum_rt_acc_noTask_l$RT <- df.C1.V.sum_rt_acc_noTask_l$RT*1000

Mplots(saveDir = traDir, curDir = curDir, expName = 'conf', df.M1.V.SDT_l,df.M1.V.sum_rt_acc_l)

# plot id-based data
CAplots(saveDir = traDir, curDir = curDir,expName = 'conf', task = 'id', df.C1.V.sum_rt_acc_l)

# plot val-based data
CAplots(saveDir = traDir, curDir = curDir,expName = 'conf', task = 'val', df.C1.V.sum_rt_acc_l)

# plot the categorization task (collapsed different tasks)
CAplots(saveDir = traDir, curDir = curDir,expName = 'conf', task = 'categ', df.C1.V.sum_rt_acc_noTask_l)

