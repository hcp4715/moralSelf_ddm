### This script is used to preprocessing the data of the moral self experiment (DDM), experiment 1. 
#
### Hu, C-P., Lan, Y., Macrae, N., & Sui. J. (2019) 
### Script author: Chuan-Peng Hu
### Email = hcp4715@gmail.com       twitter= @hcp4715
#
### Input: 
###      MS_matchTask_raw.csv -- raw data of matching task;
###      MS_categTask_raw.csv -- raw data of categorization task;
#
### Output:
###     df.M.hddm_m.csv  -- clean data of matching trials from matchign task, for HDDM analysis;
###     df.M.hddm_nm.csv -- clean data of nonmatching trials from matchign task, for HDDM analysis;
###     df.C.hddm_val    -- clean data of valence-based categorization trials, for HDDM analysis;
###     df.C.hddm_id     -- clean data of identity-based categorization trials, for HDDM analysis;
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
# ---------- 1. Initializing and prepare  -----------------------------------------------
# ---------------------------------------------------------------------------------------

curDir  <- dirname(rstudioapi::getSourceEditorContext()$path)   # get the directory for preprocessing
setwd(curDir)
source('Initial_exp7_pilot.r')  # initializing (clear global environment; load packages and functions)
curDir  <- dirname(rstudioapi::getSourceEditorContext()$path)   # get the directory for preprocessing
rootDir <- gsub('.{9}$', '', curDir)                            # get the parental folder, remove last 9 characters
traDir <- paste(rootDir,'2_trad_analysis',sep = '')             # folder for traditional analysis
ddmDir <- paste(rootDir,'4_hddm',sep = '')                      # folder for DDM analysis
exgDir <- paste(rootDir,'3_exGaussian',sep = '')                # folder for ExGaussian analysis

# ---------------------------------------------------------------------------------------
# ---------- 2. Loading data and clean the data  ----------------------------------------
# ---------------------------------------------------------------------------------------

df.M.raw <- read.csv("MS_matchTask_raw.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE) # data for matching task
df.C.raw <- read.csv("MS_categTask_raw.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE) # data for categorization task

#subTotl <- unique(df.M.raw$Subject)

# make the variables in a specified order
df.M.raw <- df.M.raw %>%
   dplyr::rename(ACC = Accuracy, Subject = SubjectID) %>%                    # rename to ACC   
   dplyr::mutate(Morality = ifelse(Morality == "moral", 'Good','Bad'),       # change condition name
                 Identity = ifelse(Identity == "self", 'Self','Other'),      # change condition name
                 Morality = factor(Morality, levels = c("Good","Bad")),      # to factor
                 Identity = factor(Identity, levels=c("Self","Other")),
                 Match    = factor(Match,    levels=c("match","nonmatch")))

df.C.raw <- df.C.raw %>%
   dplyr::filter(Task == 'self' | Task == 'moral') %>%      # exclude trials from importance-based task
   dplyr::rename(ACC = Accuracy, Subject = SubjectID) %>%                    # rename to ACC   
   dplyr::mutate(Morality = ifelse(Morality == "moral", 'Good','Bad'),       # change condition name
                 Identity = ifelse(Identity == "self", 'Self','Other'),      # change condition name
                 Task     = ifelse(Task     == "self", 'Id','Val'),          # change condition name
                 Morality = factor(Morality, levels=c("Good","Bad")),      # to factor
                 Identity = factor(Identity, levels=c("Self","Other")))

# Exclude Subject, criterion 1: procedure failure
excldSub1 <- c("7","8","2027","7035")
df.M <- df.M.raw[!(df.M.raw$Subject %in% excldSub1),]
df.C <- df.C.raw[!(df.C.raw$Subject %in% excldSub1),]

# exclude trials, criterio 2: practicing trials in matching task (first 48 trials)
subNo <- unique(df.M$Subject)
for (subj in subNo) {
      if (exists('df.M.fm')){
            df.tmp <- df.M[df.M$Subject == subj,]
            df.tmp <- df.tmp[49:nrow(df.tmp),]
            df.M.fm <- rbind(df.M.fm,df.tmp)
      } else {
            df.M.fm <- df.M[df.M$Subject == subj,]
            df.M.fm <- df.M.fm[49:nrow(df.M.fm),]
      }
}   # df.M.fm should be 14880 rows

df.M <- df.M.fm          # all the experimental trials
rm(subNo,df.M.fm,df.tmp) # remove the intermedia variables

# exlude subject criterion 2: less than 50% overall accuracy
# re-code trials without response as incorrect
## NOTE: this was not done in previous analysis, 
##       resulting the wronly excluded participants in registration
df.M.tmp <- df.M
df.M.tmp$ACC[df.M.tmp$ACC == -1] <- 0
df.C.tmp <- df.C
df.C.tmp$ACC[df.C.tmp$ACC == -1] <- 0 

# calculate the overall accuracy for matching task
df.M.acc.g <-  df.M.tmp %>%
      dplyr::group_by(Subject) %>%
      dplyr::summarise(N = length(ACC),
                     countN = sum(ACC),
                     ACC = sum(ACC)/length(ACC)) %>%
      dplyr::ungroup()

# calculate the overall accuracy for categorziation task
df.C.acc.g <-  df.C.tmp  %>%
      dplyr::group_by(Subject) %>%
      dplyr::summarise(N = length(ACC),
                       countN = sum(ACC),
                       ACC = sum(ACC)/length(ACC)) %>%
      dplyr::ungroup()

excldSub2.M <- df.M.acc.g$Subject[df.M.acc.g$ACC < 0.5]      # < 50% accuracy in matching task
excldSub2.C <- df.C.acc.g$Subject[df.C.acc.g$ACC < 0.5]      # < 50% accuracy in categorization task
df.M.valid <- df.M[!(df.M$Subject %in% excldSub2.M),]        # exclude the invalid subjects
df.C.valid <- df.C[!(df.C$Subject %in% excldSub2.M),]
df.C.valid <- df.C.valid[!(df.C.valid$Subject %in% excldSub2.C),] # one participants were excluded for categorization task
rm(df.M,df.C,df.M.tmp,df.C.tmp)

# clean data for traditional analysis
df.M1 <- df.M.valid          # percentage of trials without response: 0.039799
df.C1 <- df.C.valid          # percentage of trials without response: 0.029514
df.M1$RT <- df.M1$RT * 1000  # transfer from seconds to min seconds
df.C1$RT <- df.C1$RT * 1000  # careful about the scale of time when using HDDM or ex-Gaussian

# no response equal to wrong
sum(df.M1$ACC == -1)/length(df.M1$ACC)
sum(df.C1$ACC == -1)/length(df.C1$ACC)
df.M1$ACC[df.M1$ACC == -1] <- 0
df.C1$ACC[df.C1$ACC == -1] <- 0

# excld.trials, criterion 3: correct response within 200 ms
excld.trials.M <- df.M1[df.M1$RT <= 200,]
df.M1.V        <- df.M1[!(df.M1$RT <= 200),]        # valid trial data for match task
excld.trials.C <- df.C1[df.C1$RT <= 200,]
df.C1.V        <- df.C1[!(df.C1$RT <= 200),]        # valid trial data for categorization task
nrow(excld.trials.C) + nrow(df.C1.V) == nrow(df.C1) # if true, everything is ok

## Basic information of the data ####
df.M1.T.basic.info <- df.M.raw %>%
      dplyr::select(Subject, Age, Gender) %>%
      dplyr::distinct(Subject, Age, Gender) %>%
      dplyr::summarise(N = length(Subject),
                       N.f = sum(Gender == 'female'),
                       N.m = sum(Gender == 'male'),
                       MeanAge = round(mean(Age),2),
                       SDAge   = round(sd(Age),2))

num.excldSub2.M <- length(unique(excldSub2.M))
num.excldSub2.C <- length(unique(excldSub2.C))

# valid data for matching task
df.M1.V.basic.info <- df.M1.V %>%
      dplyr::select(Subject, Age, Gender) %>%
      dplyr::distinct(Subject, Age, Gender) %>%
      dplyr::summarise(N = length(Subject),
                       N.f = sum(Gender == 'female'),
                       N.m = sum(Gender == 'male'),
                       MeanAge = round(mean(Age),2),
                       SDAge   = round(sd(Age),2))

ratio.excld.trials.M <- nrow(excld.trials.M)/nrow(df.M1)

# valid data for categorization task
df.C1.V.basic.info <- df.C1.V %>%
      dplyr::select(Subject, Age, Gender) %>%
      dplyr::distinct(Subject, Age, Gender) %>%
      dplyr::summarise(N = length(Subject),
                       N.f = sum(Gender == 'female'),
                       N.m = sum(Gender == 'male'),
                       MeanAge = round(mean(Age),2),
                       SDAge   = round(sd(Age),2))

ratio.excld.trials.C <- nrow(excld.trials.C)/nrow(df.C1)

###########################################################
########   get the data file for hddm analysis      #######
###########################################################
# exclusion trials, criterion: trials without response or wrong key
# Note: we didn't remove the correct response less than 200 ms

df.M.hddm_m <- df.M.valid %>%                            # start from the valid data, RT in seconds.
   dplyr::filter(ACC == 1 | ACC == 0) %>%                # exclude trials without response or with wrong keys
   dplyr::filter(Match == 'match') %>%
   dplyr::select(Subject,Morality,Identity,RT,ACC) %>%   # select column
   dplyr::rename(subj_idx = Subject, val = Morality,     # rename column
                 id = Identity, rt = RT, response = ACC) # rename column


df.M.hddm_nm <- df.M.valid %>%
   dplyr::filter(ACC == 1 | ACC == 0) %>%                # exclude trials without response or with wrong keys
   dplyr::filter(Match == 'nonmatch') %>%
   dplyr::select(Subject,Morality,Identity,RT,ACC) %>%   # select column
   dplyr::rename(subj_idx = Subject, val = Morality, 
                 id = Identity, rt = RT, response = ACC) # rename column

##################################################
### preparing the data for Stimuli-code hddm  ####
##################################################
# for valence based, stim: moral = 1, immoral = 0; response: moralKey = 1; immoralKey = 0
df.C.hddm_val_stim <- df.C.valid %>%
   dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
   dplyr::filter(Task == 'Val') %>%                     # select the valence-based task
   dplyr::mutate(response = ifelse((trialType == "moral" & ACC ==1) | (trialType == "immoral" & ACC ==0), 1, 0)) %>%
   dplyr::select(Subject,Morality,Identity,trialType, response, RT) %>%  # select column
   dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, stim = trialType) %>% # rename column
   dplyr::mutate(stim = ifelse(stim == "moral",1, 0))

# for Id based, stim: self = 1, other = 0; response: selfKey = 1; otherKey = 0
df.C.hddm_id_stim <- df.C.valid %>%
   dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
   dplyr::filter(Task == 'Id') %>%                      # select the valence-based task
   dplyr::mutate(response = ifelse((trialType=="self" & ACC==1) | (trialType=="other" & ACC==0), 1, 0)) %>%
   dplyr::select(Subject,Morality,Identity,trialType, response, RT) %>%  # select column
   dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, stim = trialType) %>% # rename column
   dplyr::mutate(stim = ifelse(stim == "self",1, 0))

setwd(ddmDir)
write.csv(df.M.hddm_m,'MS_match_hddm.csv',row.names = F)
write.csv(df.M.hddm_nm,'MS_mismatch_hddm.csv',row.names = F)
write.csv(df.C.hddm_val_stim,'MS_categ_val_hddm_stim.csv',row.names = F)
write.csv(df.C.hddm_id_stim,'MS_categ_id_hddm_stim.csv',row.names = F)
setwd(curDir)
rm('df.M.hddm_m','df.M.hddm_nm','df.C.hddm_val_stim','df.C.hddm_id_stim')

###########################################################
######## get the data file for exGaussian analysis  #######
###########################################################
df.C.exG <- df.C.valid %>%
      dplyr::filter(ACC == 1 & RT > .2)    # only for correct and valid trials

setwd(exgDir)
write.csv(df.C.exG,'MS_categ_exG.csv',row.names = F)
rm(df.C.exG)

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
                    ACC = sum(ACC)/length(ACC)) %>%
   dplyr::ungroup()

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
      } else if (df.M1.V$ACC[i] == 1 & df.M1.V$Match[i] == "nonmatch" ){
            df.M1.V$sdt[i] <- "CR"
      } else if (df.M1.V$ACC[i] == 0 & df.M1.V$Match[i] == "match"){
            df.M1.V$sdt[i] <- "miss"
      } else if (df.M1.V$ACC[i] == 0 & df.M1.V$Match[i] == "nonmatch" ){
            df.M1.V$sdt[i] <- "FA"
      }
}

# calculate the number of each for each condition
#df.M1.V.SDT <- plyr::ddply(df.M1.V,.(Subject,Age, Gender, Morality, Identity,sdt), summarise, N = length(sdt))
df.M1.V.SDT <- df.M1.V %>%
   dplyr::group_by(Subject,Age, Gender, Morality, Identity,sdt) %>%
   dplyr::summarise(N = length(sdt)) %>%
   dplyr::ungroup()
df.M1.V.SDT <- df.M1.V.SDT[complete.cases(df.M1.V.SDT$sdt),] # exclude NA, which represents no-response trials

# long format to wide
df.M1.V.SDT_w <- reshape2::dcast(df.M1.V.SDT, Subject + Age + Gender + Morality + Identity ~ sdt,value.var = "N")
df.M1.V.SDT_w$miss[is.na(df.M1.V.SDT_w$miss)] <- 0 # transfer NA to 0
df.M1.V.SDT_w$FA[is.na(df.M1.V.SDT_w$FA)]     <- 0
df.M1.V.SDT_w$hitR <- df.M1.V.SDT_w$hit/(df.M1.V.SDT_w$hit + df.M1.V.SDT_w$miss)
df.M1.V.SDT_w$faR <- df.M1.V.SDT_w$FA/(df.M1.V.SDT_w$FA + df.M1.V.SDT_w$CR)

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
df.M1.V.SDT_w$dprime <- mapply(dprime, df.M1.V.SDT_w$hitR, df.M1.V.SDT_w$faR)

# transfor from long format to wide format
df.M1.V.SDT_ww <- reshape2::dcast(df.M1.V.SDT_w, Subject + Age + Gender ~ Morality + Identity ,value.var = "dprime")

# rename the column number
colnames(df.M1.V.SDT_ww)[4:7] <- paste("d", colnames(df.M1.V.SDT_ww[,4:7]), sep = "_")

df.M1.V.SDT_l <- df.M1.V.SDT_w[,c(1:5,12)] # re-order columns

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
# merge the dprime and RT data and save
df.M1.V.sum_w <- merge(df.M1.V.acc_w,  df.M1.V.SDT_ww,by = "Subject")
df.M1.V.sum_w <- merge(df.M1.V.sum_w,df.M1.V.RT.subj_w,by = 'Subject')

# merge the RT and ACC data (long-format)
df.M1.V.sum_rt_acc_l <- merge(df.M1.V.acc,df.M1.V.RT.subj,by = c("Subject","Match","Morality",'Identity'))
df.M1.V.sum_rt_acc_l <- df.M1.V.sum_rt_acc_l[order(df.M1.V.sum_rt_acc_l$Subject),]

df.M1.V.sum_rt_acc_l <- df.M1.V.sum_rt_acc_l[,c("Subject","Match","Morality",'Identity',"N.x","countN","ACC","RT")]
colnames(df.M1.V.sum_rt_acc_l) <- c("Subject","Match","Morality",'Identity',"Ntrials","corrTrials","ACC","RT")

# order the columns
df.M1.V.sum_w <- df.M1.V.sum_w[,c(colnames(df.M1.V.sum_w)[c(1,10:11,2:9,12:23)])]

####################################
#####  effects in matching  ########
####################################
# calculate the effect of self-ref and valence
df.M1.v.sum_eff_w <- df.M1.V.sum_w %>%
   dplyr::mutate(d_goodslf_goodoth = d_Good_Self - d_Good_Other,
                 d_goodslf_badslf  = d_Good_Self - d_Bad_Self,
                 d_goodoth_badoth  = d_Good_Other - d_Bad_Other,
                 d_badslf_badoth   = d_Bad_Self - d_Bad_Other,
                 RT_goodslf_goodoth = RT_match_Good_Other - RT_match_Good_Self,
                 RT_goodslf_badslf  = RT_match_Bad_Self -  RT_match_Good_Self,
                 RT_goodoth_badoth  = RT_match_Bad_Other - RT_match_Good_Other,
                 RT_badslf_badoth   = RT_match_Bad_Self -  RT_match_Bad_Other) %>%
   dplyr::select(Subject, Age, Gender,
                 d_goodslf_goodoth,  d_goodslf_badslf,  d_goodoth_badoth,  d_badslf_badoth,
                 RT_goodslf_goodoth, RT_goodslf_badslf, RT_goodoth_badoth, RT_badslf_badoth)

# write files
setwd(traDir)
write.csv(df.M1.V.sum_w,'MS_match_behav_wide.csv',row.names = F)
write.csv(df.M1.V.SDT_l,'MS_match__dprime_long.csv',row.names = F)
# write.csv(df.M1.V.sum_rt_acc_l,'MS_match__rt_acc_long.csv',row.names = F) # long-format data
setwd(curDir)


# ---------------------------------------------------------------------------------------
# ---------- 4. Categorization task: prepare the Accuracy and RT ------------------------
# ---------------------------------------------------------------------------------------

####################################
###############   ACC    ###########
####################################

df.C1.V.acc <- plyr::ddply(df.C1.V,.(Subject,Age, Gender, Task,Morality,Identity), summarise,
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

# calculate the mean differences, might be useful to correlation with questionnarie data
df.C1.V.acc_noTask_w$good_bad_slf <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Bad_Self
df.C1.V.acc_noTask_w$good_bad_oth <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Bad_Other
df.C1.V.acc_noTask_w$slf_oth_good <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Good_Other
df.C1.V.acc_noTask_w$slf_oth_bad  <- df.C1.V.acc_noTask_w$ACC_Bad_Self - df.C1.V.acc_noTask_w$ACC_Bad_Other

####################################
###############   RT     ###########
####################################

df.C1.V.RT <- df.C1.V[df.C1.V$ACC == 1,]  # exclued inaccurate data
df.C1.V.RT.subj <- summarySEwithin(df.C1.V.RT,measurevar = 'RT', withinvar = c('Subject','Task','Morality','Identity'), idvar = 'Subject',na.rm = TRUE)
df.C1.V.RT.subj_w <- dcast(df.C1.V.RT.subj, Subject ~ Task + Morality + Identity ,value.var = "RT") 

# rename the columns of RT data
colnames(df.C1.V.RT.subj_w)[2:9] <- paste("RT", colnames(df.C1.V.RT.subj_w[,2:9]), sep = "_")

# combining data form different task for analyszing interaction of val and id
df.C1.V.RT.subj_noTask <- summarySEwithin(df.C1.V.RT,measurevar = 'RT', withinvar = c('Subject','Morality','Identity'), idvar = 'Subject',na.rm = TRUE)
df.C1.V.RT.subj_noTask_w <- dcast(df.C1.V.RT.subj_noTask, Subject ~ Morality + Identity ,value.var = "RT") 

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
df.C1.v.sum_eff_w$Val_RT_goodoth_badoth <- df.C1.V.sum_w$RT_Val_Bad_Other - df.C1.V.sum_w$RT_Val_Good_Other
df.C1.v.sum_eff_w$Val_RT_badslf_badoth <- df.C1.V.sum_w$RT_Val_Bad_Self - df.C1.V.sum_w$RT_Val_Bad_Other

df.C1.v.sum_eff_w$Id_RT_goodslf_goodoth  <- df.C1.V.sum_w$RT_Id_Good_Other - df.C1.V.sum_w$RT_Id_Good_Self
df.C1.v.sum_eff_w$Id_RT_goodslf_badslf   <- df.C1.V.sum_w$RT_Id_Bad_Self - df.C1.V.sum_w$RT_Id_Good_Self
df.C1.v.sum_eff_w$Id_RT_goodoth_badoth <- df.C1.V.sum_w$RT_Id_Bad_Other - df.C1.V.sum_w$RT_Id_Good_Other
df.C1.v.sum_eff_w$Id_RT_badslf_badoth <- df.C1.V.sum_w$RT_Id_Bad_Self - df.C1.V.sum_w$RT_Id_Bad_Other

df.C1.v.sum_eff_w$Val_ACC_goodslf_goodoth <- df.C1.V.sum_w$ACC_Val_Good_Self - df.C1.V.sum_w$ACC_Val_Good_Other
df.C1.v.sum_eff_w$Val_ACC_goodslf_badslf <- df.C1.V.sum_w$ACC_Val_Good_Self - df.C1.V.sum_w$ACC_Val_Bad_Self
df.C1.v.sum_eff_w$Val_ACC_goodoth_badoth <- df.C1.V.sum_w$ACC_Val_Good_Other - df.C1.V.sum_w$ACC_Val_Bad_Other
df.C1.v.sum_eff_w$Val_ACC_badslf_badoth <- df.C1.V.sum_w$ACC_Val_Bad_Self - df.C1.V.sum_w$ACC_Val_Bad_Other

df.C1.v.sum_eff_w$Id_ACC_goodslf_goodoth  <- df.C1.V.sum_w$ACC_Id_Good_Self - df.C1.V.sum_w$ACC_Id_Good_Other
df.C1.v.sum_eff_w$Id_ACC_goodslf_badslf   <- df.C1.V.sum_w$ACC_Id_Good_Self - df.C1.V.sum_w$ACC_Id_Bad_Self
df.C1.v.sum_eff_w$Id_ACC_goodoth_badoth <- df.C1.V.sum_w$ACC_Id_Good_Other - df.C1.V.sum_w$ACC_Id_Bad_Other
df.C1.v.sum_eff_w$Id_ACC_badslf_badoth <- df.C1.V.sum_w$ACC_Id_Bad_Self - df.C1.V.sum_w$ACC_Id_Bad_Other

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
#write.csv(df.C1.V.sum_rt_acc_l,'MS_categ__rt_acc_long.csv',row.names = F)
#write.csv(df.C1.V.sum_rt_acc_noTask_l,'MS_categ__rt_acc_noTask_long.csv',row.names = F)
write.csv(df.v.sum_eff_all_w,'MS_cross_taskeffect_wide.csv',row.names = F)
setwd(curDir)


# ---------------------------------------------------------------------------------------
# ------------ 5. Plots  ----------------------------------------------------------------
# ---------------------------------------------------------------------------------------

# plot for matching task
df.M1.plot <- df.M1.V.sum_rt_acc_l %>%
   dplyr::filter(Match == 'match') %>%  # select matching data for plotting only.
   dplyr::full_join(., df.M1.V.SDT_l) %>%  
   tidyr::pivot_longer(., cols = c(RT, dprime), 
                       names_to = 'DVs', 
                       values_to = "value") %>% # to longer format
   dplyr::mutate(DVs = factor(DVs, levels = c('RT', 'dprime')),
                 # create an extra column for ploting the individual data cross different conditions.
                 Conds = mosaic::derivedFactor("0.88" = (Identity == "Self" & Morality == 'Good'), 
                                              "1.12" = (Identity == "Other" & Morality == 'Good'), 
                                              "1.88" = (Identity == "Self" & Morality == 'Bad'),
                                              "2.12" = (Identity == "Other" & Morality == 'Bad'), 
                                              method ="first", .default = NA),
                 Conds = as.numeric(as.character(Conds)))
   

df.M1.sum_p <- summarySE(df.M1.plot, measurevar = "value", groupvars = c("Identity", 'Morality',"DVs")) %>%
   dplyr::mutate(Mrl_num = ifelse(Morality == 'Good', 1, 2))

pd1 <- position_dodge(0.5)

# New facet label names for panel variable
# https://stackoverflow.com/questions/34040376/cannot-italicize-facet-labels-with-labeller-label-parsed
levels(df.M1.plot$DVs ) <- c("RT"=expression(paste("Reaction ", "times (ms)")),
                             "dprime"=expression(paste(italic("d"), ' prime')))
levels(df.M1.sum_p$DVs ) <- c("RT"=expression(paste("Reaction ", "times (ms)")),
                             "dprime"=expression(paste(italic("d"), ' prime')))

p_df_M1_sum <- df.M1.plot%>%
   ggplot(., aes(x = Morality, y = value, fill = Identity)) +
   geom_point(aes(x = Conds, y = value, group = Subject),   # plot individual points
              colour = "#000000",
              size = 3, shape = 20, alpha = 0.1)+
   geom_line(aes(x = Conds, y = value, group = Subject),         # link individual's points by transparent grey lines
             linetype = 1, size = 0.8, colour = "#000000", alpha = 0.06) +   
   geom_line(data = df.M1.sum_p, aes(x = as.numeric(Mrl_num), # plot the group means  
                                      y = value, 
                                      group = Identity, 
                                      colour = Identity), 
             linetype = 1, position = pd1, size = 2)+
   geom_point(data = df.M1.sum_p, aes(x = as.numeric(Mrl_num), # group mean
                                       y = value, 
                                       group = Identity, 
                                       colour = Identity), 
              shape = 18, position = pd1, size = 8) +
   geom_errorbar(data = df.M1.sum_p, aes(x = as.numeric(Mrl_num),  # group error bar.
                                          y = value, group = Identity, 
                                          colour = Identity,
                                          ymin = value- 1.96*se, 
                                          ymax = value+ 1.96*se), 
                 width = .05, position = pd1, size = 2, alpha = 0.75) +
   scale_colour_brewer(palette = "Dark2")+
   scale_x_continuous(breaks=c(1, 2),
                    labels=c("Good", "Bad"))+
   scale_fill_brewer(palette = "Dark2")+
   ggtitle("A. Matching task") +
   #theme_classic(base_size = 16) +
   theme_bw()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         panel.border = element_blank(),
         text=element_text(family='Times'),
         legend.title=element_blank(),
         legend.text = element_text(size =16),
         plot.title = element_text(lineheight=.8, face="bold", size = 18, margin=margin(0,0,20,0)),
         axis.text = element_text (size = 16, color = 'black'),
         axis.title = element_text (size = 16),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.line.x = element_line(color='black', size = 1),   # increase the size of font
         axis.line.y = element_line(color='black', size = 1),   # increase the size of font
         strip.text = element_text (size = 16, color = 'black'), # size of text in strips, face = "bold"
         panel.spacing = unit(3, "lines")
   ) +
   facet_wrap( ~ DVs,
               scales = "free_y", nrow = 1,
               labeller = label_parsed) 

pdf('Fig2_A_R1.pdf',  height = 6, width = 9)
p_df_M1_sum
dev.off()

### plot categorization task
df.C1.plot <- df.C1.V.sum_rt_acc_noTask_l %>%
   tidyr::pivot_longer(., cols = c(RT, ACC), 
                       names_to = 'DVs', 
                       values_to = "value") %>% # to longer format
   dplyr::mutate(DVs = factor(DVs, levels = c('RT', 'ACC')),
                 # create an extra column for ploting the individual data cross different conditions.
                 Conds = mosaic::derivedFactor("0.88" = (Identity == "Self" & Morality == 'Good'), 
                                               "1.12" = (Identity == "Other" & Morality == 'Good'), 
                                               "1.88" = (Identity == "Self" & Morality == 'Bad'),
                                               "2.12" = (Identity == "Other" & Morality == 'Bad'), 
                                               method ="first", .default = NA),
                 Conds = as.numeric(as.character(Conds)))


df.C1.sum_p <- summarySE(df.C1.plot, measurevar = "value", groupvars = c("Identity", 'Morality',"DVs")) %>%
   dplyr::mutate(Mrl_num = ifelse(Morality == 'Good', 1, 2))

levels(df.C1.plot$DVs ) <- c("RT"= expression("Reaction times (ms)"),
                             "ACC"= expression("Accuracy"))
levels(df.C1.sum_p$DVs ) <- c("RT"= expression("Reaction times (ms)"),
                              "ACC"= expression("Accuracy"))

p_df_C1_sum <- df.C1.plot%>%
   ggplot(., aes(x = Morality, y = value, fill = Identity)) +
   geom_point(aes(x = Conds, y = value, group = Subject),   # plot individual points
              colour = "#000000",
              size = 3, shape = 20, alpha = 0.1)+
   geom_line(aes(x = Conds, y = value, group = Subject),         # link individual's points by transparent grey lines
             linetype = 1, size = 0.8, colour = "#000000", alpha = 0.06) +   
   geom_line(data = df.C1.sum_p, aes(x = as.numeric(Mrl_num), # plot the group means  
                                     y = value, 
                                     group = Identity, 
                                     colour = Identity), 
             linetype = 1, position = pd1, size = 2)+
   geom_point(data = df.C1.sum_p, aes(x = as.numeric(Mrl_num), # group mean
                                      y = value, 
                                      group = Identity, 
                                      colour = Identity), 
              shape = 18, position = pd1, size = 8) +
   geom_errorbar(data = df.C1.sum_p, aes(x = as.numeric(Mrl_num),  # group error bar.
                                         y = value, group = Identity, 
                                         colour = Identity,
                                         ymin = value- 1.96*se, 
                                         ymax = value+ 1.96*se), 
                 width = .05, position = pd1, size = 2, alpha = 0.75) +
   scale_x_continuous(breaks=c(1, 2),
                      labels=c("Good", "Bad"))+
   scale_colour_brewer(palette = "Dark2")+
   scale_fill_brewer(palette = "Dark2")+
   ggtitle("B. Categorization task (combined)") +
   theme_bw()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         panel.border = element_blank(),
         text=element_text(family='Times'),
         legend.title=element_blank(),
         legend.text = element_text(size =16),
         plot.title = element_text(lineheight=.8, face="bold", size = 18, margin=margin(0,0,20,0)),
         axis.text = element_text (size = 16, color = 'black'),
         axis.title = element_text (size = 16),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.line.x = element_line(color='black', size = 1),   # increase the size of font
         axis.line.y = element_line(color='black', size = 1),   # increase the size of font
         strip.text = element_text (size = 16, color = 'black'), # size of text in strips, face = "bold"
         panel.spacing = unit(3, "lines")
   ) +
   facet_wrap( ~ DVs,
               scales = "free_y", nrow = 1) 

pdf('Fig2_B_R1.pdf',  height = 6, width = 9)
p_df_C1_sum
dev.off()