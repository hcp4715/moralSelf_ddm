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
# ---------- 5. Questionnaire: prepare and save  ----------------------------------------
# ---------- 6. Plots -------------------------------------------------------


# ---------------------------------------------------------------------------------------
# ---------- 1. Initializing and prepare               ----------------------------------
# ---------------------------------------------------------------------------------------

source('Initial.r')  # initializing (clear global environment; load packages and functions)
curDir  <- dirname(rstudioapi::getSourceEditorContext()$path)   # get the directory for preprocessing
rootDir <- gsub('.{7}$', '', curDir)                      # get the parental folder
traDir <- paste(rootDir,'traditional_analysis',sep = '')  # folder for traditional analsysi
ddmDir <- paste(rootDir,'hddm',sep = '')                  # folder for DDM analysis
exgDir <- paste(rootDir,'exGaussian',sep = '')            # folder for ExGaussian analysis

# ---------------------------------------------------------------------------------------
# ---------- 2. Loading data and clean the data        ----------------------------------
# ---------------------------------------------------------------------------------------

df.M.raw <- read.csv("MS_matchTask_raw.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE) # data for matching task
df.C.raw <- read.csv("MS_categTask_raw.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE) # data for categorization task

subTotl <- unique(df.M.raw$Subject)

# make the variables in a specified order
df.M.raw$Morality <- factor(df.M.raw$Morality, levels=c("Good","Bad"))    
df.M.raw$Identity <- factor(df.M.raw$Identity, levels=c("Self","Other"))  
df.M.raw$Match    <- factor(df.M.raw$Match,    levels=c("match","nonmatch"))
df.C.raw$Morality <- factor(df.C.raw$Morality, levels=c("Good","Bad"))    
df.C.raw$Identity <- factor(df.C.raw$Identity, levels=c("Self","Other"))

# Exclude Subject, criterion 1: procedure failure
excldSub1 <- c("7","8","2027","7035")
df.M <- df.M.raw[!(df.M.raw$Subject %in% excldSub1),]
df.C <- df.C.raw[!(df.C.raw$Subject %in% excldSub1),]

# exclude trials, criterio 1: no importance-based trials
df.C <- df.C[df.C$Task != 'importance',]  # exclude the results from importance task, should be 17856 rows

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
df.M.tmp <- df.M
df.M.tmp$ACC[df.M.tmp$ACC == -1] <- 0
df.C.tmp <- df.C
df.C.tmp$ACC[df.C.tmp$ACC == -1] <- 0 # this was not done in previous analysis, resulting the wronly excluded participants in registration

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
df.M1 <- df.M.valid
df.C1 <- df.C.valid
df.M1$RT <- df.M1$RT * 1000 # transfer from seconds to min seconds
df.C1$RT <- df.C1$RT * 1000 # careful about the scale of time when using HDDM or ex-Gaussian

# no response equal to wrong
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

df.M.hddm_m <- df.M.valid %>%                          # start from the valid data, RT in seconds.
      dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
      dplyr::filter(Match == 'match') %>%
      dplyr::select(Subject,Morality,Identity,RT,ACC) %>%  # select column
      dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, response = ACC) # rename column

df.M.hddm_nm <- df.M.valid %>%
      dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
      dplyr::filter(Match == 'nonmatch') %>%
      dplyr::select(Subject,Morality,Identity,RT,ACC) %>%  # select column
      dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, response = ACC) # rename column

df.C.hddm_val <- df.C.valid %>%
      dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
      dplyr::filter(Task == 'Val') %>%                     # select the valence-based task
      dplyr::select(Subject,Morality,Identity,RT,ACC) %>%  # select column
      dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, response = ACC) # rename column

df.C.hddm_id <- df.C.valid %>%
      dplyr::filter(ACC == 1 | ACC == 0) %>%               # exclude trials without response or with wrong keys
      dplyr::filter(Task == 'Id') %>%                      # select the valence-based task
      dplyr::select(Subject,Morality,Identity,RT,ACC) %>%  # select column
      dplyr::rename(subj_idx = Subject, val = Morality, id = Identity, rt = RT, response = ACC) # rename column

setwd(ddmDir)
write.csv(df.M.hddm_m,'MS_match_hddm.csv',row.names = F)
write.csv(df.M.hddm_nm,'MS_mismatch_hddm.csv',row.names = F)
write.csv(df.C.hddm_val,'MS_categ_val_hddm.csv',row.names = F)
write.csv(df.C.hddm_id,'MS_categ_id_hddm.csv',row.names = F)
setwd(curDir)
rm('df.M.hddm_m','df.M.hddm_nm','df.C.hddm_val','df.C.hddm_id')

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
                       ACC = sum(ACC)/length(ACC))

# transform to wide-format
df.M1.V.acc_w <- dcast(df.M1.V.acc, Subject ~ Match + Morality + Identity,value.var = "ACC")

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
df.M1.V.SDT <- plyr::ddply(df.M1.V,.(Subject,Age, Gender, Morality, Identity,sdt), summarise, N = length(sdt))
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
df.M1.V.SDT_w$dprime <- mapply(dprime,df.M1.V.SDT_w$hitR,df.M1.V.SDT_w$faR)

# transfor from long format to wide format
df.M1.V.SDT_ww <- dcast(df.M1.V.SDT_w, Subject + Age + Gender ~ Morality + Identity ,value.var = "dprime")
# rename the column number
colnames(df.M1.V.SDT_ww)[4:7] <- paste("d", colnames(df.M1.V.SDT_ww[,4:7]), sep = "_")

df.M1.V.SDT_l <- df.M1.V.SDT_w[,c(1:5,12)]

####################################
############      RT     ###########
####################################

df.M1.V.RT <- df.M1.V[df.M1.V$ACC == 1,]

# get the summary results
df.M1.V.RT.subj <- summarySEwithin(df.M1.V.RT,measurevar = 'RT', withinvar = c('Subject','Match','Morality','Identity'),idvar = 'Subject', na.rm = TRUE)
df.M1.V.RT.subj_w <- dcast(df.M1.V.RT.subj, Subject ~ Match + Morality + Identity ,value.var = "RT") 

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

df.C1.V.acc <- plyr::ddply(df.C1.V,.(Subject,Age, Gender, Task,Morality,Identity), summarise,
                           N = length(ACC),
                           countN = sum(ACC),
                           ACC = sum(ACC)/length(ACC))

# wide-format
df.C1.V.acc_w  <- dcast(df.C1.V.acc, Subject ~ Task + Morality + Identity ,value.var = "ACC")
# rename the column number
colnames(df.C1.V.acc_w)[2:9] <- paste("ACC", colnames(df.C1.V.acc_w[,2:9]), sep = "_")

# combing data from diff task for analyzing the interaction btw val and id
df.C1.V.acc_noTask  <-  plyr::ddply(df.C1.V,.(Subject, Morality, Identity), summarise,
                                    N = length(ACC),
                                    countN = sum(ACC),
                                    ACC = sum(ACC)/length(ACC))

df.C1.V.acc_noTask_w <- dcast(df.C1.V.acc_noTask, Subject ~ Morality + Identity,value.var = "ACC")
# rename the column number
colnames(df.C1.V.acc_noTask_w)[2:5] <- paste("ACC", colnames(df.C1.V.acc_noTask_w[,2:5]), sep = "_")

# calculate the mean differences(?)
df.C1.V.acc_noTask_w$good_bad_slf <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Bad_Self
df.C1.V.acc_noTask_w$good_bad_oth <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Bad_Other
df.C1.V.acc_noTask_w$slf_oth_good <- df.C1.V.acc_noTask_w$ACC_Good_Self - df.C1.V.acc_noTask_w$ACC_Good_Other
df.C1.V.acc_noTask_w$slf_oth_bad <- df.C1.V.acc_noTask_w$ACC_Bad_Self - df.C1.V.acc_noTask_w$ACC_Bad_Other

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
write.csv(df.C1.V.sum_rt_acc_l,'MS_categ__rt_acc_long.csv',row.names = F)
write.csv(df.C1.V.sum_rt_acc_noTask_l,'MS_categ__rt_acc_noTask_long.csv',row.names = F)
write.csv(df.v.sum_eff_all_w,'MS_cross_taskeffect_wide.csv',row.names = F)
setwd(curDir)

# ---------------------------------------------------------------------------------------
# ------------ 5. Questionnaire: prepare and save  --------------------------------------
# ---------------------------------------------------------------------------------------

# load the data
df.quest <- read.csv("exp7_survey_monkey_raw_c.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE)

# calculate the reliability of each questionnair
# get the name and key for each scale
slfEsteemNames <- c( "slfEsteem_1", "slfEsteem_2", "slfEsteem_3", "slfEsteem_4", "slfEsteem_5",
                     "slfEsteem_6", "slfEsteem_7", "slfEsteem_8", "slfEsteem_9", "slfEsteem_10")
slfEsteemKeys  <- c(1,2,-3,4,-5,6,7,-8,-9,-10)        # for current dataset
slfEsteemKeys2 <- c( "slfEsteem_1", "slfEsteem_2", "-slfEsteem_3", "slfEsteem_4", "-slfEsteem_5",
                     "slfEsteem_6", "slfEsteem_7", "-slfEsteem_8", "-slfEsteem_9", "-slfEsteem_10")
slfEsteemAlpha <- psych::alpha(df.quest[,slfEsteemNames],keys = slfEsteemKeys) # alpha
print(slfEsteemAlpha$total[2])
slfEsteemOmega <- psych::omega(df.quest[,slfEsteemNames]) # omega
print(c(slfEsteemOmega$omega_h,slfEsteemOmega$omega.tot))

# calculate the score for the trait self-esteem
slfEsteemScore <- psych::scoreItems(slfEsteemKeys2,df.quest[,slfEsteemNames], totals = F, min = 1, max = 4)
df.quest$slfEsteem <- slfEsteemScore$scores

# calculate for moral identity
moralIdNames <- c("MoralId_1","MoralId_2","MoralId_3","MoralId_4","MoralId_5","MoralId_6",
                  "MoralId_7","MoralId_8","MoralId_9","MoralId_10","MoralId_11","MoralId_12",
                  "MoralId_13","MoralId_14","MoralId_15")
moralIdKeys <- c(1,2,3,4,-5,6,7,8,9,10,11,12,13,14,15)

moralIdAlpha <- psych::alpha(df.quest[,moralIdNames],keys = moralIdKeys) # alpha
print(moralIdAlpha$total[2])
moralIdOmega <- psych::omega(df.quest[,moralIdNames]) # omega
print(c(moralIdOmega$omega_h,moralIdOmega$omega.tot))

moralIdinNames <- c("MoralId_1","MoralId_2","MoralId_5","MoralId_8","MoralId_10",
                    "MoralId_11","MoralId_12", "MoralId_13","MoralId_14")
moralIdinKeys <- c(1,2,-3,4,5,6,7,8,9)
moralIdinKeys2 <- c("MoralId_1","MoralId_2","-MoralId_5","MoralId_8","MoralId_10",
                    "MoralId_11","MoralId_12", "MoralId_13","MoralId_14")
moralIdinAlpha <- psych::alpha(df.quest[,moralIdinNames],keys = moralIdinKeys) # alpha
print(moralIdinAlpha$total[2])
moralIdinOmega <- psych::omega(df.quest[,moralIdinNames]) # omega
print(c(moralIdinOmega$omega_h,moralIdinOmega$omega.tot))

# caculate the score for internal moral identity
moralIdinScore <- psych::scoreItems(moralIdinKeys2,df.quest[,moralIdinNames], totals = F, min = -2, max = 2)
df.quest$moralIdinScore <- moralIdinScore$scores

moralIdexNames <- c("MoralId_3","MoralId_4","MoralId_6",
                    "MoralId_7","MoralId_9","MoralId_15")
moralIdexKeys <- c(1,2,3,4,5,6)
moralIdexAlpha <- psych::alpha(df.quest[,moralIdexNames],keys = moralIdexKeys) # alpha
print(moralIdexAlpha$total[2])
moralIdexOmega <- psych::omega(df.quest[,moralIdexNames]) # omega
print(c(moralIdexOmega$omega_h,moralIdexOmega$omega.tot))

# calculate the score for external moral identity
moralIdexScore <- psych::scoreItems(moralIdexNames,df.quest[,moralIdexNames], totals = F, min = -2, max = 2)
df.quest$moralIdexScore <- moralIdexScore$scores

# find the intersection between behavioral data and questionnarie data.
tmp <- merge(df.v.sum_eff_all_w,df.quest[,c("subjID","slfEsteem","moralIdinScore","moralIdexScore")],by.x = "Subject",by.y = "subjID",all = TRUE)
tmp <- tmp[complete.cases(tmp),]

# write files to an upper-level folder
setwd(traDir)
write.csv(tmp,'MS_cross_taskeffect_wide_w_scale.csv',row.names = F)
setwd(curDir)

# ---------------------------------------------------------------------------------------
# ------------ 6. Plots  ----------------------------------------------------------------
# ---------------------------------------------------------------------------------------

## Matching task, plots are save to 'saveDir'
Mplots(saveDir = traDir, curDir = curDir, expName = 'exp7', df.M1.V.SDT_l,df.C1.V.sum_rt_acc_l)

# plot id-based data
CAplots(saveDir = traDir, curDir = curDir,expName = 'exp7', task = 'id', df.C1.V.sum_rt_acc_l)

# plot val-based data
CAplots(saveDir = traDir, curDir = curDir,expName = 'exp7', task = 'val', df.C1.V.sum_rt_acc_l)

# plot the categorization task (collapsed different tasks)
CAplots(saveDir = traDir, curDir = curDir,expName = 'exp7', task = 'categ', df.C1.V.sum_rt_acc_noTask_l)

