

curDir  <- dirname(rstudioapi::getSourceEditorContext()$path)   # get the directory for preprocessing
setwd(curDir)
source('Initial_exp7.r')  # initializing (clear global environment; load packages and functions)
curDir  <- dirname(rstudioapi::getSourceEditorContext()$path)   # get the directory for preprocessing

# ---------------------------------------------------------------------------------------
# ---------- 2. Loading data and clean the data        ----------------------------------
# ---------------------------------------------------------------------------------------

df.M.raw <- read.csv("MS_matchTask_raw.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE) # data for matching task
df.C.raw <- read.csv("MS_categTask_raw.csv",header = TRUE, sep = ',', stringsAsFactors=FALSE) # data for categorization task

# make the variables in a specified order
df.M.raw <- df.M.raw %>%
        dplyr::rename(ACC = Accuracy, Subject = SubjectID) %>%                                         # rename to ACC   
        dplyr::mutate(Morality = ifelse(Morality == "moral", 'Good','Bad')) %>%   # change condition name
        dplyr::mutate(Identity = ifelse(Identity == "self", 'Self','Other')) %>%  # change condition name
        dplyr::mutate(Morality = factor(Morality, levels = c("Good","Bad")),      # to factor
                      Identity = factor(Identity, levels=c("Self","Other")),
                      Match    = factor(Match,    levels=c("match","nonmatch")))

df.C.raw <- df.C.raw %>%
        dplyr::filter(Task == 'self' | Task == 'moral') %>%  # exclude trials, criterio 1: no importance-based trials
        dplyr::rename(ACC = Accuracy, Subject = SubjectID) %>%                                         # rename to ACC   
        dplyr::mutate(Morality = ifelse(Morality == "moral", 'Good','Bad')) %>%   # change condition name
        dplyr::mutate(Identity = ifelse(Identity == "self", 'Self','Other')) %>%  # change condition name
        dplyr::mutate(Task     = ifelse(Task     == "self", 'Id','Val')) %>%          # change condition name
        dplyr::mutate(Morality = factor(Morality, levels = c("Good","Bad")),      # to factor
                      Identity = factor(Identity, levels=c("Self","Other")))

# Exclude Subject, criterion 1: procedure failure
excldSub1 <- c("7","8","2027","7035")
df.M <- df.M.raw[!(df.M.raw$Subject %in% excldSub1),]
df.C <- df.C.raw[!(df.C.raw$Subject %in% excldSub1),]

# exclude trials, criterion 1: practicing trials in matching task (first 48 trials)
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
## NOTE: this was not done in previous analysis, resulting the wronly excluded participants in registration
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
df.M1 <- df.M.valid
df.C1 <- df.C.valid
df.M1$RT <- df.M1$RT * 1000 # transfer from seconds to min seconds
df.C1$RT <- df.C1$RT * 1000 # careful about the scale of time when using HDDM or ex-Gaussian

# no response equal to wrong
df.C1$ACC[df.C1$ACC == -1] <- 0

df.C1_a <- df.C1 %>% dplyr::filter(Bin < 4)  # first half
df.C1_b <- df.C1 %>% dplyr::filter(Bin > 3)  # second half

df.C1 <- df.C1 %>%
        dplyr::mutate(Finger = ifelse(ResponseKey == "Y" | ResponseKey == "H", "left", 
                                      ifelse(ResponseKey == "U" | ResponseKey == "J", "right", NA))) %>%
        dplyr::mutate(Finger = ifelse(ACC == 1 | (ACC = 0 & ResponseKey = NA), Finger, 
                                      ifelse(ACC = 0 & ResponseKey = 'left', 'right',
                                             ifelse(ACC = 0 & ResponseKey ='right', 'left'))))

# calculating the ACC 
df.C1_a.acc <- df.C1_a %>%
        dplyr::group_by(Subject,Age, Gender, Task,Morality,Identity) %>%
        dplyr::summarise(N = length(ACC),
                         corrN = sum(ACC),
                         meanACC = sum(ACC)/length(ACC)) %>%                                    # calculate the counts for each 
        dplyr::ungroup() %>%
        dplyr::mutate(Morality = factor(Morality, levels = c("Good","Bad")),
                      Identity = factor(Identity, levels = c("Self","Other")))

afex::aov_ez("Subject", "meanACC", df.C1_a.acc, within = c('Task',"Morality", "Identity"))

df.C1_b.acc <- df.C1_b %>%
        dplyr::group_by(Subject,Age, Gender, Task,Morality,Identity) %>%
        dplyr::summarise(N = length(ACC),
                         corrN = sum(ACC),
                         meanACC = sum(ACC)/length(ACC)) %>%                                    # calculate the counts for each 
        dplyr::ungroup() %>%
        dplyr::mutate(Morality = factor(Morality, levels = c("Good","Bad")),
                      Identity = factor(Identity, levels = c("Self","Other")))

afex::aov_ez("Subject", "meanACC", df.C1_b.acc, within = c('Task',"Morality", "Identity"))

df.C1_a.RT <- df.C1_a %>%
        dplyr::filter(ACC == 1)  # exclued inaccurate data

df.C1_a.RT.subj <- summarySEwithin(df.C1_a.RT, measurevar = 'RT', withinvar = c('Subject','Task','Morality','Identity'), idvar = 'Subject',na.rm = TRUE)

afex::aov_ez("Subject", "RT", df.C1_a.RT.subj, within = c('Task',"Morality", "Identity"))

df.C1_b.RT <- df.C1_b %>%
        dplyr::filter(ACC == 1)  # exclued inaccurate data

df.C1_b.RT.subj <- summarySEwithin(df.C1_b.RT, measurevar = 'RT', withinvar = c('Subject','Task','Morality','Identity'), idvar = 'Subject',na.rm = TRUE)
afex::aov_ez("Subject", "RT", df.C1_b.RT.subj, within = c('Task',"Morality", "Identity"))

df.C1_a.V.sum_rt_acc_l <- merge(df.C1_a.acc, df.C1_a.RT.subj,by = c("Subject","Task","Morality",'Identity'))
df.C1_a.V.sum_rt_acc_l <- df.C1_a.V.sum_rt_acc_l[order(df.C1_a.V.sum_rt_acc_l$Subject),]
df.C1_a.V.sum_rt_acc_l <- df.C1_a.V.sum_rt_acc_l %>%
        dplyr::rename(ACC = meanACC)

df.C1_b.V.sum_rt_acc_l <- merge(df.C1_b.acc, df.C1_b.RT.subj,by = c("Subject","Task","Morality",'Identity'))
df.C1_b.V.sum_rt_acc_l <- df.C1_b.V.sum_rt_acc_l[order(df.C1_b.V.sum_rt_acc_l$Subject),]
df.C1_b.V.sum_rt_acc_l <- df.C1_b.V.sum_rt_acc_l %>%
        dplyr::rename(ACC = meanACC)

### plot id-based data
CAplots(saveDir = curDir, curDir = curDir,expName = 'conf_b', task = 'id', df.C1_b.V.sum_rt_acc_l)

### plot val-based data
CAplots(saveDir = curDir, curDir = curDir,expName = 'conf_b', task = 'val', df.C1_b.V.sum_rt_acc_l)

### plot the categorization task (collapsed different tasks)
#CAplots(saveDir = curDir, curDir = curDir,expName = 'conf_a', task = 'categ', df.C1_a.V.sum_rt_acc_l)
