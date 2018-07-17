##
# This script transfer the raw data and save them to de-identified raw data
# create at 26 June. 2018
# by @Chuan-Peng Hu
# 

# clear the working memory
Sys.setlocale("LC_ALL", "English")  # set local encoding to English
Sys.setenv(LANG = "en") # set the feedback language to English
options(scipen = 999)   # force R to output in decimal instead of scientifc notion
options(digits=5)       # limit the number of reporting
rm(list = setdiff(ls(), lsf.str()))  # remove all data but keep functions


# set the directories
curDir <- getwd()   # directory for preprocessing
rootDir <- gsub('.{7}$', '', curDir)
rawDir <- paste(curDir,'/data/',sep = '')
traDir <- paste(rootDir,'tradAnal',sep = '')
ddmDir <- paste(rootDir,'hddmMod',sep = '')

# get the file names for files contain "data_exp7_rep_match_"
fNames.m <- list.files(path = rawDir, pattern = "data_exp7_rep_match_*")
# add the directory information to the filenames
fNames.m <- paste(rawDir,fNames.m, sep = '')

# read these files and combine them into one file, get the raw data for matching task
df.L <- do.call("rbind",lapply(as.list(fNames.m),FUN=function(files){read.table(files, header=TRUE, sep="",stringsAsFactors=F)}))

# do the same to data from categorization, get the raw data for categorization task
fNames.c <- list.files(path = rawDir, pattern = "data_exp7_rep_categ_*")
# remove invalid data data_exp7_rep_categ_7338.out
# fNames.c <- fNames.c[fNames.c != "data_exp7_rep_categ_7338.out"]
fNames.c <- paste(rawDir,fNames.c, sep = '')
df.T    <- do.call("rbind",lapply(as.list(fNames.c),FUN=function(files){read.table(files, header=TRUE, sep="",stringsAsFactors=F)}))


# render the data in numeric format for future analysis
cols.num <- c('Sub',"Age","Block","Bin","Trial","RT","ACC")
df.L[cols.num] <- sapply(df.L[cols.num], as.numeric)  
df.T[cols.num] <- sapply(df.T[cols.num], as.numeric)  
df.L <- df.L[!is.na(df.L$ACC),]
df.T <- df.T[!is.na(df.T$ACC),]

# get the independent varialbes
# moral valence
df.L$Morality[grepl("moral", df.L$Shape, fixed=TRUE)]   <- "Good"
df.L$Morality[grepl("immoral", df.L$Shape, fixed=TRUE)] <- "Bad"

# self-referential
df.L$Identity[grepl("Self", df.L$Shape, fixed=TRUE)]    <- "Self"
df.L$Identity[grepl("Other", df.L$Shape, fixed=TRUE)]   <- "Other"

df.T$Morality[grepl("moral", df.T$Shape, fixed=TRUE)]   <- "Good"
df.T$Morality[grepl("immoral", df.T$Shape, fixed=TRUE)] <- "Bad"

df.T$Identity[grepl("Self", df.T$Shape, fixed=TRUE)]    <- "Self"
df.T$Identity[grepl("Other", df.T$Shape, fixed=TRUE)]   <- "Other"


df.T$Task[df.T$Task == 'self']            <- 'Id'
df.T$Task[df.T$Task == 'moral']           <- 'Val'

# rename columns
colnames(df.L)[colnames(df.L)=="Sub"] <- "Subject"
colnames(df.T)[colnames(df.T)=="Sub"] <- "Subject"

# order the variables
df.L$Morality <- factor(df.L$Morality, levels = c("Good","Bad"))    
df.L$Identity <- factor(df.L$Identity, levels = c("Self","Other"))  
df.L$Match    <- factor(df.L$Match,    levels = c("match","mismatch"))

# make the variables in a specified order
df.T$Morality <- factor(df.T$Morality, levels=c("Good","Bad"))    
df.T$Identity <- factor(df.T$Identity, levels=c("Self","Other"))

# change all gender from number to string
df.L$Sex[df.L$Sex == '1'] <- 'male'
df.L$Sex[df.L$Sex == '2'] <- 'female'
df.T$Sex[df.T$Sex == '1'] <- 'male'
df.T$Sex[df.T$Sex == '2'] <- 'female'

df.L1 <- df.L[,c("Subject", "Age", "Sex", "Match", "Morality", "Identity", "ACC","RT")]
df.T1 <- df.T[,c("Subject", "Age", "Sex", "Task", "Morality", "Identity", "ACC","RT")]

write.csv(df.L1,'MS_rep_matchingTask_raw.csv',row.names = F)
write.csv(df.T1,'MS_rep_categTask_raw.csv',row.names = F)
