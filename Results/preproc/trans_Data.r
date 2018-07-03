# this script is used for transfer raw data to usable and sharable data

# initializing
source('Initial.r')

# set directory
curDir  <- dirname(rstudioapi::getSourceEditorContext()$path)   # directory for preprocessing
rootDir <- gsub('.{7}$', '', curDir)
rawDir  <- paste(curDir,'/rawData/',sep = '')
traDir  <- paste(rootDir,'Traditional_Analysis',sep = '')
ddmDir  <- paste(rootDir,'HDDM',sep = '')

setwd(rawDir)
file_list.match <- list.files(path = ".", pattern = "*_learn*")
file_list.categ <- list.files(path = ".", pattern = "*_test*")

# merge the data from each participant
for (fileL in file_list.match){
  # if the merged dataset doesn't exist, creat it
  if (!exists("df.M")){
    df.M <- read.csv(fileL,header = TRUE, sep = ' ', stringsAsFactors=FALSE)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("df.M")){
    temp_df.L <- read.csv(fileL,header = TRUE, sep = ' ', stringsAsFactors=FALSE)
    df.M <- rbind(df.M,temp_df.L)
    rm(temp_df.L)
  }
}
df.M <- subset(df.M,SubjectID != "SubjectID")  # there should be 18672 obs. of 19 variables.
# write.csv(df.M,"df.Learn.csv", row.names=FALSE)

for (fileT in file_list.categ){
  # if the merged dataset doesn't exist, creat it
  if (!exists("df.C")){
    df.C <- read.csv(fileT,header = TRUE, sep = ' ', stringsAsFactors=FALSE)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("df.C")){
    temp_df.T <- read.csv(fileT,header = TRUE, sep = ' ', stringsAsFactors=FALSE)
    df.C <- rbind(df.C,temp_df.T)
    rm(temp_df.T)
  }
}
df.C <- subset(df.C,SubjectID != "SubjectID")   # there should be 29808 obs. of 19 variables
# write.csv(df.C,"df.Test.csv", row.names=FALSE)  # save the data

setwd(curDir)
# render the data in numeric format for future analysis
cols.num <- c("Age","Block","Bin","Trial","RT","Accuracy")
df.M[cols.num] <- sapply(df.M[cols.num], as.numeric)  
df.C[cols.num] <- sapply(df.C[cols.num], as.numeric)  

# rename columns
colnames(df.M)[colnames(df.M)=="Accuracy"]  <- "ACC"
colnames(df.C)[colnames(df.C)=="Accuracy"]  <- "ACC"
colnames(df.M)[colnames(df.M)=="SubjectID"] <- "Subject"
colnames(df.C)[colnames(df.C)=="SubjectID"] <- "Subject"

df.M$Subject <- factor(df.M$Subject)
df.C$Subject <- factor(df.C$Subject)

# Change the name for different levels of independent variables
df.M$Morality[df.M$Morality == 'moral']   <- 'Good'
df.M$Morality[df.M$Morality == 'immoral'] <- 'Bad'
df.M$Identity[df.M$Identity == 'self']    <- 'Self'
df.M$Identity[df.M$Identity == 'other']   <- 'Other'

# make the variables in a specified order
df.M$Morality <- factor(df.M$Morality, levels = c("Good","Bad"))    
df.M$Identity <- factor(df.M$Identity, levels = c("Self","Other"))  
df.M$Match    <- factor(df.M$Match,    levels = c("match","nonmatch"))

# Change the name for different levels of independent variables
df.C$Morality[df.C$Morality == 'moral']   <- 'Good'
df.C$Morality[df.C$Morality == 'immoral'] <- 'Bad'
df.C$Identity[df.C$Identity == 'self']    <- 'Self'
df.C$Identity[df.C$Identity == 'other']   <- 'Other'
df.C$Task[df.C$Task == 'Self']            <- 'Id'
df.C$Task[df.C$Task == 'self']            <- 'Id'
df.C$Task[df.C$Task == 'moral']           <- 'Val'

# make the variables in a specified order
df.C$Morality <- factor(df.C$Morality, levels=c("Good","Bad"))    
df.C$Identity <- factor(df.C$Identity, levels=c("Self","Other"))

df.M1 <- df.M[,c("Subject", "Age", "Gender", "Match", "Morality", "Identity", "ACC","RT")]
df.C1 <- df.C[,c("Subject", "Age", "Gender", "Task", "Morality", "Identity", "ACC","RT")]

write.csv(df.M1,'MS_matchingTask_raw.csv',row.names = F)
write.csv(df.C1,'MS_categTask_raw.csv',row.names = F)
