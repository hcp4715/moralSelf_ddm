# this script is used for initializing the analysis
# preparing necessary functions used in current analysis

Sys.setlocale("LC_ALL", "English")  # set local encoding to English
Sys.setenv(LANG = "en") # set the feedback language to English
options(scipen = 999)   # force R to output in decimal instead of scientifc notion
options(digits=5)       # limit the number of reporting
#rm(list = setdiff(ls(), lsf.str()))  # remove all data but keep functions
rm(list = ls())
pkgTest <- function(x)
{
        if (!require(x,character.only = TRUE))
        {
                install.packages(x,dep = TRUE)
                if(!require(x,character.only = TRUE)) stop("Package not found")
        }
}

pkgNeeded <- (c("plyr","tidyverse","ggplot2","ez", "bootES","MBESS", 
                "psych","corrplot","readr", "Hmisc","RColorBrewer"))

lapply(pkgNeeded,pkgTest)
rm('pkgNeeded') # remove the variable 'pkgNeeded';

# run the geo_flat_violin.r, which is from:https://gist.githubusercontent.com/
# benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R
source("geom_flat_violin.R")

# Save some time and stor APA format-related code in an object so you can easily
# use it in multiple plots
windowsFonts(Times=windowsFont("TT Times New Roman")) # explicit mapping to "times"
apatheme = theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              text=element_text(family='Times'),
              legend.title=element_blank(),
              legend.text = element_text(size =16),
              #legend.position='top',
              plot.title = element_text(lineheight=.8, face="bold", size = 16),
              axis.text = element_text (size = 16, color = 'black'),
#              axis.text.x = element_text(angle = 45, vjust = 0.5),   # x-axis's label font
              axis.title = element_text (size = 16),
              axis.title.x = element_text(margin=margin(10,0,0,0)),  # increase the sapce betwen title and x axis
              axis.title.y = element_text(margin=margin(0,12,0,0)),  # increase the space between title and y axis
              axis.line.x = element_line(color='black', size = 1),   # increase the size of font
              axis.line.y = element_line(color='black', size = 1))   # increase the size of font

# define the d prime function
dprime <- function(hit,fa) {
        qnorm(hit) - qnorm(fa)
}

## below is the code from blog, and adapted from A C Del Re from email
d.sgpp <- function(m.1,m.2,sd.1,sd.2,n,r=.5)
{
        # m.1 = mean at time 1
        # m.2 = mean at time 2
        # sd.1 = standard dev. at time 1
        # sd.2 = standard dev. at time 2
        # n = sample size
        # r = correlation between time 1 and 2
        s.within <- (sqrt((sd.1^2 + sd.2^2)-2*r*sd.1*sd.2))/(sqrt(2*(1-r))) 
        d <- ((m.1 - m.2)/s.within)
        var.d <- 2*(1-r) * (1/n + d^2/(2*n))
        out <- cbind(d, var.d)
        return(out)
}


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
        
        datac <- plyr::rename(datac,c("mean" = measurevar))
        
        datac$se <- datac$sd /sqrt(datac$N)   # calculate standard error of the mean
        
        # Confidence interval mltiplier for standard error
        # calculate t-statistic for confidence interval:
        # e.g., if conf.interval is .95, use .975 (above/below), and use df = N-1
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
#


########### define a function for the plots ##########
#### For categorization task
CAplots <- function(saveDir = traDir, curDir = curDir,expName = 'exp7', task = 'id', inData){
      inData$Identity <- factor(inData$Identity,levels = c("Self","Other"))
      inData$Morality <- factor(inData$Morality,levels = c("Good","Bad"))
      #inData$Morality[inData$Morality == "Good"] <- 1
      #inData$Morality[inData$Morality == "Bad"]  <- 2
      if(task == 'val'){              # valence-based categorization
            ACCdata <- inData %>%
                  select(Subject,Task,Morality,Identity,ACC) %>% 
                  filter(Task == "Val")
            rtData <- inData %>%
                  select(Subject,Task,Morality,Identity,RT) %>% 
                  filter(Task == "Val")
            
          } else if (task == 'id'){   # id-based categorization
            ACCdata <- inData %>%
                  select(Subject,Task,Morality,Identity,ACC) %>% 
                  filter(Task == "Id")
            rtData <- inData %>%
                  select(Subject,Task,Morality,Identity,RT) %>% 
                  filter(Task == "Id")
          }else{                         #  combined for experiment 1
            ACCdata <- inData %>%
                  select(Subject,Morality,Identity,ACC)
            rtData <- inData %>%
                  select(Subject,Morality,Identity,RT)
            
      }

    P.acc <- ggplot(ACCdata,aes(x = Morality, 
                                y = ACC, fill = Identity))+
          geom_flat_violin(aes(fill = Identity),position = position_nudge(x = 0.1, y = 0),
                           adjust = 1.5, trim = FALSE, alpha = 0.5,color = NA) +
          geom_dotplot(aes(x = Morality,y = ACC, color = Identity), 
                       binaxis='y', binwidth = 0.0125, stackdir='center', dotsize= 0.5,position = position_dodge(0.15)) +
          geom_boxplot(aes(x = Morality,  y = ACC,fill = Identity),outlier.shape = NA,
                       alpha = 0.5, width = 0.1,  color = "black",
                       position = position_dodge(0.15))+ 
          scale_color_brewer(palette = "Dark2")+
          scale_fill_brewer(palette = "Dark2")+
          ylab("Accuracy")+
          #scale_x_discrete(breaks = c(1,2),labels = c("Good","Bad")) +
          scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
          apatheme
    
    fileName = paste0('p_',expName,'_',task,'_ACC','.pdf')
    ggsave(fileName, P.acc, scale = 1,height = 6, width = 6, dpi = 300, family = "Times",path = saveDir)
    
    
    P.rt <- ggplot(rtData,aes(x = Morality, y = RT, fill = Identity))+
          geom_flat_violin(aes(fill = Identity),position = position_nudge(x = 0.1, y = 0),
                           adjust = 1.5, trim = FALSE, alpha = 0.5,color = NA) +
          #geom_point(aes(x = as.numeric(Morality)-0.15,y = RT, color = Identity), 
          #           position = position_jitter(width = 0.02),size = 1, shape = 20)+
          geom_dotplot(aes(x = Morality,y = RT, color = Identity), 
                       binaxis='y', binwidth = 8, stackdir='center', dotsize= 0.5,position = position_dodge(0.15)) + 
          geom_boxplot(aes(x = Morality,  y = RT,fill = Identity),outlier.shape = NA,
                       alpha = 0.5, width = 0.1,  color = "black",
                       position = position_dodge(0.15))+ 
          scale_color_brewer(palette = "Dark2")+
          scale_fill_brewer(palette = "Dark2")+
          ylab("Reaction Times")+
          #scale_x_discrete(breaks = c(1,2),labels = c("Good","Bad")) +
          scale_y_continuous(expand = c(0, 0),limits = c(200,1000))+
          apatheme
    
    fileName = paste0('p_',expName,'_',task,'_RT','.pdf')
    ggsave(fileName, P.rt, scale = 1,height = 6, width = 6, dpi = 300, family = "Times",path = saveDir)
    
    fileName = paste0('p_',expName,'_',task,'.tiff')
    setwd(saveDir)
    tiff(fileName, width = 9, height = 6, units = 'in', res = 300)
    p_dprime_match <- multiplot(P.rt,P.acc,cols = 2)
    dev.off()
    setwd(curDir)
    return(multiplot(P.rt,P.acc,cols = 2))
}
 
#### For Match task
Mplots <- function(saveDir = traDir, curDir = curDir, expName = 'exp7', dData,rtData){
      dData$Identity <- factor(dData$Identity,levels = c("Self","Other"))
      dData$Morality <- factor(dData$Morality,levels = c("Good","Bad"))
      rtData$Identity <- factor(rtData$Identity,levels = c("Self","Other"))
      rtData$Morality <- factor(rtData$Morality,levels = c("Good","Bad"))

      P.dprime <- ggplot(dData,aes(x = Morality, y = dprime, fill = Identity)) +
            geom_flat_violin(aes(fill = Identity),position = position_nudge(x = 0.1, y = 0),
                             adjust = 1.5, trim = FALSE, alpha = 0.5,color = NA) +
            geom_dotplot(aes(x = Morality,y = dprime, color = Identity), 
                         binaxis='y', binwidth = 0.1, stackdir='center', dotsize= 0.5,
                         position = position_dodge(0.15)) +
            geom_boxplot(aes(x = Morality,  y = dprime,fill = Identity),outlier.shape = NA,
                         alpha = 0.5, width = 0.1,  color = "black",
                         position = position_dodge(0.15)) + 
            scale_color_brewer(palette = "Dark2") +
            scale_fill_brewer(palette = "Dark2") +
            ylab("d prime") +
           # scale_x_discrete(breaks = c(1,2),labels = c("Good","Bad")) +
            scale_y_continuous(expand = c(0, 0), limits = c(-1,5)) +
            apatheme
      fileName = paste0('p_',expName,'_match_dprime','.pdf')
      ggsave(fileName, P.dprime, scale = 1,height = 6, width = 6, dpi = 300, family = "Times",path = saveDir)
      
      P.rt <- ggplot(rtData,aes(x = Morality, y = RT, fill = Identity))+
            geom_flat_violin(aes(fill = Identity),position = position_nudge(x = 0.1, y = 0),
                             adjust = 1.5, trim = FALSE, alpha = 0.5,color = NA) +
            #geom_point(aes(x = as.numeric(Morality)-0.15,y = RT, color = Identity), 
            #           position = position_jitter(width = 0.02),size = 1, shape = 20)+
            geom_dotplot(aes(x = Morality,y = RT, color = Identity), 
                         binaxis='y', binwidth = 8, stackdir='center', dotsize= 0.5,position = position_dodge(0.15)) + 
            geom_boxplot(aes(x = Morality,  y = RT,fill = Identity),outlier.shape = NA,
                         alpha = 0.5, width = 0.1,  color = "black",
                         position = position_dodge(0.15))+ 
            scale_color_brewer(palette = "Dark2")+
            scale_fill_brewer(palette = "Dark2")+
            ylab("Reaction Times")+
            #scale_x_discrete(breaks = c(1,2),labels = c("Good","Bad")) +
            scale_y_continuous(expand = c(0, 0),limits = c(200,1000))+
            apatheme
      fileName = paste0('p_',expName,'_match_RT','.pdf')
      ggsave(fileName, P.rt, scale = 1,height = 6, width = 6, dpi = 300, family = "Times",path = saveDir)
  
      fileName = paste0('p_',expName,'_match_','.tiff')
      setwd(saveDir)
      tiff(fileName, width = 12, height = 6, units = 'in', res = 300)
      p_dprime_match <- multiplot(P.rt,P.dprime,cols = 2)
      dev.off()
      setwd(curDir)
      return(multiplot(P.rt,P.dprime,cols = 2))
}
