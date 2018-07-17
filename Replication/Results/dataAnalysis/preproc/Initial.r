# this script is used for initializing the analysis
# preparing necessary functions used in current analysis

Sys.setlocale("LC_ALL", "English")  # set local encoding to English
Sys.setenv(LANG = "en") # set the feedback language to English
options(scipen = 999)   # force R to output in decimal instead of scientifc notion
options(digits=5)       # limit the number of reporting
rm(list = setdiff(ls(), lsf.str()))  # remove all data but keep functions

pkgTest <- function(x)
{
        if (!require(x,character.only = TRUE))
        {
                install.packages(x,dep = TRUE)
                if(!require(x,character.only = TRUE)) stop("Package not found")
        }
}

pkgNeeded <- (c("tidyverse","ggplot2", "reshape2","ez", "bootES","MBESS", 
                "BayesFactor","psych","corrplot",'plyr','ggstatsplot',"readr",
                "tidyr","Hmisc","RColorBrewer"))

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
              legend.position='top',
              plot.title = element_text(size = 16),
              axis.text = element_text (size = 16, color = 'black'),
              axis.title = element_text (size = 16),
              axis.title.y = element_text(margin=margin(0,16,0,0)),  # increase the space between title and y axis
              axis.title.x = element_text(margin=margin(16,0,0,0)),  # increase the sapce betwen title and x axis
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
raincloud_theme <-  theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16,colour = 'black'),
  legend.text=element_text(size=16,colour = 'black'),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16,colour = 'black'),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

# define a function for the plots
MSplots <- function(saveDir = traDir, curDir = curDir, task = 'match',type = 'ACC', inData){
  if(type == 'ACC'){
    if(task == 'val'){
      ACCdata1 <- inData %>%
        select(Subject,Task,Morality,Identity,ACC) %>% 
        filter(Task == "Val"& Morality == 'Good')
      ACCdata2 <- inData %>%
        select(Subject,Task,Morality,Identity,ACC) %>% 
        filter(Task == "Val"& Morality == 'Bad')
    } else if (task == 'id'){
      ACCdata1 <- inData %>%
        select(Subject,Task,Morality,Identity,ACC) %>% 
        filter(Task == "Id"& Morality == 'Good')
      ACCdata2 <- inData %>%
        select(Subject,Task,Morality,Identity,ACC) %>% 
        filter(Task == "Id"& Morality == 'Bad')
    }else{
      ACCdata1 <- inData %>%
        select(Subject,Morality,Identity,ACC) %>% 
        filter(Morality == 'Good')
      ACCdata2 <- inData %>%
        select(Subject,Morality,Identity,ACC) %>% 
        filter(Morality == 'Bad')
    }

    ACCdata1$Identity <- factor(ACCdata1$Identity,levels = c("Self","Other"))
    ACCdata2$Identity <- factor(ACCdata2$Identity,levels = c("Self","Other"))
    p1 <- ggplot(data = ACCdata1, aes(y = ACC, x = Identity,fill = Identity)) +
      geom_flat_violin(position = position_nudge(x = .1, y = 0)) +
      geom_point(aes(y = ACC,color = Identity), position = position_jitter(width = .1), size = 1) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      guides(fill = FALSE) +guides(color = FALSE)+
     # theme_bw() +
      raincloud_theme+scale_y_continuous(limits = c(0,1))+labs(x = "Good",y = "Accuracy")
    
    p2 <- ggplot(data = ACCdata2, aes(y = ACC, x = Identity,fill = Identity)) +
      geom_flat_violin(position = position_nudge(x = .1, y = 0)) +
      geom_point(aes(y = ACC,color = Identity), position = position_jitter(width = .1), size = 1) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      guides(fill = FALSE) +guides(color = FALSE)+
      #theme_bw() +
      raincloud_theme+scale_y_continuous(limits = c(0,1),breaks = NULL)+labs(x = "Bad",y = "")
  }else if(type == 'dprime'){
    Ddata1 <- inData %>%
      select(Subject,Morality,Identity,dprime) %>% 
      filter(Morality == "Good")
    Ddata1$Identity <- factor(Ddata1$Identity,levels = c("Self","Other"))
    p1 <- ggplot(data = Ddata1, aes(y = dprime, x = Identity,fill = Identity)) +
      geom_flat_violin(position = position_nudge(x = .1, y = 0)) +
      geom_point(aes(y = dprime,color = Identity), position = position_jitter(width = .1), size = 1) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      guides(fill = FALSE) +guides(color = FALSE)+
      #theme_bw() +
      raincloud_theme+scale_y_continuous(limits = c(0,4))+labs(x = "Good",y = "dprime")
    Ddata2 <- inData %>%
      select(Subject,Morality,Identity,dprime) %>% 
      filter(Morality == "Bad")
    Ddata2$Identity <- factor(Ddata2$Identity,levels = c("Self","Other"))
    p2 <- ggplot(data = Ddata2, aes(y = dprime, x = Identity,fill = Identity)) +
      geom_flat_violin(position = position_nudge(x = .1, y = 0)) +
      geom_point(aes(y = dprime,color = Identity), position = position_jitter(width = .1), size = 1) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      guides(fill = FALSE) +guides(color = FALSE)+
      # theme_bw() +
      raincloud_theme+scale_y_continuous(limits = c(0,4),breaks = NULL)+labs(x = "Bad",y = "")
    
  }else if(type == 'RT'){
    if(task == 'match'){
      RTdata1 <- inData %>%
        select(Subject,Match,Morality,Identity,RT) %>% 
        filter(Match == "match" & Morality == "Good")
      RTdata2 <- inData %>%
        select(Subject,Match,Morality,Identity,RT) %>% 
        filter(Match == "match" & Morality == "Bad")
    }else if(task == 'val'){
      RTdata1 <- inData %>%
        select(Subject,Task,Morality,Identity,RT) %>% 
        filter(Task == "Val"& Morality == 'Good')
      RTdata2 <- inData %>%
        select(Subject,Task,Morality,Identity,RT) %>% 
        filter(Task == "Val"& Morality == 'Bad')
    }else if (task == 'id'){
      RTdata1 <- inData %>%
        select(Subject,Task,Morality,Identity,RT) %>% 
        filter(Task == "Id"& Morality == 'Good')
      RTdata2 <- inData %>%
        select(Subject,Task,Morality,Identity,RT) %>% 
        filter(Task == "Id"& Morality == 'Bad')
    }else {
      RTdata1 <- inData %>%
        select(Subject,Morality,Identity,RT) %>% 
        filter(Morality == 'Good')
      RTdata2 <- inData %>%
        select(Subject,Morality,Identity,RT) %>% 
        filter(Morality == 'Bad')
    }
    
    RTdata1$Identity <- factor(RTdata1$Identity,levels = c("Self","Other"))
    RTdata2$Identity <- factor(RTdata2$Identity,levels = c("Self","Other"))
    
    p1 <- ggplot(data = RTdata1, aes(y = RT, x = Identity,fill = Identity)) +
      geom_flat_violin(position = position_nudge(x = .1, y = 0)) +
      geom_point(aes(y = RT,color = Identity), position = position_jitter(width = .1), size = 1) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      guides(fill = FALSE) +guides(color = FALSE)+
      # theme_bw() +
      raincloud_theme+scale_y_continuous(limits = c(300,900))+labs(x = "Good",y = "Reaction times(ms)")
    
    p2 <- ggplot(data = RTdata2, aes(y = RT, x = Identity,fill = Identity)) +
      geom_flat_violin(position = position_nudge(x = .1, y = 0)) +
      geom_point(aes(y = RT,color = Identity), position = position_jitter(width = .1), size = 1) +
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
      guides(fill = FALSE) +guides(color = FALSE)+
      #theme_bw() +
      raincloud_theme+scale_y_continuous(limits = c(300,900),breaks = NULL)+labs(x = "Bad",y = "")
  }
  
  fileName = paste0('p_',task,'_',type,'.tiff')
  setwd(saveDir)
  tiff(fileName, width = 9, height = 6, units = 'in', res = 300)
  p_dprime_match <- multiplot(p1,p2,cols = 2)
  dev.off()
  setwd(curDir)
  return(multiplot(p1,p2,cols = 2))
}
