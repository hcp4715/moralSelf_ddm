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

pkgNeeded <- (c("tidyverse", "reshape2", 'plyr', "readr","Hmisc","RColorBrewer"))

lapply(pkgNeeded,pkgTest)
rm('pkgNeeded') # remove the variable 'pkgNeeded';

# Save some time and stor APA format-related code in an object so you can easily
# use it in multiple plots
windowsFonts(Times=windowsFont("TT Times New Roman")) # explicit mapping to "times"

# define the d prime function
dprime <- function(hit,fa) {
        qnorm(hit) - qnorm(fa)
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