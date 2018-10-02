library(readr)
library(cowplot)
library(tidyverse)
library(ggstatsplot)
source("R_rainclouds.R")
#setwd("C:\\Users\\Lanronics\\Desktop\\?½??ļ???\\moral\\exp7r")
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
#y_scale <- scale_y_continuous(limits = range(ACCdata1$ACC))
ACCdata2 <- read_csv("MS_rep_match__rt_acc_long.csv") %>%
  select(Subject,Match,Morality,Identity,ACC) %>% 
  filter(Match == "match" & Morality == "Bad")
ACCdata2$Identity <- factor(ACCdata2$Identity,levels = c("Self","Other"))
p2 <- ggbetweenstats(data = ACCdata2, 
                     x = Identity, 
                     y = ACC,  
                     xlab = "Bad", 
                     ylab = "",
                     type = FALSE) + theme_classic() + theme(legend.position = "none",
                                                             axis.line.x=element_line(size=1.2),
                                                             axis.title.x =element_text(size=14),
                                                             axis.text.x=element_text(size=12,color = "black"),
                                                             axis.ticks.x=element_line(size=1.2)) +scale_y_continuous(breaks = NULL)
y_scale <- scale_y_continuous(limits = range(ACCdata2$ACC))
ACCdata1 <- read_csv("MS_rep_match__rt_acc_long.csv") %>%
  select(Subject,Match,Morality,Identity,ACC) %>% 
  filter(Match == "match" & Morality == "Good")
ACCdata1$Identity <- factor(ACCdata1$Identity,levels = c("Self","Other"))
p1 <- ggbetweenstats(data = ACCdata1, 
                     x = Identity, 
                     y = ACC,  
                     xlab = "Good", 
                     ylab = "Accuracy",
                     type = FALSE) + theme_classic() + theme(legend.position = "none",
                                                             axis.line.x=element_line(size=1.2),
                                                             axis.line.y=element_line(size=1.2),
                                                             axis.title.x =element_text(size=14), 
                                                             axis.title.y=element_text(size=14),
                                                             axis.text.x=element_text(size=12, color = "black"),
                                                             axis.text.y=element_text(size=12, color = "black"),
                                                             axis.ticks.x=element_line(size=1.2),
                                                             axis.ticks.y=element_line(size=1.2))+y_scale

p3 <- multiplot(p1,p2,cols = 2)
#grouped_ggbetweenstats(data = ACCdata1, x = Morality, y = ACC, grouping.var = Identity,ylab = "Reaction times(ms)")
ACCdata <- rbind(ACCdata1,ACCdata2)

#new geom_flat_violin plot
ACCdata$Morality[ACCdata$Morality == "Good"] = 2
ACCdata$Morality[ACCdata$Morality == "Bad"] = 1
p4 <- ggplot(ACCdata,aes(x = Morality, 
                          y = ACC, 
                          fill = Identity))+geom_flat_violin(aes(fill = Identity),
                                                             position = position_nudge(x = 0.1, y = 0),
                                                             adjust = 1.5, trim = FALSE, 
                                                             alpha = 0.5,color = NA)+geom_point(aes(x = as.numeric(Morality)-0.15, 
                                                                                                    y = ACC, 
                                                                                                    color = Identity), 
                                                                                                position = position_jitter(width = 0.02),
                                                                                                size = 1, 
                                                                                                shape = 20)+geom_boxplot(aes(x = Morality, 
                                                                                                                             y = ACC, 
                                                                                                                             fill = Identity),
                                                                                                                         outlier.shape = NA, 
                                                                                                                         alpha = 0.5, 
                                                                                                                         width = 0.1, 
                                                                                                                         color = "black",
                                                                                                                         position = position_dodge(0.15))+scale_color_brewer(palette = "Dark2")+scale_fill_brewer(palette = "Dark2")+ylab("Accuracy")+scale_x_discrete(breaks = c(1,2),labels = c("Bad","Good"))

