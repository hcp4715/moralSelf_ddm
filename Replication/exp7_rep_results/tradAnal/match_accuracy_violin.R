library(tidyverse)
library(ggstatsplot)
setwd("C:\\Users\\Lanronics\\Desktop\\moral\\exp7r")
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
ACCdata1 <- read_csv("exp7_rep_match__rt_acc_long.csv") %>%
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
                                                             axis.ticks.y=element_line(size=1.2))
#y_scale <- scale_y_continuous(limits = range(ACCdata1$ACC))
ACCdata2 <- read_csv("exp7_rep_match__rt_acc_long.csv") %>%
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
multiplot(p1,p2,cols = 2)
#grouped_ggbetweenstats(data = ACCdata1, x = Morality, y = ACC, grouping.var = Identity,ylab = "Reaction times(ms)")
