rm(list = ls())
Sys.setlocale("LC_ALL", "English")  # set local encoding to English
Sys.setenv(LANG = "en") # set the feedback language to English
options(scipen = 999)   # force R to output in decimal instead of scientifc notion
options(digits=5)       # limit the number of reporting

# get the working directory of the script
curDir <- dirname(rstudioapi::getSourceEditorContext()$path) #Get the directory ofcurrent script
setwd(curDir)

library(ggplot2)
library(bootES)

# APA format-related code for plots
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
            plot.title = element_text(lineheight=.8, face="bold", size = 16,hjust = 0.5),
            axis.text = element_text (size = 16, color = 'black'),
            #              axis.text.x = element_text(angle = 45, vjust = 0.5),   # x-axis's label font
            axis.title = element_text (size = 18),
            axis.title.x = element_text(margin=margin(10,0,0,0)),  # increase the sapce betwen title and x axis
            axis.title.y = element_text(margin=margin(0,10,0,0)),  # increase the space between title and y axis
            axis.line.x = element_line(color='black', size = 1),   # increase the size of font
            axis.line.y = element_line(color='black', size = 1))   # increase the size of font

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

# load data
ct.df <- read.csv('MS_cross_taskeffect_wide_combined.csv',header = T)
bootES(ct.df[,c('d_goodslf_goodoth', 'Id_ACC_goodslf_goodoth')], R = 2000, effect.type = 'r',plot = T)
bootES(ct.df[,c('RT_goodslf_goodoth', 'Id_RT_goodslf_goodoth')], R = 2000, effect.type = 'r',plot = T)
# plot the RT & ACC between matching-task and id-based categorization task
ct.p1 <- ggplot(ct.df,aes(d_goodslf_goodoth,Id_ACC_goodslf_goodoth)) +
      geom_point(colour = 'grey', size = 3) +
      geom_smooth(method="lm", colour = 'black',size = 1.5) +
      ggtitle('Good-self vs. Good-other') +
      scale_x_continuous(name = 'Matching task (d prime)',
                         breaks = seq(-2,4,1),
                         limits = c(-2,4)) +
      scale_y_continuous(name = 'Identity-based categorization (accuracy)',
                         breaks = seq(-0.2,0.2,0.1),
                         limits = c(-0.2,0.2)) +
      apatheme


ct.p2 <- ggplot(ct.df,aes(RT_goodslf_goodoth,Id_RT_goodslf_goodoth)) +
      geom_point(colour = 'grey', size = 3) +
      geom_smooth(method="lm", colour = 'black',size = 1.5) +
      ggtitle('Good-self vs. Good-other') +
      #ylab("Identity-based categorization (RT)") + 
      #xlab() +
      scale_x_continuous(name ='Matching task (RT)', breaks = seq(-200,200,50)) +
      scale_y_continuous(name ='Identity-based categorization (RT)', breaks = seq(-100,150,50)) +
      apatheme

ct.p3 <- ggplot(ct.df,aes(d_badslf_badoth,Id_ACC_badslf_badoth)) +
      geom_point(colour = 'grey', size = 3) +
      geom_smooth(method="lm", colour = 'black',size = 1.5) +
      ggtitle('Bad-self vs. Bad-other') +
      ylab("Identity-based categorization (accuracy)") + 
      #xlab() +
      scale_x_continuous(name = 'Matching task (d prime)', breaks = c(-2,-1,0,1,2,3)) +
      apatheme


ct.p4 <- ggplot(ct.df,aes(RT_badslf_badoth,Id_RT_badslf_badoth)) +
      geom_point(colour = 'grey', size = 3) +
      geom_smooth(method="lm", colour = 'black',size = 1.5) +
      ggtitle('Bad-self vs. Bad-other') +
      scale_x_continuous(name ='Matching task (RT)', breaks = seq(-200,200,50)) +
      scale_y_continuous(name ='Identity-based categorization (RT)', breaks = seq(-100,150,50)) +
     # scale_x_continuous(breaks = c(-2,-1,0,1,2)) +
      apatheme
pdf("scatterplot_crosstask.pdf",width = 16, height = 12)
multiplot(ct.p1,ct.p2,ct.p3,ct.p4,cols = 2)
dev.off()

