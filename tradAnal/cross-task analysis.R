library(tidyverse)
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

# read data
M_ACC_data <- read_csv("exp7_rep_match__dprime_long.csv") %>%
  select(Subject,Morality,Identity,dprime) %>%
  unite(Type,Morality,Identity) %>%
  spread(key = Type, value = dprime) %>%
  mutate("self-referential" = Good_Self - Good_Other, "valence" = Good_Self - Bad_Self)
M_RT_data<- read_csv("exp7_rep_match__rt_acc_long.csv") %>%
  filter(Match == "match") %>%
  select(Subject,Morality,Identity,RT) %>%
  unite(Type,Morality,Identity) %>%
  spread(key = Type, value = RT) %>%
  mutate("self-referential" = Good_Self - Good_Other, "valence" = Good_Self - Bad_Self)
C_Id_RT_data <- read_csv("exp7_rep_categ__rt_acc_long.csv") %>%
  filter(Task == "Id") %>%
  select(Subject,Morality,Identity,RT) %>%
  unite(Type,Morality,Identity) %>%
  spread(key = Type, value = RT) %>%
  mutate("self-referential" = Good_Self - Good_Other, "valence" = Good_Self - Bad_Self)
C_Id_ACC_data <- read_csv("exp7_rep_categ__rt_acc_long.csv") %>%
  filter(Task == "Id") %>%
  select(Subject,Morality,Identity,ACC) %>%
  unite(Type,Morality,Identity) %>%
  spread(key = Type, value = ACC) %>%
  mutate("self-referential" = Good_Self - Good_Other, "valence" = Good_Self - Bad_Self)
C_Val_RT_data <- read_csv("exp7_rep_categ__rt_acc_long.csv") %>%
  filter(Task == "Val") %>%
  select(Subject,Morality,Identity,RT) %>%
  unite(Type,Morality,Identity) %>%
  spread(key = Type, value = RT) %>%
  mutate("self-referential" = Good_Self - Good_Other, "valence" = Good_Self - Bad_Self)
C_Val_ACC_data <- read_csv("exp7_rep_categ__rt_acc_long.csv") %>%
  filter(Task == "Val") %>%
  select(Subject,Morality,Identity,ACC) %>%
  unite(Type,Morality,Identity) %>%
  spread(key = Type, value = ACC) %>%
  mutate("self-referential" = Good_Self - Good_Other, "valence" = Good_Self - Bad_Self)

p1 <- ggplot()+geom_point(aes(M_ACC_data$`self-referential`,C_Val_ACC_data$`self-referential`))+geom_smooth(aes(M_ACC_data$`self-referential`,C_Val_ACC_data$`self-referential`),method = 'lm')+labs(x = "Matching task", y = "Accuracy\n\nMoral categorization task",title = "(A) Effect of self bias")+theme_classic()
cor1 <- cor.test(M_ACC_data$`self-referential`,C_Val_ACC_data$`self-referential`,alternative = "greater")
p2 <- ggplot()+geom_point(aes(M_ACC_data$`self-referential`,C_Id_ACC_data$`self-referential`))+geom_smooth(aes(M_ACC_data$`self-referential`,C_Id_ACC_data$`self-referential`),method = 'lm')+labs(x = "Matching task", y = "\n\n\nIdentity categorization task",title = " ")+theme_classic()
cor2 <- cor.test(M_ACC_data$`self-referential`,C_Id_ACC_data$`self-referential`,alternative = "greater")
p3 <- ggplot()+geom_point(aes(M_RT_data$`self-referential`,C_Val_RT_data$`self-referential`))+geom_smooth(aes(M_RT_data$`self-referential`,C_Val_RT_data$`self-referential`),method = 'lm')+labs(x = "Matching task", y = "Reaction times\n\nMoral categorization task")+theme_classic()
cor3 <- cor.test(M_RT_data$`self-referential`,C_Val_RT_data$`self-referential`,alternative = "greater")
p4 <- ggplot()+geom_point(aes(M_RT_data$`self-referential`,C_Id_RT_data$`self-referential`))+geom_smooth(aes(M_RT_data$`self-referential`,C_Id_RT_data$`self-referential`),method = 'lm')+labs(x = "Matching task", y = "\n\n\nIdentity categorization task")+theme_classic()
cor4 <- cor.test(M_RT_data$`self-referential`,C_Id_RT_data$`self-referential`,alternative = "greater")
p5 <- ggplot()+geom_point(aes(M_ACC_data$valence,C_Val_ACC_data$valence))+geom_smooth(aes(M_ACC_data$valence,C_Val_ACC_data$valence),method = 'lm')+labs(x = "Matching task", y = "Accuracy\n\nMoral categorization task",title = "(B) Effect of valence bias")+theme_classic()
cor5 <- cor.test(M_ACC_data$valence,C_Val_ACC_data$valence,alternative = "greater")
p6 <- ggplot()+geom_point(aes(M_ACC_data$valence,C_Id_ACC_data$valence))+geom_smooth(aes(M_ACC_data$valence,C_Id_ACC_data$valence),method = 'lm')+labs(x = "Matching task", y = "\n\n\nIdentity categorization task",title = " ")+theme_classic()
cor6 <- cor.test(M_ACC_data$valence,C_Id_ACC_data$valence,alternative = "greater")
p7 <- ggplot()+geom_point(aes(M_RT_data$valence,C_Val_RT_data$valence))+geom_smooth(aes(M_RT_data$valence,C_Val_RT_data$valence),method = 'lm')+labs(x = "Matching task", y = "Reaction times\n\nMoral categorization task")+theme_classic()
cor7 <- cor.test(M_RT_data$valence,C_Val_RT_data$valence,alternative = "greater")
p8 <- ggplot()+geom_point(aes(M_RT_data$valence,C_Id_RT_data$valence))+geom_smooth(aes(M_RT_data$valence,C_Id_RT_data$valence),method = 'lm')+labs(x = "Matching task", y = "\n\n\nIdentity categorization task")+theme_classic()
cor8 <- cor.test(M_RT_data$valence,C_Id_RT_data$valence,alternative = "greater")

multiplot(p1,p3,p2,p4,cols = 2)
# draw self-reference plot
multiplot(p5,p7,p6,p8,cols = 2)
# draw moral valence plot