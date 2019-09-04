# using ggdag to plot the theoretical models

library(ggdag)
library(tidyverse)

# M1: self and valence are parallel

#M1 <- ggdag::dagify(Perception ~ Self_Relevance,
#                    Perception ~ Valence,
#                    labels = c("Perception" = "percept",
#                               "Self_Relevance" = "self\n relevance",
#                               "Valence" = "valence")) %>%
#        ggdag::ggdag(.,text = FALSE, use_labels = 'label')
        #ggdag::ggdag(M1, text = FALSE, use_labels = 'label')
M1 <- ggdag::collider_triangle(x = "Self_Relevance",
                         y = "Valence",
                         m = "Percept") %>%
        ggdag::ggdag(.,text = FALSE, use_labels = 'label')

# M2: self influence perception via valence
M2 <- ggdag::dagify( Valence ~ Self_Relevance,
                    Perception ~ Valence,
                    labels = c("Perception" = "percept",
                               "Self_Relevance" = "self\n relevance",
                               "Valence" = "valence")) %>%
        ggdag::ggdag(. ,text = FALSE, use_labels = 'label')