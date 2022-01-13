




### ########################################################## ###
#####               SCRIPT  FOR  VISUALIZATION                ####
### ########################################################## ###

require(ggplot2)
param.z = 30
setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z)))
referencia <- readRDS(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z),"/SCC_vs_SPM.RDS"))
SCC_vs_SPM <- readRDS(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z),"/SCC_vs_SPM_complete.RDS"))
table <- SCC_vs_SPM # For commodity
table$region <- factor(table$region,      
                       levels = c("w32", "w79", "w214", "w271", "w413", "wroiAD"))
table$Roi <- as.numeric(table$Roi)
table$Roi <- table$Roi*10
table$Roi <- as.character(table$Roi)
attach(table)
View(table)


########################################################################################
########################################################################################

# Graph1: Sensibility for Only wroi for zoom analysis:

graph1 <- ggplot(data = table[table$region == "wroiAD", ], 
          aes(x = Roi, y = sens)) + 
          geom_boxplot(aes(fill = method)) + 
          xlab("Level of Induced Hypoactivity") + 
          ylab("Sensibility") + ggtitle(paste0("Comparison of SCC's and SPM's levels of sensibility for wroiAD in z=", as.numeric(param.z))) +
          guides(fill = guide_legend(title = "Legend"))

graph1

graph2 <- ggplot(data = table[table$region == "wroiAD", ], 
          aes(x = Roi, y = esp)) + 
          geom_boxplot(aes(fill = method)) + 
          # coord_cartesian(ylim = c(0, 100)) +
          xlab("Level of Induced Hypoactivity") + 
          ylab("Especificity") + ggtitle(paste0("Comparison of SCC's and SPM's levels of especificity for wroiAD in z=", as.numeric(param.z))) +
          guides(fill = guide_legend(title = "Legend"))

graph2

########################################################################################
########################################################################################

# Graph2: Sensibility for All regions & All ROIs

graph2 <- ggplot(data = table, 
                 aes(x = Roi, y = sens)) + 
                 geom_boxplot(outlier.colour = NULL, aes(fill = method), outlier.size = 1) +
                 xlab("Induced Hypoactivity (%)") + 
                 ylab("Sensibility (%)") + 
                 # ggtitle(paste0("Comparison of sensibility by regions in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Methods")) + 
                 facet_wrap( ~ region, nrow = 3) 

graph2 


graph22 <- ggplot(data = table, 
                 aes(x = Roi, y = esp)) + 
                 geom_boxplot(outlier.colour = NULL, aes(fill = method), outlier.size = 1) +
                 xlab("Induced Hypoactivity (%)") + 
                 ylab("Specificity (%)") + 
                 # ggtitle(paste0("Comparison of especificity by regions in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Methods")) + 
                 facet_wrap( ~ region, nrow = 3) 

graph22 


########################################################################################
########################################################################################

# Academic Publication Theme:

install.packages("envalysis")
library(envalysis)
graph22 + theme_publish() 

ggsave(
  "academic_plot.png",
  plot = last_plot(),
  scale = 2)

ggsave(filename = "academic_plot.png", 
       plot = last_plot(), 
       width = 210, 
       height = 297, 
       units = "mm")



# latex Table:

library(tidyverse)

referencia[referencia$region == "w32", ]$region <- as.character("ROI1")

referencia[referencia$region == "w79", ]$region <- as.character("ROI2")

referencia[referencia$region == "w214", ]$region <- as.character("ROI3")

referencia[referencia$region == "w271", ]$region <- as.character("ROI4")

referencia[referencia$region == "w413", ]$region <- as.character("ROI5")

referencia[referencia$region == "wroiAD", ]$region <- as.character("ROI6")

referencia$roi <- referencia$roi*10
referencia$roi <- as.character(referencia$roi)

for (i in 1:nrow(referencia)) {
referencia$sens[i] <- paste0(round(as.numeric(referencia$sensMEAN[i]), 2), as.character("?"), round(as.numeric(referencia$sensSD[i]),2)) 
}

for (i in 1:nrow(referencia)) {
referencia$esp[i] <- paste0(round(as.numeric(referencia$espMEAN[i]), 2), as.character("?"), round(as.numeric(referencia$espSD[i]),2)) 
}


my_data <- as_tibble(referencia)[ , c(1, 2, 3, 8, 9)]

my_data <- my_data %>% arrange(method, region)

print(xtable(my_data, type = "latex"), file = "my_table.tex")

write.table(my_data, "ourTable.txt", quote=FALSE, eol="\\\\\n", sep=" & ")
