############################## ################### ############## ### 
##
## Script name: SCRIPT FOR GRAPHS AND VISUALIZATION
##
## Purpose of script: Visualize data obtained from MASTERSCRIPT of the
## SCC Neuroimage project.
##
## Date Created: 2022-01-13
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://messy-dataset.xyz
##
## Notes:
##       
##   
############################## ################### ############## ### 

####  
# PREAMBLE: ----
#### 

#* Set working directory: ----

setwd("~/GitHub/SCCneuroimage") 

#* Tune Options: ----
options(scipen = 6, digits = 4) # View outputs in non-scientific notation
memory.limit(30000000)     # This is needed on some PCs to increase memory allowance

#* Load up packages: ---- 

# install.packages("envalysis")
library(envalysis); library(tidyverse)
require(ggplot2); require(gridExtra)

#* Hyperparameter: ----

param.z = 30

####  
# PART 1: Load Data ----------
####  

referencia <- readRDS(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z),"/results/SCC_vs_SPM.RDS"))
table <- readRDS(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z),"/results/SCC_vs_SPM_complete.RDS"))
  
  
####  
# PART 2: Modify Data Structure ----------
####  

table$region <- factor(table$region,      
                       levels = c("w32", "w79", "w214", "w271", "w413", "wroiAD"))
table$Roi <- as.numeric(table$Roi)
table$Roi <- table$Roi*10
table$Roi <- as.character(table$Roi)
attach(table)
View(table) # Quite beautiful


####  
# PART 3: Graphs and Plots ----------
####  


#* Graph1 & Graph 2: Sens & Esp for wROI ----

graph1 <- ggplot(data = table[table$region == "wroiAD", ], 
          aes(x = Roi, y = sens)) + 
          coord_cartesian(ylim = c(0, 100)) +
          geom_boxplot(aes(fill = method)) + 
          xlab("Level of Induced Hypoactivity") + 
          ylab("Sensibility") + ggtitle(paste0("Comparison of SCC's and SPM's levels of sensibility for wroiAD in z=", as.numeric(param.z))) +
          guides(fill = guide_legend(title = "Legend"))

graph2 <- ggplot(data = table[table$region == "wroiAD", ], 
          aes(x = Roi, y = esp)) + 
          geom_boxplot(aes(fill = method)) + 
          coord_cartesian(ylim = c(0, 100)) + # necessary to specify because in most of cases the scale ranges from 90% to 100% only
          xlab("Level of Induced Hypoactivity") + 
          ylab("Especificity") + ggtitle(paste0("Comparison of SCC's and SPM's levels of especificity for wroiAD in z=", as.numeric(param.z))) +
          guides(fill = guide_legend(title = "Legend"))

plot1 <- graph1 + theme_publish(); plot2 <- graph2 + theme_publish()
grid.arrange(plot1, plot2, ncol = 2)

ggsave(filename = "academic_plot.png", 
       plot = last_plot(), 
       width = 210, 
       height = 297, 
       units = "mm")


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


