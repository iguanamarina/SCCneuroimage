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
install.packages("lemon")
library(envalysis); library(tidyverse); library(lemon)
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

table2 <- rbind(table[table$method == "SCC", ], table[table$method == "SPM", ]) # just two methods (no strongSPM)


####  
# PART 3: Graphs and Plots ----------
####  


#* Graph 1 & Graph 2: Sens & Esp for wROI ----

setwd(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z), "/Figures")) 

graph1 <- ggplot(data = table2[table2$region == "wroiAD", ], aes(x = Roi, y = sens)) + 
          coord_cartesian(ylim = c(0, 100)) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
          geom_boxplot(aes(fill = method)) + 
          xlab("Level of Induced Hypoactivity (%)") + 
          ylab("Sensibility (%)") + 
          # ggtitle(paste0("SCC's and SPM's sensibilities for wroiAD at z=", as.numeric(param.z))) +
          guides(fill = guide_legend(title = "Legend")) +
          scale_fill_brewer(palette = "Set1") 
          
graph2 <- ggplot(data = table2[table2$region == "wroiAD", ], 
          aes(x = Roi, y = esp)) + 
          geom_boxplot(aes(fill = method, col = method)) + 
          guides(col = FALSE) +
          coord_cartesian(ylim = c(0, 100)) + # necessary to specify because in most of cases the scale ranges from 90% to 100% only
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
          xlab("Level of Induced Hypoactivity (%)") + 
          ylab("Specificity (%)") + 
          guides(fill = guide_legend(title = "Legend")) +
          scale_color_brewer(palette = "Set1") +
          scale_fill_brewer(palette = "Set1") 
          
plot1 <- graph1 + theme_publish(); plot2 <- graph2 + theme_publish()
print <- grid.arrange(plot1, plot2, ncol = 2)

ggsave(filename = paste0("sens_esp_wroiAD_", as.numeric(param.z), ".png"), 
       plot = print, 
       width = 24, 
       height = 18, 
       units = "cm",
       dpi = 600)


#* Graph 3 & 4: Sensibility for All regions & All ROIs ----

# The use of the package "lemon" and functions for rewrapping
# keeping axis lines in every graph comes mostly from the following page:
# https://cran.r-project.org/web/packages/lemon/vignettes/facet-rep-labels.html

setwd(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z), "/Figures")) 

graph3 <- ggplot(data = table2, 
                 aes(x = Roi, y = sens)) + 
                 scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                 geom_boxplot(outlier.colour = NULL, aes(fill = method), outlier.size = 1, lwd = 0.25) +
                 xlab("Hypoactivity (%)") + 
                 ylab("Sensitivity (%)") + 
                 # ggtitle(paste0("Sensibility by region in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Legend")) + 
                 facet_wrap( ~ region, ncol = 2) +
                 facet_rep_wrap( ~ region, repeat.tick.labels = TRUE) +
                 # facet_grid(rows= vars(V9), cols = vars(region)) +
                 coord_capped_cart(bottom = 'both', left = 'both') +
                 theme(panel.border = element_blank(), axis.line = element_line()) +
                 theme(legend.text = element_text(size = 14)) +
                 theme_publish(base_size = 12) +
                 scale_fill_brewer(palette = "Set1") 
 
graph3

plot3 <- grid.arrange(graph3)
 
  ggsave(filename = paste0("sens_ALL_", as.numeric(param.z), ".png"), 
         plot = plot3, 
         width = 28.95, 
         height = 18.3, 
         units = "cm",
         dpi = 600)

  
graph4 <- ggplot(data = table2, 
                 aes(x = Roi, y = esp)) + 
                 scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                 geom_boxplot(outlier.colour = NULL, aes(fill = method, col = method), outlier.size = 1, lwd = 0.25) +
                 guides(col = FALSE) +
                 xlab("Hypoactivity (%)") + 
                 ylab("Specificity (%)") + 
                 coord_cartesian(ylim = c(0, 100)) +
                 # ggtitle(paste0("Specificity by region in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Legend")) + 
                 facet_wrap( ~ region, nrow = 3) +
                 scale_fill_brewer(palette = "Set1") +
                 scale_color_brewer(palette = "Set1")

  plot4 <- graph4 + theme_publish(base_size = 8); plot4
  plot4 <- grid.arrange(plot4)
  
  ggsave(filename = paste0("esp_ALL_", as.numeric(param.z), ".png"), 
         plot = plot4, 
         width = 22, 
         height = 26, 
         units = "cm",
         dpi = 600)

  
  plot_together <- grid.arrange(plot3, plot4, ncol = 2)
  
  ggsave(filename = paste0("sens_esp_ALL_", as.numeric(param.z), ".png"), 
         plot = plot_together, 
         width = 22, 
         height = 26, 
         units = "cm",
         dpi = 600)

#* Graph 5 & 6: PPV and NPV for wROI ----

setwd(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z), "/Figures")) 

graph5 <- ggplot(data = table2, 
                 aes(x = Roi, y = ppv)) + 
                 scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                 geom_boxplot(aes(fill = method)) +
                 xlab("Hypoactivity (%)") + 
                 ylab("Positive Predictive Value (%)") + 
                 # ggtitle(paste0("Sensibility by region in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Legend")) + 
                 facet_wrap( ~ region, ncol = 2) +
                 facet_rep_wrap( ~ region, repeat.tick.labels = TRUE) +
                 coord_capped_cart(bottom='both', left='both') +
                 theme(panel.border=element_blank(), axis.line=element_line()) +
                 theme_publish(base_size = 8) +
                 scale_fill_brewer(palette = "Set1")   +
                 coord_cartesian(ylim = c(0, 100)) # I get a warning message related to this
                                                   # temporal solution was to bring it down here

graph5

graph6 <- ggplot(data = table2, 
          aes(x = Roi, y = npv)) + 
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
          coord_cartesian(ylim = c(0, 100)) +
          geom_boxplot(aes(fill = method)) + 
          xlab("Level of Induced Hypoactivity (%)") + 
          ylab("Negative Predicting Value") + 
          ggtitle(" ") +
          guides(fill = guide_legend(title = "Legend")) +
          scale_fill_brewer(palette = "Set1") +
          scale_color_brewer(palette = "Set1")

graph6

plot6 <- graph6 + theme_publish()
graph5 <- grid.arrange(graph5)

  ggsave(filename = paste0("PPV_wroiAD_", as.numeric(param.z), ".png"), 
         plot = graph5, 
         width = 22, 
         height = 26, 
         units = "cm",
         dpi = 600)


#* Graph 7: PPV and NPV All regions & All ROIs ----

graph7 <- ggplot(data = table2, 
                 aes(x = Roi, y = ppv)) +
                 scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                 geom_boxplot(outlier.colour = NULL, aes(fill = method), outlier.size = 1, lwd = 0.25) +
                 xlab("Hypoactivity (%)") + 
                 ylab("Positive Predictive Value (%)") + 
                 coord_cartesian(ylim = c(0, 100)) +
                 # ggtitle(paste0("PPV by region in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Methods")) + 
                 facet_wrap( ~ region, nrow = 3) +
                 scale_fill_brewer(palette = "Set1") +
                 scale_color_brewer(palette = "Set1")

graph7 

graph8 <- ggplot(data = table2, 
                 aes(x = Roi, y = npv)) +
                 scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
                 geom_boxplot(outlier.colour = NULL, aes(fill = method, col = method), outlier.size = 1, lwd = 0.25) +
                 xlab("Hypoactivity (%)") + 
                 ylab("Negative Predictive Value (%)") + 
                 coord_cartesian(ylim = c(0, 100)) +
                 # ggtitle(paste0("NPV by region in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Methods")) + 
                 guides(col = FALSE) +
                 facet_wrap( ~ region, nrow = 3) +
                 scale_fill_brewer(palette = "Set1") +
                 scale_color_brewer(palette = "Set1")

graph8 

plot7 <- graph7 + theme_publish(base_size = 8)
plot7 <- grid.arrange(plot7, ncol = 1)

plot8 <- graph8 +  theme_publish(base_size = 8)
plot8 <- grid.arrange(plot8, ncol = 1)

plot9 <- grid.arrange(plot7, plot8, ncol = 2)

  ggsave(filename = paste0("PPV_ALL_", as.numeric(param.z), ".png"), 
         plot = plot7, 
         width = 22, 
         height = 26, 
         units = "cm",
         dpi = 600)

    ggsave(filename = paste0("NPV_ALL_", as.numeric(param.z), ".png"), 
         plot = plot8, 
         width = 22, 
         height = 26, 
         units = "cm",
         dpi = 600)

    ggsave(filename = paste0("PPV_NPV_ALL_", as.numeric(param.z), ".png"), 
         plot = plot9, 
         width = 22, 
         height = 26, 
         units = "cm",
         dpi = 600)

