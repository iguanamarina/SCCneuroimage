############################## ################### ############## ### 
##
## Script name: SCRIPT FOR SENS-ESP TABLE
##
## Purpose of script: To create a table with sensibilities, specificities, PPV, NPV
## and associated means and standard deviations.
##
## Date Created: 2022-01-25
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

library(dplyr)

#* Load up functions: ----

SCC_vs_SPM <- readRDS("~/GitHub/SCCneuroimage/z30/results/SCC_vs_SPM.RDS")


####  
# PART 1: Modify Data Structure ----------
####  

SCC_vs_SPM <- rbind(SCC_vs_SPM[SCC_vs_SPM$method == "SCC", ], SCC_vs_SPM[SCC_vs_SPM$method == "SPM", ]) # just two methods (no strongSPM)

SCC_vs_SPM$region <-   factor(SCC_vs_SPM$region,      
                       levels = c("w32", "w79", "w214", "w271", "w413", "wroiAD"))

SCC_vs_SPM$roi <- as.numeric(SCC_vs_SPM$roi)
SCC_vs_SPM$roi <- SCC_vs_SPM$roi*10
SCC_vs_SPM$roi <- as.character(SCC_vs_SPM$roi)


####  
# PART 2: dplyr ----------
####  

  
table <- SCC_vs_SPM %>%
  mutate(Sensibility = paste(round(sensMEAN, 2), round(sensSD, 2), sep = "\u00B1")) %>%
  mutate(Specificity = paste(round(espMEAN, 2), round(espSD, 2), sep = "\u00B1")) %>%
  mutate(PPV = paste(round(ppvMEAN, 2), round(ppvSD, 2), sep = "\u00B1")) %>%
  mutate(NPV = paste(round(npvMEAN, 2), round(npvSD, 2), sep = "\u00B1")) %>%
  arrange(method, roi, region) %>%
  dplyr::select(method, roi, region, Sensibility, Specificity, PPV, NPV)
table

