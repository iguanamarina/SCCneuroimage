############################## ################### ############## ### 
##
## Script name: MASTERSCRIPT (for simulations)
##
## Purpose of script: Script for the calculation of sensibility, specificity, predictive value... 
## and other measures for comparing SCC and SPM using simulated data.
##
## Date Created: 2022-01-10
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://messy-dataset.xyz
##
## Notes: This script is self-contained and nothing should be necessary but for
## pressing enter again and again (after hyperparameter setup). However, results
## may differ with slight changes.
##   
############################## ################### ############## ### 


####  
# PREAMBLE: ----
#### 


#* Set working directory ----

setwd("~/GitHub/SCCneuroimage")

#* Tune Options ----
options(scipen = 6, digits = 4) # View outputs in non-scientific notation
memory.limit(30000000)     # This is needed on some PCs to increase memory allowance

#* Install Packgs ----

# install.packages("remotes")
# remotes::install_github("skgrange/threadr")

#* Load up packages ---- 

library(gamair);library(oro.nifti);library(memisc);library(devtools);library(remotes);library(readr);library(imager);library(itsadug);library(ggplot2);library(contoureR);library(fields);library(BPST);library(Triangulation);library(ImageSCC); library(tidyr); library(dplyr);library(stringr); library(threadr);library(memisc)

#* Load up functions ----

load("~/GitHub/SCCneuroimage/Functions/f.clean.RData") # Function for NiFTi -> df
load("~/GitHub/SCCneuroimage/Functions/my_points.RData") # Function for getting relevant points from SCC

#* Hyper-parameters ----

param.z = 35


####  
# PART 1: CONTOURS OF NEURO-DATA ----------
####


#* Basic data for contours ----
  
Data = f.clean("new_mask") 
Y <- subset(Data, Data$z == param.z) 
Y <- Y[1:9919,4] 
Y <- as.matrix(Y)
Y = t(Y) 
Y[is.nan(Y)] <- 0
SCC <- Y; rm(Y); rm(Data)


#* Coordinates ---- 

# These are adapted to my current nifti size, in the future it would be strategic to develop
# this code for a more generic setup

x <- rep(1:91, each = 109, length.out = 9919) 
y <- rep(1:109, length.out = 9919)
Z <- cbind(as.matrix(x),as.matrix(y)); Z
dat <- cbind(Z,t(SCC))
dat <- as.data.frame(dat)
dat[is.na(dat)] <- 0
sum(is.na(dat$pet)) # should be = 0
rownames(dat) <- NULL 
rm(x); rm(y)

# Get contour for the area where values change from 0 to 1 (it's a template)
df = getContourLines(dat[1:9919,], levels = c(0)) 
ggplot(df, aes(x, y, colour = z)) + geom_path() 
contour = df 
rm(df); rm(dat)

f.contour <- function(x){
  
  aa <- contour[contour$GID == x,] # We keep GID==x (0,1,2,3...)
  a <- aa[,5:6] # Then we keep just the coordinates 
  print(a) # and then print in order to loop and make a list
}

coord <- list()

for (i in 0:max(contour$GID)) { #change contour to any other name previously assigned if necessary
   
  coord[[i + 1]] <- f.contour(i)
  rownames(coord[[ i + 1 ]]) <- NULL
}


### Test the results are coherent:

head(coord[[1]],10); plot(coord[[1]])   # external boundaries
head(coord[[2]],10); points(coord[[2]]) # first hole
head(coord[[3]],10); points(coord[[3]]) # second hole
# (...) 

# Most of times, if there is no 'holes' in that brain slice, some of these lines, loops and stuff
# won't be necessary.


####
# PART 2: TRIANGULATION PARAMETERS ----------
####

#* Get coordinates in pckg Triangulation format: ----

VT = TriMesh(coord[[1]], n = 8) 

# n = Triangulation degree of fineness (8 is recommended)
# However, higher values can be used as Arias-López et al. (2021) suggests computing times for
# higher n values are still low and sensible.

head(VT$V,10);head(VT$Tr,10) 

# Create it first if it doesn't exist
setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z)))
save(VT, file = paste0("contour", as.character(param.z), ".RData"))
     

####
# PART 3: CREATE SCC MATRIXES FOR CONTROLS ------
####


# *Load data and Create DB ----

setwd("~/GitHub/SCCneuroimage/PETimg_masked for simulations")

number <- paste0("C", 1:25)
name <- paste0("masked_swww", number, "_tripleNormEsp_w00_rrec_OSEM3D_32_it1")
# Only 25 files are controls and they follow the above defined structure

database_CN <- data.frame(CN_number = integer(),z = integer(), x = integer(), y = integer(), pet = integer())

for (i in 1:length(name)) {
  
  temporal <- f.clean(name[i])
  CN_number <- rep(number[i], length.out = nrow(Z))
  temporal <- cbind(CN_number,temporal)
  database_CN <- rbind(database_CN,temporal)
}

nrow(database_CN[database_CN$pet < 0, ]) # No negative values
rm(temporal); rm(CN_number)


# *Create SCC Matrix ----

# Working on a functional data setup requires for the data to be in a concrete format which
# usually implies a long line of data points so that each row represents a function

SCC_CN <- matrix(nrow = length(name), ncol = nrow(Z))

for (i in 1:length(number)) {

  Y <- subset(database_CN, database_CN$CN_number == number[i] & database_CN$z == param.z) 
  Y <- Y[1:9919, 5] 
  Y <- as.matrix(Y)
  Y = t(Y) 
  Y[is.nan(Y)] <- 0
  SCC_CN[i, ] <- Y
  
}

# Sometimes R really doesn't want to remove NA so this might be necessary:

  # na.zero <- function(x) {
  #     x[is.na(x)] <- 0
  # }
  # SCC_CN <-  apply(SCC_CN, 2, na.zero)


setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z)))
save(SCC_CN, file = "SCC_CN.RData") # SCC matrix for Controls


####
# PART 4: CREATE SCC MATRIXES FOR PATHOLOGICAL -------
####

# This section is pretty much the same as the section on Controls, but now with much more info


# *Load data and Create DB ----

setwd("~/GitHub/SCCneuroimage/PETimg_masked for simulations")
number <- paste0("C", 1:25)
region <- c("roiAD", "w32", "w79", "w214", "w271", "w413")
roi <- c(1, 2, 4, 6, 8)

for (i in 1:length(region)) {
  
  for (j in 1:length(roi)) {
    
    setwd("~/GitHub/SCCneuroimage/PETimg_masked for simulations")
    
    # IMPORTANTE: tuning parameters
    
    REGION = region[i]
    ROI = roi[j] 
    name <- paste0("masked_swww", number, "_tripleNormEsp_", REGION, "_0_", ROI, "_rrec_OSEM3D_32_it1")
    
    database_AD <- data.frame(AD_number = integer(),z = integer(), x = integer(), y = integer(), pet = integer())
    
    for (k in 1:length(name)) {
      
      temporal <- f.clean(name[k])
      AD_number <- rep(paste0(number[k], "_", as.character(REGION), "_", as.character(ROI)), length.out = nrow(Z))
      temporal <- cbind(AD_number, temporal)
      database_AD <- rbind(database_AD, temporal)
    }

# *Create SCC Matrix ---- 
       
    SCC_matrix <- matrix(nrow = length(name), ncol = nrow(Z))
    
    for (k in 1:length(number)) {
    
      Y <- subset(database_AD, AD_number == paste0(number[k], "_", as.character(REGION), "_", as.character(ROI)) &
                  z == param.z)
      Y <- Y[1:9919, 5] 
      Y <- as.matrix(Y)
      Y = t(Y) 
      Y[is.nan(Y)] <- 0
      SCC_matrix[k, ] <- Y
    }
    
    setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z)))
    save(SCC_matrix, file = paste0("SCC_", REGION, "_", ROI, ".RData"))

  }
}


####
# PART 5: COMPUTE SCCs PER SE -------
####


#* Preliminary stuff ----

setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z)))

# In order to be consistent we use common names Brain.V and Brain.Tr. 
# From here onwards most of the names follow the ones provided by Wang et al (2019)

Brain.V <- VT[[1]]
Brain.Tr <- VT[[2]]

V.est = as.matrix(Brain.V)
# Brain.v <- cbind(Brain.V[,2],Brain.V[,1]) # In case you need to transpose the data
Tr.est = as.matrix(Brain.Tr)
V.band = as.matrix(Brain.V)
Tr.band = as.matrix(Brain.Tr) 


#* The Loop: Fasten your seat belts ----

for (i in 1:length(region)) {
  
  for (j in 1:length(roi)) {
    
    setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z)))
    
    # Response Variable:
    
    name_CN <- paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/SCC_CN.RData")
    name_AD <- paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/SCC_", 
                      region[i], "_", roi[j], ".RData")

    SCC_CN <- threadr::read_rdata(name_CN)
    SCC_AD <- threadr::read_rdata(name_AD)
    
    #** Mean Average Normalization:  ----
    
    for (k in 1:nrow(SCC_CN)) {
  
    temp <- SCC_CN[k, ]
    mean <- mean(as.numeric(temp), na.rm = T)
    SCC_CN[k, ] <- (temp/mean)
    }

    for (k in 1:nrow(SCC_AD)) {
    
    temp <- SCC_AD[k, ]
    mean <- mean(as.numeric(temp), na.rm = T)
    SCC_AD[k, ] <- (temp/mean)
    }
    
    #** Other parameters for SCC computation ----
    # Following Wang et al's recomendations:
    
    d.est = 5 # degree of spline for mean function  5
    d.band = 2 # degree of spline for SCC  2
    r = 1 # smoothing parameter  1
    lambda = 10^{seq(-6,3,0.5)} # penalty parameters
    alpha.grid = c(0.10,0.05,0.01) # vector of confidence levels

    #** Construction of SCC's ----
    
    setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/results"))
    
    if (file.exists(paste0("SCC_COMP_", region[i], "_", roi[j] ,".RData")) == TRUE)  {
      
      print("Nice!") # I already have this files so this way the lines below (very time-consuming)
                     # run only if these files are not present at current folder
      
    } else  {
    
      SCC_COMP = scc.image(Ya = SCC_AD, Yb = SCC_CN, Z = Z, d.est = d.est, d.band = d.band, r = r,
                           V.est.a = V.est, Tr.est.a = Tr.est,
                           V.band.a = V.band, Tr.band.a = Tr.band,
                           penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                           adjust.sigma = TRUE)    
      save(SCC_COMP, file = paste0("SCC_COMP_", region[i], "_", roi[j] ,".RData"))
      
    }
    
  }
  
}


####
# PART 6: SCC EVALUATION ----
####


# This section requires loading the TRUE points with differences in PET activity to test
# against them. These are called ROI_something. Basically, we will get the TRUE POINTS
# and then test them against SCC, SPM, SPMstrong... and get some metrics on each method's afficiency.
# Enjoy.


#* Theoretical ROIs ----

region <- c("w32", "w79", "w214", "w271", "w413", "wroiAD")

# With this loop we get a series of tables with the exact ROI points as provided by
# Santiago de Compostela's General Hospital:

for (i in 1:length(region)) { 

    for (k in 1:length(number)) {

      setwd("~/GitHub/SCCneuroimage/roisNormalizadas/tables")
      
      if (file.exists(paste0("ROItable_", region[i], "_", number[k], ".RDS")) == TRUE)  {
         print("Nice!")
        } else  {
          
      setwd("~/GitHub/SCCneuroimage/roisNormalizadas")
      roi_table <- data.frame(group = integer(),
                              z = integer(),
                              x = integer(),
                              y = integer())
      
      # Name of file:
      name <- paste0("wwwx", region[i], "_redim_crop_squ_flipLR_newDim_", number[k])
      # Clean data for that file:
      temporal <- f.clean(name)
      # Add meta-data:
      group <- rep(paste0(as.character(region[i]), "_number", number[k]), length.out = nrow(Z))
      # Merge with main data.frame:
      temporal <- cbind(group, temporal)
      # Solve problem with zeros and NA's:
      suelto <- temporal[, "pet"]
      suelto[is.nan(suelto)] <- 0
      temporal[, "pet"] <- as.data.frame(suelto)
      # Final:
      roi_table <- rbind(roi_table, temporal)

      setwd("~/GitHub/SCCneuroimage/roisNormalizadas/tables")
      saveRDS(roi_table, file = paste0("ROItable_", region[i], "_", number[k], ".RDS"))
        
      }
      
    }
  
  }


#* Theoretical Points (True according to ROIs) ----
  
ROI_data <- list() # First load data

  for (i in 1:length(region)) {
    for (k in 1:length(number)) {
      
    setwd("~/GitHub/SCCneuroimage/roisNormalizadas/tables")
    ROI_data[[paste0(as.character(region[i]), "_", as.character(number[k]))]] <- readRDS(
      paste0("~/GitHub/SCCneuroimage/roisNormalizadas/tables/ROItable_", 
             region[i], "_", number[k], ".RDS"))
  
    }
  }  


T_points <- list() # Now filter by Z 

for (i in 1:length(region)) {
  
  for (j in 1:length(number)) {
    
    T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]] <- ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]][ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]]$z == param.z & ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]]$pet == 1, 3:4]
    
    T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]] <- unite(as.data.frame(T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]]), newcol, c(y,x), remove = T)
        
    }
  }
    
rm(ROI_data) # No longer necessary


#* Hypothetical points according to SCCs ---- 


  # Hypothetical points and sens/esp for different methods go all together
  # in a big loop. First, you need to perform the analysis with SPM and all the 
  # different setups and export SPM results to a binary file named "binary.nii"
  # Pay attention to the folders so that you understand where things come from.


SCC_vs_SPM_complete <-   data.frame(method = integer(),
                                    region = integer(),
                                    roi = integer(),
                                    number = integer(),
                                    sens = integer(),
                                    esp = integer(),
                                    PPV = integer(),
                                    NPV = integer())

SCC_vs_SPM <- data.frame(method = integer(),
                         region = integer(),
                         roi = integer(),
                         sensMEAN = integer(),
                         sensSD = integer(),
                         espMEAN = integer(),
                         espSD = integer(),
                         ppvMEAN = integer(),
                         ppvSD = integer(),
                         npvMEAN = integer(),
                         npvSD = integer())

# These above are the final tables

x <- rep(1:91, each = 109, length.out = 9919) 
y <- rep(1:109, length.out = 9919)
total.coords <- data.frame(y, x) ###!!!!!
total.coords <- unite(as.data.frame(total.coords), newcol, c(y, x), remove = T)
rm(x); rm(y)


# NOW FOR REAL, FASTEN YOUR SEAT BELT, THIS IS GOING TO BE AS HARD AND FAST AS A TECHNNO RAVE 

for (k in 1:length(roi)) {
  
  region <- c("w32", "w79", "w214", "w271", "w413", "roiAD")
  
  H_points <- list()
  
  setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/results"))
  
  for (i in 1:length(region)) {
  
    load(paste0("SCC_COMP_", region[i], "_", roi[k], ".RData"))
    H_points[[as.character(region[i])]] <- my_points(SCC_COMP, 2)[[1]] # 2 = alpha 0'95
    H_points[[as.character(region[i])]] <- unite(as.data.frame(H_points[[i]]), newcol, c(row, col), remove = T)
    
  }
     
  rm(list = ls(pattern = "^SCC_COMP")) # Just beautiful
  
  #* SENS/ESP for SCC ----

  region <- c("w32", "w79", "w214", "w271", "w413", "wroiAD")
  names(H_points)[6] <- "wroiAD" # because of reasons
  
  for (i in 1:length(region)) {  
    
    SCC_sens_esp <- data.frame(region = integer(), group = integer(), sens = integer(), esp = integer(), ppv = integer(), npv = integer())
    
    for (j in 1:length(number)) {
    
      inters <- inner_join(H_points[[as.character(region[i])]],
                           T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])
      
      sensibilitySCC <- nrow(inters)/nrow(T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])*100 
      # p(correctly identify hypoactive pixel)
      
      true_neg <- anti_join(total.coords,
                            T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])
      hypo_neg <- anti_join(total.coords,
                            H_points[[as.character(region[i])]])
      
      anti_inters <- inner_join(true_neg, hypo_neg)
      
      specificitySCC <- nrow(anti_inters)/nrow(true_neg)*100
      # p(correctly identify healthy pixel)
    
      # PPV	= TP/(TP+FP) -> probability of having the disease after a positive test result
      
      FalsePositive <- inner_join(H_points[[as.character(region[i])]], 
                                  true_neg)
      FalseNegative <- inner_join(hypo_neg, 
                                  T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])  
      
      PPV = (nrow(inters)/(nrow(inters) + nrow(FalsePositive)))*100
      
      # NPV	= TN/(FN+TN) -> probability of not having the disease after a negative test result
      
      NPV = (nrow(anti_inters)/(nrow(anti_inters) + nrow(FalseNegative)))*100
      
      temp <- data.frame(region = region[i], group = number[j], sens = sensibilitySCC, esp = specificitySCC, ppv = PPV, npv = NPV)  
      SCC_sens_esp <- rbind(SCC_sens_esp, temp)
      
    }
    
    means <- data.frame(region = region[i], group = "MEAN", sens = mean(SCC_sens_esp$sens, na.rm = TRUE), 
                                                            esp = mean(SCC_sens_esp$esp, na.rm = TRUE),
                                                            ppv = mean(SCC_sens_esp$ppv, na.rm = TRUE),
                                                            npv = mean(SCC_sens_esp$npv, na.rm = TRUE))
    sds <- data.frame(region = region[i], group = "SD", sens = sd(SCC_sens_esp$sens, na.rm = TRUE), 
                                                        esp = sd(SCC_sens_esp$esp, na.rm = TRUE),
                                                        ppv = sd(SCC_sens_esp$ppv, na.rm = TRUE),
                                                        npv = sd(SCC_sens_esp$npv, na.rm = TRUE))
    SCC_sens_esp <- rbind(SCC_sens_esp, means, sds)
    
    setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/results", "/ROI", roi[k]))
    write_csv(SCC_sens_esp, paste0("sens_esp_SCC_", region[i], "_", roi[k], ".csv"), na = "NA", append = FALSE)
    
  }
  
  setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/results", "/ROI", roi[k]))
  
  sens_esp_SCC_w32 <- read_csv(paste0("sens_esp_SCC_w32_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_w79 <- read_csv(paste0("sens_esp_SCC_w79_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_w214 <- read_csv(paste0("sens_esp_SCC_w214_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_w271 <- read_csv(paste0("sens_esp_SCC_w271_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_w413 <- read_csv(paste0("sens_esp_SCC_w413_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_wroiAD <- read_csv(paste0("sens_esp_SCC_wroiAD_", roi[k], ".csv"), na = "NA")
 
  #* SENS/ESP for SPM ----
  
  # This part requires you to have previously performed these analysis in SPM
  # using an uncorrected p-value of 0'05 and no pixel threshold and exporting
  # a binary file ("binary.nii") into folders following the below defined name format
  
  for (i in 1:length(region)) {  
    
    setwd(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z), "/SPM", "/ROI", i, "_", region[i], "_0", roi[k]))
    binary <- f.clean("binary.nii")
    H_points_SPM <- binary[binary$z == as.numeric(param.z) & binary$pet == 1, 2:3]
    H_points_SPM <- unite(as.data.frame(H_points_SPM), newcol, c(y, x), remove = T)
    
    SPM_sens_esp <- data.frame(region = integer(), group = integer(), sens = integer(), esp = integer(), ppv = integer(), npv = integer())
    
    for (j in 1:length(number)) {
    
      inters <- inner_join(H_points_SPM,
                           T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])
    
      sensibilitySPM <- nrow(inters)/nrow(T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])*100 
      # p(correctly identify hypoactive pixel)
      
      true_neg <- anti_join(total.coords,
                            T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])
      hypo_neg <- anti_join(total.coords,
                            H_points_SPM)
      
      anti_inters <- inner_join(true_neg, hypo_neg)
      
      specificitySPM <- nrow(anti_inters)/nrow(true_neg)*100
      # p(correctly identify healthy pixel)
      
      # PPV	= TP/(TP+FP) -> probability of having the disease after a positive test result
      
      FalsePositive <- inner_join(H_points_SPM, 
                                  true_neg)
      FalseNegative <- inner_join(hypo_neg, 
                                  T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])  
      
      PPV = (nrow(inters)/(nrow(inters) + nrow(FalsePositive)))*100
      
      # NPV	= TN/(FN+TN) -> probability of not having the disease after a negative test result
      
      NPV = (nrow(anti_inters)/(nrow(anti_inters) + nrow(FalseNegative)))*100
      
      temp <- data.frame(region = region[i], group = number[j], sens = sensibilitySPM, esp = specificitySPM, ppv = PPV, npv = NPV)  
      SPM_sens_esp <- rbind(SPM_sens_esp, temp)
  
    }
    
    means <- data.frame(region = region[i], group = "MEAN", sens = mean(SPM_sens_esp$sens, na.rm = TRUE), 
                                                            esp = mean(SPM_sens_esp$esp, na.rm = TRUE),
                                                            ppv = mean(SPM_sens_esp$ppv, na.rm = TRUE),
                                                            npv = mean(SPM_sens_esp$npv, na.rm = TRUE))
    sds <- data.frame(region = region[i], group = "SD", sens = sd(SPM_sens_esp$sens, na.rm = TRUE), 
                                                        esp = sd(SPM_sens_esp$esp, na.rm = TRUE),
                                                        ppv = sd(SPM_sens_esp$ppv, na.rm = TRUE),
                                                        npv = sd(SPM_sens_esp$npv, na.rm = TRUE))
    SPM_sens_esp <- rbind(SPM_sens_esp, means, sds)
  
    setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/results", "/ROI", roi[k]))
    
    write_csv(SPM_sens_esp, paste0("sens_esp_SPM_", region[i], "_", roi[k], ".csv"), na = "NA", append = FALSE)
  
  }
  
  setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/results", "/ROI", roi[k]))
  
  sens_esp_SPM_w32 <- read_csv(paste0("sens_esp_SPM_w32_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_w79 <- read_csv(paste0("sens_esp_SPM_w79_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_w214 <- read_csv(paste0("sens_esp_SPM_w214_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_w271 <- read_csv(paste0("sens_esp_SPM_w271_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_w413 <- read_csv(paste0("sens_esp_SPM_w413_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_wroiAD <- read_csv(paste0("sens_esp_SPM_wroiAD_", roi[k], ".csv"), na = "NA")
  
  

     
  #* SENS/ESP for strong SPM (FWE 0'05 + 100 voxels threshold) ----
  
  # Same as the previous section, here you need to carry out SPM analysis 
  # by yourself and then save the resulting binary file as "binary.nii" in 
  # folders following below-shown path structure. In this case, SPM analysis
  # is performed with FWE 0'05 and 100 voxels of threshold.
  
  
  for (i in 1:length(region)) {  
    
    setwd(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z), "/SPMstrong", "/ROI", i, "_", region[i], "_0", roi[k]))
    binary <- f.clean("binary.nii")
    H_points_SPMstrong <- binary[binary$z == as.numeric(param.z) & binary$pet == 1, 2:3]
    H_points_SPMstrong <- unite(as.data.frame(H_points_SPMstrong), newcol, c(y, x), remove = T)
    
    SPMstrong_sens_esp <- data.frame(region = integer(), group = integer(), sens = integer(), esp = integer(), ppv = integer(), npv = integer())
    
    for (j in 1:length(number)) {
    
      inters <- inner_join(H_points_SPMstrong,
                           T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])
    
      sensibilitySPMstrong <- nrow(inters)/nrow(T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])*100 
      
      true_neg <- anti_join(total.coords,
                            T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])
      hypo_neg <- anti_join(total.coords,
                            H_points_SPMstrong)
      
      anti_inters <- inner_join(true_neg, hypo_neg)
      
      specificitySPMstrong <- nrow(anti_inters)/nrow(true_neg)*100
      # p(correctly identify healthy pixel)
      
      # PPV	= TP/(TP+FP) -> probability of having the disease after a positive test result
      
      FalsePositive <- inner_join(H_points_SPMstrong, 
                                  true_neg)
      FalseNegative <- inner_join(hypo_neg, 
                                  T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]])  
      
      PPV = (nrow(inters)/(nrow(inters) + nrow(FalsePositive)))*100
      
      # NPV	= TN/(FN+TN) -> probability of not having the disease after a negative test result
      
      NPV = (nrow(anti_inters)/(nrow(anti_inters) + nrow(FalseNegative)))*100
      
      
      temp <- data.frame(region = region[i], group = number[j], sens = sensibilitySPMstrong, esp = specificitySPMstrong, ppv = PPV, npv = NPV)  
      SPMstrong_sens_esp <- rbind(SPMstrong_sens_esp, temp)
  
    }
    
    means <- data.frame(region = region[i], group = "MEAN", sens = mean(SPMstrong_sens_esp$sens, na.rm = TRUE), 
                                                            esp = mean(SPMstrong_sens_esp$esp, na.rm = TRUE),
                                                            ppv = mean(SPMstrong_sens_esp$ppv, na.rm = TRUE),
                                                            npv = mean(SPMstrong_sens_esp$npv, na.rm = TRUE))
    sds <- data.frame(region = region[i], group = "SD", sens = sd(SPMstrong_sens_esp$sens, na.rm = TRUE), 
                                                        esp = sd(SPMstrong_sens_esp$esp, na.rm = TRUE),
                                                        ppv = sd(SPMstrong_sens_esp$ppv, na.rm = TRUE),
                                                        npv = sd(SPMstrong_sens_esp$npv, na.rm = TRUE))
    SPMstrong_sens_esp <- rbind(SPMstrong_sens_esp, means, sds)
  
    setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/results", "/ROI", roi[k]))
    
    write_csv(SPMstrong_sens_esp, paste0("sens_esp_SPMstrong_", region[i], "_", roi[k], ".csv"), na = "NA", append = FALSE)
  
  }
  
  setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z), "/results", "/ROI", roi[k]))
  
  sens_esp_strongSPM_w32 <- read_csv(paste0("sens_esp_SPMstrong_w32_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_w79 <- read_csv(paste0("sens_esp_SPMstrong_w79_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_w214 <- read_csv(paste0("sens_esp_SPMstrong_w214_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_w271 <- read_csv(paste0("sens_esp_SPMstrong_w271_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_w413 <- read_csv(paste0("sens_esp_SPMstrong_w413_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_wroiAD <- read_csv(paste0("sens_esp_SPMstrong_wroiAD_", roi[k], ".csv"), na = "NA")
  
  
  #* Create a Complete List: ----
  
  listSCC <- ls(pattern = "^sens_esp_SCC")
  
  for (i in 1:length(listSCC)) {
    
    data <- get(listSCC[[i]])[1:25, ]
    method <- rep("SCC", times = 25)
    Roi <- rep(as.character(roi[k]), times = 25) 
    tempSCC <- cbind(as.data.frame(method), data[, 1], as.data.frame(Roi), data[, 2:6])  
    
    SCC_vs_SPM_complete <- rbind(SCC_vs_SPM_complete, tempSCC)
  
    }
  
  
  listSPM <- ls(pattern = "^sens_esp_SPM")
  
  for (i in 1:length(listSPM)) {
    
    data <- get(listSPM[[i]])[1:25, ]
    method <- rep("SPM", times = 25)
    Roi <- rep(as.character(roi[k]), times = 25)  
    tempSPM <- cbind(as.data.frame(method), data[, 1], as.data.frame(Roi), data[, 2:6])  
    
    SCC_vs_SPM_complete <- rbind(SCC_vs_SPM_complete, tempSPM)
  
  }
  
  
  listSPMstrong <- ls(pattern = "^sens_esp_strongSPM")
  
  for (i in 1:length(listSPMstrong)) {
    
    data <- get(listSPMstrong[[i]])[1:25, ]
    method <- rep("SPMstrong", times = 25)
    Roi <- rep(as.character(roi[k]), times = 25)  
    tempSPMstrong <- cbind(as.data.frame(method), data[, 1], as.data.frame(Roi), data[, 2:6])  
    
    SCC_vs_SPM_complete <- rbind(SCC_vs_SPM_complete, tempSPMstrong)
  
  }
  
  
  #* Create a Final Reduced List: ----
  
  
  for (i in 1:length(listSCC)) {
    
    data <- get(listSCC[[i]])[26:27, ]
    temp <- data.frame(method = "SCC",
                       region = as.character(data$region[1]), 
                       roi = roi[k], 
                       sensMEAN = data$sens[1], 
                       sensSD = data$sens[2], 
                       espMEAN = data$esp[1], 
                       espSD = data$esp[2],
                       ppvMEAN = data$ppv[1],
                       ppvSD = data$ppv[2],
                       npvMEAN = data$npv[1],
                       npvSD = data$npv[2])    
    
    SCC_vs_SPM <- rbind(SCC_vs_SPM, temp)
  
    }
  
  
  for (i in 1:length(listSPM)) {
    
    data <- get(listSPM[[i]])[26:27, ]
    temp <- data.frame(method = "SPM",
                       region = as.character(data$region[1]), 
                       roi = roi[k], 
                       sensMEAN = data$sens[1], 
                       sensSD = data$sens[2], 
                       espMEAN = data$esp[1], 
                       espSD = data$esp[2],
                       ppvMEAN = data$ppv[1],
                       ppvSD = data$ppv[2],
                       npvMEAN = data$npv[1],
                       npvSD = data$npv[2]) 
    
    SCC_vs_SPM <- rbind(SCC_vs_SPM, temp)
  
    }
  
  
  for (i in 1:length(listSPMstrong)) {
    
    data <- get(listSPMstrong[[i]])[26:27, ]
    temp <- data.frame(method = "SPMstrong",
                       region = as.character(data$region[1]), 
                       roi = roi[k], 
                       sensMEAN = data$sens[1], 
                       sensSD = data$sens[2], 
                       espMEAN = data$esp[1], 
                       espSD = data$esp[2],
                       ppvMEAN = data$ppv[1],
                       ppvSD = data$ppv[2],
                       npvMEAN = data$npv[1],
                       npvSD = data$npv[2])    
    
    SCC_vs_SPM <- rbind(SCC_vs_SPM, temp)
  
    }

}

View(SCC_vs_SPM)
View(SCC_vs_SPM_complete)

setwd(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z), "/results"))
saveRDS(SCC_vs_SPM, file = paste0("SCC_vs_SPM", ".RDS"))
saveRDS(SCC_vs_SPM_complete, file = paste0("SCC_vs_SPM_complete", ".RDS"))

# NOW WE ALREADY HAVE THE TWO FILES WE NEED WITH SENSIBILITY, SPECIFICITY, 
# POSITIVE PREDICTIVE VALUE, AND NEGATIVE PREDICTIVE VALUE. LET'S MOVE TO ANOTHER
# SCRIPT FOR VISUALIZATION