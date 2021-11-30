#### ########################################################## ####
#### ########################################################## ####
######                * MASTER SCRIPT FOR Z'S*                 #####
#### ########################################################## ####
#### ########################################################## ####

library(gamair);library(oro.nifti);library(memisc);library(devtools);library(remotes);library(readr);library(imager);library(itsadug);library(ggplot2);library(contoureR);library(fields);library(BPST);library(Triangulation);library(ImageSCC)

setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS")

# HIPERPARAMETER:

param.z = 30

### ########################################################## ###
#####                  * LOAD    TEMPLATES *                  #### 
### ########################################################## ###

f.clean <- function(name) {
  
  ## Load Data
  
  file <- readNIfTI(fname = name, verbose = FALSE, warn = -1, reorient = TRUE, call = NULL, read_data = TRUE)
  namex <- as.character(name)
  n = img_data(file)
  n = to.data.frame(n)

  dataframe <- data.frame(z = integer(),x = integer(),y = integer(),pet = integer()) 
  
  # Loop for 91 slices of Z in the NiFtI image -> move to dataframe
  
  for (i in seq(1:91)) {
    
    n_lim = n[n$Var2==i,] # Select just one Z slice
    n_lim$Var1=NULL
    n_lim$Var2=NULL
    
    z <- rep(i, length.out = 9919)
    x <- rep(1:91, each = 109, length.out = 9919) 
    y <- rep(1:109, length.out = 9919)
    
    attach(n_lim)
    pet <- c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`,`13`,`14`,
             `15`,`16`,`17`,`18`,`19`,`20`,`21`,`22`,`23`,`24`,`25`,`26`,
             `27`,`28`,`29`,`30`,`31`,`32`,`33`,`34`,`35`,`36`,`37`,`38`,
             `39`,`40`,`41`,`42`,`43`,`44`,`45`,`46`,`47`,`48`,`49`,`50`,
             `51`,`52`,`53`,`54`,`55`,`56`,`57`,`58`,`59`,`60`,`61`,`62`,
             `63`,`64`,`65`,`66`,`67`,`68`,`69`,`70`,`71`,`72`,`73`,`74`,
             `75`,`76`,`77`,`78`,`79`,`80`,`81`,`82`,`83`,`84`,`85`,`86`,
             `87`,`88`,`89`,`90`,`91`)
    detach(n_lim)
    
    temp0 = data.frame(z,x,y,pet) # temporal dataframe
    temp1 <- print(temp0) 
    dataframe <- rbind(dataframe,temp1)
  }

  print(dataframe) # Necessary for assigning an object name

}

Data = f.clean("new_mask") 

### ########################################################## ###
#####                *CONTOURS OF NEURO-DATA*                 ####
### ########################################################## ###

Y <- subset(Data, Data$z == param.z) 
Y <- Y[1:9919,4] 
Y <- as.matrix(Y)
Y = t(Y) 
Y[is.nan(Y)] <- 0
SCC <- Y; rm(Y)

# Z are the coordinates where data is measured:

x <- rep(1:91, each = 109, length.out = 9919) 
y <- rep(1:109,length.out = 9919)
Z <- cbind(as.matrix(x),as.matrix(y)); Z

dat <- cbind(Z,t(SCC))

dat <- as.data.frame(dat)
dat[is.na(dat)] <- 0
sum(is.na(dat$pet)) # should be = 0
head(dat)

rownames(dat) <- NULL 

df = getContourLines(dat[1:9919,], levels = c(0)) 

ggplot(df, aes(x, y, colour = z)) + geom_path() 

contour = df 
head(contour);str(contour)

f.contour <- function(x){
  
  aa <- contour[contour$GID == x,] # We keep GID==x (0,1,2,3...)
  a <- aa[,5:6] # Then we keep just the coordinates 
  print(a) # and then print in order to loop and make a list
}

coord <- list()

for (i in 0:max(contour$GID)) { #change contour30 to any other name previously assigned if necessary
   
  coord[[i + 1]] <- f.contour(i)
  rownames(coord[[ i + 1 ]]) <- NULL
}


### Test the results are coherent:

head(coord[[1]],10); plot(coord[[1]])   # external boundaries
head(coord[[2]],10); points(coord[[2]]) # first hole
head(coord[[3]],10); points(coord[[3]]) # second hole
# (...) 


### ########################################################## ###
#####             *TRIANGULATION PARAMETERS*                  ####
### ########################################################## ###

VT = TriMesh(coord[[1]], 8) # Triangulation degree of fineness: tuning parameter

head(VT$V,10);head(VT$Tr,10) 

# Create it first if it doesn't exist
setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.character(param.z)))
save(VT, file = paste0("contour", as.character(param.z), ".RData"))
     


#################################################
#  THIS PART IS TO CREATE SCC MATRIXES FOR CONTROLS
#################################################


### ################################################## ###
#####           *FUNCIÓN CARGA DE DATOS*              ####
### ################################################## ###

setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/masked")

number <- paste0("C", 1:25)
name <- paste0("masked_swww", number, "_tripleNormEsp_w00_rrec_OSEM3D_32_it1")


### ################################################## ###
#####             *CREATE   DATABASE*                ##### 
### ################################################## ###

database_CN <- data.frame(CN_number = integer(),z = integer(), x = integer(), y = integer(), pet = integer())

for (i in 1:length(name)) {
  
  temporal <- f.clean(name[i])
  CN_number <- rep(number[i], length.out = 9919)
  temporal <- cbind(CN_number,temporal)
  database_CN <- rbind(database_CN,temporal)
}

nrow(database_CN[database_CN$pet < 0, ]) # No negative values
rm(temporal)


### ################################################## ###
#####                *SCC   MATRIX*                  #####
### ################################################## ###

SCC_CN <- matrix(nrow = 25, ncol = 9919)

for (i in 1:length(number)) {

  Y <- subset(database_CN, database_CN$CN_number == number[i] & database_CN$z == param.z) 
  Y <- Y[1:9919,5] 
  Y <- as.matrix(Y)
  Y = t(Y) 
  Y[is.nan(Y)] <- 0
  SCC_CN[i, ] <- Y
  
}


# Sometimes its necessary:

  # na.zero <- function(x) {
  #     x[is.na(x)] <- 0
  # }
  # 
  # 
  # SCC_CN <-  apply(SCC_CN, 2, na.zero)


setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.character(param.z)))
save(SCC_CN, file = "SCC_CN.RData")


##########################################################################
#  SCRIPT TO CREATE SCC MATRIXES FOR PATHOLOGICAL 
#  THIS SCRIPT PRODUCES SCC MATRICES PRESENT IN SIMULACIONES/SCC_matrix
#  AFTER THIS SCRIPT GO TO  ->  "SCRIPT FOR SCCs.R" 
##########################################################################



### ################################################## ###
#####           *FUNCIÓN CARGA DE DATOS*              ####
### ################################################## ###

setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/masked")

number <- paste0("C", 1:25)
region <- c("roiAD", "w32", "w79", "w214", "w271", "w413")
roi <- c(1, 2, 4, 6, 8)

for (i in 1:length(region)) {
  
  for (j in 1:length(roi)) {
    
    setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/masked")
    
    # IMPORTANTE: tuning parameters
    
    REGION = region[i]
    ROI = roi[j] 
    name <- paste0("masked_swww", number, "_tripleNormEsp_", REGION, "_0_", ROI, "_rrec_OSEM3D_32_it1")
    
    
    ### ################################################## ###
    #####             *CREATE   DATABASE*                ##### This can be reduced with functions (!!)
    ### ################################################## ###
    
    database_AD <- data.frame(AD_number = integer(),z = integer(), x = integer(), y = integer(), pet = integer())
    
    for (k in 1:length(name)) {
      
      temporal <- f.clean(name[k])
      AD_number <- rep(paste0(number[k], "_", as.character(REGION), "_", as.character(ROI)), length.out = 9919)
      temporal <- cbind(AD_number, temporal)
      database_AD <- rbind(database_AD, temporal)
      
    }
    
    
    SCC_matrix <- matrix(nrow = length(name), ncol = 9919)
    
    for (k in 1:length(number)) {
    
      Y <- subset(database_AD, AD_number == paste0(number[k], "_", as.character(REGION), "_", as.character(ROI)) &
                            z == param.z)
      Y <- Y[1:9919, 5] 
      Y <- as.matrix(Y)
      Y = t(Y) 
      Y[is.nan(Y)] <- 0
      SCC_matrix[k, ] <- Y
      
    }
    
    setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.character(param.z)))
    save(SCC_matrix, file = paste0("SCC_", REGION, "_", ROI, ".RData"))

  }
}


### ########################################################## ###
#####            **   SCRIPT YA PARA SCCs  **                 ####
### ########################################################## ###

# install.packages("remotes")
# remotes::install_github("skgrange/threadr")

library(gamair);library(oro.nifti);library(remotes);library(readr);library(imager);library(itsadug);library(ggplot2);library(contoureR);library(fields);library(BPST);library(Triangulation);library(ImageSCC); library(tidyr); library(dplyr);library(stringr); library(threadr);library(memisc)

setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.character(param.z)))

# In order to be consistent we use common names Brain.V and Brain.Tr. From here onwards most of the names follow the ones provided by Wang et al (2019)

Brain.V <- VT[[1]]
Brain.Tr <- VT[[2]]

head(Brain.V);head(Brain.Tr)


V.est = as.matrix(Brain.V)
# Brain.v<-cbind(Brain.V[,2],Brain.V[,1]) # In case you need to transpose the data
Tr.est = as.matrix(Brain.Tr)

V.band = as.matrix(Brain.V)
Tr.band = as.matrix(Brain.Tr) 


region <- c("roiAD", "w32", "w79", "w214", "w271", "w413")
roi <- c(1, 2, 4, 6, 8) 
number <- paste0("C", 1:25)


### HERE STARTS THE LOOP: LADIES AND GENTLEMAN...FASTEN YOUR SEATBELTS ###

## TIEMPOS -> START AT:

for (i in 1:length(region)) {
  
  for (j in 1:length(roi)) {
    
    setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.character(param.z)))
    
    # Response Variable:
    
    name_CN <- paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.character(param.z), "/SCC_CN.RData")
    name_AD <- paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.character(param.z), "/SCC_", 
                      region[i], "_", roi[j], ".RData")

    SCC_CN <- threadr::read_rdata(name_CN)
    SCC_AD <- threadr::read_rdata(name_AD)
    
      # Mean Average Normalization: TEST WHETHER THIS HELPS (!!)
      
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
    
    
    ### ########################################################## ###
    #####        *OTHER PARAMETERS FOR SCC ESTIMATION*            ####
    ### ########################################################## ###
    
    # Following Wang et al's recomendations:
    
    d.est = 5 # degree of spline for mean function  5
    d.band = 2 # degree of spline for SCC  2
    r = 1 # smoothing parameter  1
    lambda = 10^{seq(-6,3,0.5)} # penalty parameters
    alpha.grid = c(0.10,0.05,0.01) # vector of confidence levels
    x <- rep(1:91, each = 109, length.out = 9919) 
    y <- rep(1:109,length.out = 9919)
    Z <- cbind(as.matrix(x),as.matrix(y)); Z
    
    ### ########################################################## ###
    #####               *CONSTRUCTION OF SCC'S*                   ####
    ### ########################################################## ###
    
    setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.character(param.z), "/SCC", param.z, "results"))
    
    SCC_COMP = scc.image(  Ya = SCC_AD, Yb = SCC_CN, Z = Z, d.est = d.est, d.band = d.band, r = r,
                           V.est.a = V.est, Tr.est.a = Tr.est,
                           V.band.a = V.band, Tr.band.a = Tr.band,
                           penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                           adjust.sigma = TRUE)    
    
    save(SCC_COMP, file = paste0("SCC_COMP_", region[i], "_", roi[j] ,".RData"))
    
  }
  
}


# This function will be necessary:

    
my_points <- function(aa){
  
  Z.band <- matrix(aa$Z.band, ncol = 2) # Positions
  z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Separated positions
  n1 <- length(z1); n2 <- length(z2) # Lengths of those positions
  scc <- matrix(NA,n1*n2,2) # Matrix with value of that SCC in that position
  ind.inside.band <- aa$ind.inside.cover # Keep only regions inside triangulation
  scc[ind.inside.band,] <- aa$scc[,,2] # Assign SCC to those areas (IMPORTANT_: 1,2,3 FOR DIFFERENT ALPHAS)
  scc.limit <- c(min(scc[,1],na.rm = TRUE), max(scc[,2],na.rm = TRUE)) # LIMITS: minimum of inferior, maximum of superior
   
  scc.l.mtx <- matrix(scc[,1], nrow = n2, ncol = n1) # Lower SCC for each location. 
  scc.u.mtx <- matrix(scc[,2], nrow = n2, ncol = n1) # Upper SCC for each location.
  scc.l.mtx[scc.l.mtx < 0] = NA 
  # The ones that work just fine are substituded by a NA as we don't want to represent them, only positive LowerSCCs are represented
  scc.u.mtx[scc.u.mtx > 0] = NA 
  # The ones that work just fine are substituded by a NA as we don't want to represent them, only positive UpperSCCs are represented
  
  points.P <- which(scc.l.mtx > 0, arr.ind = TRUE) # Points where mean difference is positive (first image is stronger)
  points.N <- which(scc.u.mtx < 0, arr.ind = TRUE) # Points where mean difference is negative (second image is stronger)
  
  pointers <- list(points.P, points.N)
  print(pointers)  
}


### ########################################################## ###
#####                 SCRIPT  FOR  EVALUATION                 ####
### ########################################################## ###

library(gamair);library(oro.nifti);library(memisc);library(devtools);library(remotes);library(readr);library(imager);library(itsadug);library(ggplot2);library(contoureR);library(fields);library(BPST);library(Triangulation);library(ImageSCC); library(tidyr); library(dplyr);library(stringr)


### ########################################################## ###
#####             THESE ARE THE THEORETICAL ROIs              ####
### ########################################################## ###


# TUNING PARAMETERS - Combinatory of regions and induced hypoactivities:

region <- c("w32", "w79", "w214", "w271", "w413", "wroiAD")
number <- paste0("C", 1:25)


# DOUBLE LOOP REGION-ROI: with this we get a series of tables with the exact ROI points as provided by CHUS.

for (i in 1:length(region)) { # RUN JUST THE FIRST TIME

    for (k in 1:length(number)) {

      setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/roisNormalizadas")
      
      roi_table <- data.frame(group = integer(),
                              z = integer(),
                              x = integer(),
                              y = integer())
      
      # Name of file:
      name <- paste0("wwwx", region[i], "_redim_crop_squ_flipLR_newDim_", number[k])
      # Clean data for that file:
      temporal <- f.clean(name)
      # Add meta-data:
      group <- rep(paste0(as.character(region[i]), "_number", number[k]), length.out = 9919)
      # Merge with main data.frame:
      temporal <- cbind(group, temporal)
      # Solve problem with zeros and NA's:
      suelto <- temporal[, "pet"]
      suelto[is.nan(suelto)] <- 0
      temporal[, "pet"] <- as.data.frame(suelto)
      # Final:
      roi_table <- rbind(roi_table, temporal)

      setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/roisNormalizadas/tables")
      saveRDS(roi_table, file = paste0("ROItable_", region[i], "_", number[k], ".RDS"))
  
    }
  
  }


  # Theoretical Points (True according to ROIs):
  

ROI_data <- list() # Cargar las databases

  for (i in 1:length(region)) {
    for (k in 1:length(number)) {
      
    setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/roisNormalizadas/tables")
    ROI_data[[paste0(as.character(region[i]), "_", as.character(number[k]))]] <- readRDS(
      paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/roisNormalizadas/tables/ROItable_", 
             region[i], "_", number[k], ".RDS"))
  
    }
  }  


T_points <- list() # Theoretical points filtered by Z 

for (i in 1:length(region)) {
  
  for (j in 1:length(number)) {
    
    T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]] <- ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]][ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]]$z == param.z & ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]]$pet == 1, 3:4]
    
    T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]] <- unite(as.data.frame(T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]]), newcol, c(y,x), remove = T)
        
    }
  }
    
rm(ROI_data) # No longer necessary


### ############################################################################# ###
#####             THESE ARE THE HYPOTHETICAL POINTS DETECTED BY SCCs             #### 
### ############################################################################# ###

roi <- c(1, 2, 4, 6, 8)

SCC_vs_SPM_complete <-   data.frame(method = integer(),
                                    region = integer(),
                                    roi = integer(),
                                    number = integer(),
                                    sens = integer(),
                                    esp = integer())

SCC_vs_SPM <- data.frame(method = integer(),
                         region = integer(),
                         roi = integer(),
                         sensMEAN = integer(),
                         sensSD = integer(),
                         espMEAN = integer(),
                         espSD = integer())

x <- rep(1:91, each = 109, length.out = 9919) 
y <- rep(1:109, length.out = 9919)
total.coords <- data.frame(y, x) ###!!!!!
total.coords <- unite(as.data.frame(total.coords), newcol, c(y, x), remove = T)

for (k in 1:length(roi)) {
  
  # Hypotetical Points (according to SCCs): 
  
  region <- c("w32", "w79", "w214", "w271", "w413", "roiAD")
  
  H_points <- list()
  
  setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.numeric(param.z) ,"/SCC", as.numeric(param.z), "results"))
  
  for (i in 1:length(region)) {
  
    load(paste0("SCC_COMP_", region[i], "_", roi[k], ".RData"))
    H_points[[as.character(region[i])]] <- my_points(SCC_COMP)[[1]] 
    H_points[[as.character(region[i])]] <- unite(as.data.frame(H_points[[i]]), newcol, c(row, col), remove = T)
    
  }
     
  rm(list = ls(pattern = "^SCC_COMP")) # Just beautiful
  
  
  ### ########################################################## ###
  #####                   SENS/ESP   FOR   SCC                  ####
  ### ########################################################## ###
  
  region <- c("w32", "w79", "w214", "w271", "w413", "wroiAD")
  names(H_points)[6] <- "wroiAD" # because of reasons
  
  for (i in 1:length(region)) {  
    
    SCC_sens_esp <- data.frame(region = integer(), group = integer(), sens = integer(), esp = integer())
    
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
    
      temp <- data.frame(region = region[i], group = number[j], sens = sensibilitySCC, esp = specificitySCC)  
      SCC_sens_esp <- rbind(SCC_sens_esp, temp)
      
    }
    
    means <- data.frame(region = region[i], group = "MEAN", sens = mean(SCC_sens_esp$sens, na.rm = TRUE), 
                                                            esp = mean(SCC_sens_esp$esp, na.rm = TRUE))
    sds <- data.frame(region = region[i], group = "SD", sens = sd(SCC_sens_esp$sens, na.rm = TRUE), 
                                                        esp = sd(SCC_sens_esp$esp, na.rm = TRUE))
    SCC_sens_esp <- rbind(SCC_sens_esp, means, sds)
    
    setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z), "/ROI", roi[k]))
    
    write_csv(SCC_sens_esp, paste0("sens_esp_SCC_", region[i], "_", roi[k], ".csv"), na = "NA", append = FALSE)
    
  }
  
  setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z), "/ROI", roi[k]))
  
  sens_esp_SCC_w32 <- read_csv(paste0("sens_esp_SCC_w32_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_w79 <- read_csv(paste0("sens_esp_SCC_w79_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_w214 <- read_csv(paste0("sens_esp_SCC_w214_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_w271 <- read_csv(paste0("sens_esp_SCC_w271_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_w413 <- read_csv(paste0("sens_esp_SCC_w413_", roi[k], ".csv"), na = "NA")
  sens_esp_SCC_wroiAD <- read_csv(paste0("sens_esp_SCC_wroiAD_", roi[k], ".csv"), na = "NA")
  
  
  ### ########################################################## ###
  #####                  SENS/ESP  FOR  SPM                     #### 
  ### ########################################################## ###
  
  
  for (i in 1:length(region)) {  
    
    setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SPM/z", as.numeric(param.z), "/ROI", i, "_", region[i], "_0", roi[k]))
    binary <- f.clean("binary.nii")
    H_points_SPM <- binary[binary$z == as.numeric(param.z) & binary$pet == 1, 2:3]
    H_points_SPM <- unite(as.data.frame(H_points_SPM), newcol, c(y, x), remove = T)
    
    SPM_sens_esp <- data.frame(region = integer(), group = integer(), sens = integer(), esp = integer())
    
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
      
      temp <- data.frame(region = region[i], group = number[j], sens = sensibilitySPM, esp = specificitySPM)  
      SPM_sens_esp <- rbind(SPM_sens_esp, temp)
  
    }
    
    means <- data.frame(region = region[i], group = "MEAN", sens = mean(SPM_sens_esp$sens, na.rm = TRUE), 
                                                            esp = mean(SPM_sens_esp$esp, na.rm = TRUE))
    sds <- data.frame(region = region[i], group = "SD", sens = sd(SPM_sens_esp$sens, na.rm = TRUE), 
                                                        esp = sd(SPM_sens_esp$esp, na.rm = TRUE))
    SPM_sens_esp <- rbind(SPM_sens_esp, means, sds)
  
    setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z), "/ROI", roi[k]))
    
    write_csv(SPM_sens_esp, paste0("sens_esp_SPM_", region[i], "_", roi[k], ".csv"), na = "NA", append = FALSE)
  
  }
  
  setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z), "/ROI", roi[k]))
  
  sens_esp_SPM_w32 <- read_csv(paste0("sens_esp_SPM_w32_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_w79 <- read_csv(paste0("sens_esp_SPM_w79_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_w214 <- read_csv(paste0("sens_esp_SPM_w214_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_w271 <- read_csv(paste0("sens_esp_SPM_w271_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_w413 <- read_csv(paste0("sens_esp_SPM_w413_", roi[k], ".csv"), na = "NA")
  sens_esp_SPM_wroiAD <- read_csv(paste0("sens_esp_SPM_wroiAD_", roi[k], ".csv"), na = "NA")
  
  
  ### ########################################################## ###
  #####                SENS/ESP  FOR  STRONG SPM                ####     FWE 0'05 + 100 voxels threshold
  ### ########################################################## ###
  
  
  for (i in 1:length(region)) {  
    
    setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SPMstrong/z", as.numeric(param.z), "/ROI", i, "_", region[i], "_0", roi[k]))
    binary <- f.clean("binary.nii")
    H_points_SPMstrong <- binary[binary$z == as.numeric(param.z) & binary$pet == 1, 2:3]
    H_points_SPMstrong <- unite(as.data.frame(H_points_SPMstrong), newcol, c(y, x), remove = T)
    
    SPMstrong_sens_esp <- data.frame(region = integer(), group = integer(), sens = integer(), esp = integer())
    
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
      
      temp <- data.frame(region = region[i], group = number[j], sens = sensibilitySPMstrong, esp = specificitySPMstrong)  
      SPMstrong_sens_esp <- rbind(SPMstrong_sens_esp, temp)
  
    }
    
    means <- data.frame(region = region[i], group = "MEAN", sens = mean(SPMstrong_sens_esp$sens, na.rm = TRUE), 
                                                            esp = mean(SPMstrong_sens_esp$esp, na.rm = TRUE))
    sds <- data.frame(region = region[i], group = "SD", sens = sd(SPMstrong_sens_esp$sens, na.rm = TRUE), 
                                                        esp = sd(SPMstrong_sens_esp$esp, na.rm = TRUE))
    SPMstrong_sens_esp <- rbind(SPMstrong_sens_esp, means, sds)
  
    setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z), "/ROI", roi[k]))
    
    write_csv(SPMstrong_sens_esp, paste0("sens_esp_SPMstrong_", region[i], "_", roi[k], ".csv"), na = "NA", append = FALSE)
  
  }
  
  setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z), "/ROI", roi[k]))
  
  sens_esp_strongSPM_w32 <- read_csv(paste0("sens_esp_SPMstrong_w32_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_w79 <- read_csv(paste0("sens_esp_SPMstrong_w79_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_w214 <- read_csv(paste0("sens_esp_SPMstrong_w214_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_w271 <- read_csv(paste0("sens_esp_SPMstrong_w271_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_w413 <- read_csv(paste0("sens_esp_SPMstrong_w413_", roi[k], ".csv"), na = "NA")
  sens_esp_strongSPM_wroiAD <- read_csv(paste0("sens_esp_SPMstrong_wroiAD_", roi[k], ".csv"), na = "NA")
  
  
  ## Create a Complete List: 
  
  listSCC <- ls(pattern = "^sens_esp_SCC")
  
  for (i in 1:length(listSCC)) {
    
    data <- get(listSCC[[i]])[1:25, ]
    method <- rep("SCC", times = 25)
    Roi <- rep(as.character(roi[k]), times = 25) 
    tempSCC <- cbind(as.data.frame(method), data[, 1], as.data.frame(Roi), data[, 2:4])  
    
    SCC_vs_SPM_complete <- rbind(SCC_vs_SPM_complete, tempSCC)
  
    }
  
  
  listSPM <- ls(pattern = "^sens_esp_SPM")
  
  for (i in 1:length(listSPM)) {
    
    data <- get(listSPM[[i]])[1:25, ]
    method <- rep("SPM", times = 25)
    Roi <- rep(as.character(roi[k]), times = 25)  
    tempSPM <- cbind(as.data.frame(method), data[, 1], as.data.frame(Roi), data[, 2:4])  
    
    SCC_vs_SPM_complete <- rbind(SCC_vs_SPM_complete, tempSPM)
  
  }
  
  
  listSPMstrong <- ls(pattern = "^sens_esp_strongSPM")
  
  for (i in 1:length(listSPMstrong)) {
    
    data <- get(listSPMstrong[[i]])[1:25, ]
    method <- rep("SPMstrong", times = 25)
    Roi <- rep(as.character(roi[k]), times = 25)  
    tempSPMstrong <- cbind(as.data.frame(method), data[, 1], as.data.frame(Roi), data[, 2:4])  
    
    SCC_vs_SPM_complete <- rbind(SCC_vs_SPM_complete, tempSPMstrong)
  
  }
  
  
  ## Create a Final Reduced List:
  
  
  for (i in 1:length(listSCC)) {
    
    data <- get(listSCC[[i]])[26:27, ]
    temp <- data.frame(method = "SCC",
                       region = as.character(data$region[1]), 
                       roi = roi[k], 
                       sensMEAN = data$sens[1], 
                       sensSD = data$sens[2], 
                       espMEAN = data$esp[1], 
                       espSD = data$esp[2])    
    
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
                       espSD = data$esp[2])    
    
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
                       espSD = data$esp[2])    
    
    SCC_vs_SPM <- rbind(SCC_vs_SPM, temp)
  
    }

}

View(SCC_vs_SPM)
View(SCC_vs_SPM_complete)

setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Preliminary results/z", as.numeric(param.z)))
saveRDS(SCC_vs_SPM, file = paste0("SCC_vs_SPM", ".RDS"))
saveRDS(SCC_vs_SPM_complete, file = paste0("SCC_vs_SPM_complete", ".RDS"))


### ########################################################## ###
#####               SCRIPT  FOR  VISUALIZATION                ####
### ########################################################## ###

require(ggplot2)

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
                 xlab("Level of Induced Hypoactivity (%)") + 
                 ylab("Sensibility (%)") + ggtitle(paste0("Comparison of sensibility by regions in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Methods")) + 
                 facet_wrap( ~ region) 

graph2 


graph22 <- ggplot(data = table, 
                 aes(x = Roi, y = esp)) + 
                 geom_boxplot(outlier.colour = NULL, aes(fill = method), outlier.size = 1) +
                 xlab("Level of Induced Hypoactivity (%)") + 
                 ylab("Especificity (%)") + ggtitle(paste0("Comparison of especificity by regions in z=", as.numeric(param.z))) +
                 guides(fill = guide_legend(title = "Methods")) + 
                 facet_wrap( ~ region) 

graph22 


########################################################################################
########################################################################################

# Academic Publication Theme:

install.packages("envalysis")
library(envalysis)
graph2 + theme_publish() 

ggsave(
  "academic_plot.png",
  plot = last_plot(),
  scale = 2)



