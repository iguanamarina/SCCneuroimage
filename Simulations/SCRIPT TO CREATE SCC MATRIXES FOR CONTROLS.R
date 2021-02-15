#################################################
#  SCRIPT TO CREATE SCC MATRIXES FOR CONTROLS
#################################################

library(gamair);library(oro.nifti);library(memisc);library(devtools);library(remotes);library(readr);library(imager);library(itsadug);library(ggplot2);library(contoureR);library(fields);library(BPST);library(Triangulation);library(ImageSCC)


### ################################################## ###
#####           *FUNCIÓN CARGA DE DATOS*              ####
### ################################################## ###


setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/masked")

number <- paste0("C",1:25)
name <- paste0("masked_swww",number,"_tripleNormEsp_w00_rrec_OSEM3D_32_it1")

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

test <- f.clean(name[1])
View(test)


### ################################################## ###
#####             *CREATE   DATABASE*                 ####
### ################################################## ###


database <- data.frame(CN_number=integer(),z=integer(),x=integer(),y=integer(),pet=integer())
#create data.frame to include data

for (i in 1:length(name)) {
  
  temporal <- f.clean(name[i])
  CN_number <- rep(number[i], length.out = 9919)
  temporal <- cbind(CN_number,temporal)
  database <- rbind(database,temporal)
}

nrow(database[database$pet<0,]) # No negative values


### ################################################## ###
#####                *SCC   MATRIX*                   ####
### ################################################## ###

SCC_CN <- matrix(nrow = 25, ncol = 9919)

for (i in 1:length(number)) {

  Y <- subset(database, database$CN_number == number[i] & database$z == 35) 
  Y <- Y[1:9919,5] 
  Y <- as.matrix(Y)
  Y = t(Y) 
  Y[is.nan(Y)] <- 0
  SCC_CN[i, ] <- Y
  
}

setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix")
save(SCC_CN, file = "SCC_CN.RData")
