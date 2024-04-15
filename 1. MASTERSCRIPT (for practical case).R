############################## ################### ############## ### 
##
## Script name: MASTERSCRIPT (for practical case)
##
## Purpose of script: Script to import NIFTI files to R, create a database, get PPT names,
## create triangulation parameters, and carry on with SCCs calculation. This script asumes
## that Matlab NIFTI pre-processing has already been carried out correctly.
##
## Date Created: 2022-01-10
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://juan-arias.xyz
##   
############################## ################### ############## ### 


### ################################################## ###
#####                 *PREAMBLE*                      ####
####         Installation of necessary packgs.         ###
### ################################################## ###

# List of required packages
packages <- c("mgcv", "gamair", "oro.nifti", "memisc")

# Check and install each package if it is not installed, then load it
for (packageName in packages) {
  if (!packageName %in% rownames(installed.packages())) {
    install.packages(packageName)
  }
  library(packageName, character.only = TRUE)
}

# Clean up temporary objects from the environment
rm(packages)


### ################################################## ###
#####            *IMPORT NIFTI FILES*                 ####
####            f.clean  (PPT,z,x,y,pet)               ###
### ################################################## ###

# Working directory path needs to be set according to user
warning("Ensure the working directory path is set correctly for your system. This path needs to be the directory containing your .hdr/.img NIfTI files, mask, and Demographics .csv file. If you forked this project into your GitHub Folder, then you shouldn't have to change anything.")

# Modify this path according to your local environment
setwd("~/GitHub/SCCneuroimage/PETimg_masked for practical case")

# Read demographics data
demo <- read.csv2("Demographics.csv")

# Load function 'NeuroSCC::neuroCleaner' to clean and load data from NIFTI files
# This function does: read NIFTI image, transform it to dataframe, preserve slice Z and organize the table.

  neuroCleaner <- function(name, demo) {
    # Load data into a dataframe
    file <- oro.nifti::readNIfTI(fname = name, verbose = FALSE, warn = -1, reorient = TRUE, call = NULL, read_data = TRUE)
    n <- memisc::to.data.frame(img_data(file))
  
    # Get File Name
    namex <- as.character(name)
  
    # Get Dimensions of File
    xDim <- file@dim_[2]
    yDim <- file@dim_[3]
    dim <- xDim * yDim
  
    # Prepare base data.frame where data from the loop will be integrated
    dataframe <- data.frame(z = integer(), x = integer(), y = integer(), pet = integer())
  
      # Loop for every slice in that Z; then attach to dataframe
      for (i in seq(1:xDim)) {
        n_lim <- n[n$Var2 == i, ] # Select just one Z slice
        n_lim$Var1 <- NULL
        n_lim$Var2 <- NULL
    
        z <- rep(i, length.out = dim)
        x <- rep(1:xDim, each = yDim, length.out = dim)
        y <- rep(1:yDim, length.out = dim)
    
        pet <- unlist(n_lim) # Convert n_lim to a vector and avoid using attach
    
        temp <- data.frame(z, x, y, pet) # Temporal dataframe
        dataframe <- rbind(dataframe, temp) # Sum new data with previous data
      }
  
    # Check if demographic data contains required columns
    requiredColumns <- c("PPT", "Group", "Sex", "Age")
    if (!all(requiredColumns %in% names(demo))) {
      stop("Demographic data must contain columns: PPT, Group, Sex, and Age")
    }
  
    # Extract demographic data
    demog <- demo[demo$PPT == namex, ]
    if (nrow(demog) == 0) {
      stop("No demographic data found for the given participant.")
    }
  
      # Replicate demographic data for each pixel
      PPT <- rep(demog$PPT, length.out = dim)
      group <- rep(demog$Group, length.out = dim)
      sex <- rep(demog$Sex, length.out = dim)
      age <- rep(demog$Age, length.out = dim)
    
      temp2 <- data.frame(PPT, group, sex, age)
      dataframe <- cbind(temp2, dataframe)
  
    return(dataframe) # Return the dataframe
  }

# Example of conversion from NIFTI to R dataframe:
example = neuroCleaner("003_S_1059", demo = demo)


### ################################################## ###
#####             *CREATE   DATABASE*                 ####
### ################################################## ###

# (!) You can skip this part until the end of MEAN AVERAGE NORMALIZATION where
# there is code to directly load 'database.RDS' from the repo.

# Get full list of files with extension .img
files <- list.files(path = getwd(), pattern = "*.img", full.names = FALSE, recursive = FALSE) 
files <- sub(pattern = "\\.img$", replacement = "", x = files) 

# Print the name of those files, check that everything is alright
print(files)

# Create database to 'rbind' results from neuroCleaner
database <- data.frame(PPT = integer(), group = integer(), sex = integer(), age = integer(), z = integer(), x = integer(), y = integer(), pet = numeric())

# Loop to include every PPT in the dataframe
for (i in seq_along(files)) {
  temporal <- neuroCleaner(files[i], demo = demo)
  database <- rbind(database, temporal)
}

# Convert negatives and zeros to NaN in the 'pet' column
database$pet[database$pet <= 0] <- NaN 


### ####################################################### ###
#####           *MEAN AVERAGE NORMALIZATION*               ####
### ####################################################### ###

# Initialize 'pet_normal'
database$pet_normal <- NA  

for (i in seq_along(files)) {
  temp <- database[database$PPT == files[i], 'pet']
  mean_value <- mean(as.numeric(temp), na.rm = TRUE)
  database$pet_normal[database$PPT == files[i]] <- temp / mean_value
}

# Remove old 'pet' column
database <- database[, -which(colnames(database) == "pet")]

# Export complete database 
# setwd("~/GitHub/SCCneuroimage")
# saveRDS(database, file = "Auxiliary Files/database.RDS") # export as .RDS IF DESIRABLE

# (!) LOAD THE DATABASE HERE IF YOU WANT TO SAVE TIME
database <- readRDS("~/GitHub/SCCneuroimage/Auxiliary Files/database.RDS")


### ########################################################## ###
#####      SIMULTANEOUS CONFIDENCE CORRIDORS (preamble)       ####
### ########################################################## ###

# This part of the code needs a series of packages which are only available on GitHub
# Besides, we will also need some extra packages not load yet:

# List of CRAN packages
packagesCRAN <- c("devtools", "remotes", "readr", "imager", "itsadug", "ggplot2", "contoureR", "fields")

# Install CRAN packages if not already installed and load them
for (pkg in packagesCRAN) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Set environment to not stop installations due to warning messages
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)

# List of GitHub packages
packagesGithub <- c("funstatpackages/BPST", "funstatpackages/Triangulation", "funstatpackages/ImageSCC")

# Install and load GitHub packages
for (pkg in packagesGithub) {
  remotes::install_github(pkg)
}

# Load all packages
library(devtools); library(remotes); library(readr); library(imager); library(itsadug); library(ggplot2); library(contoureR); library(fields); library(BPST); library(Triangulation); library(ImageSCC)

# Clean up temporary objects from the environment
rm(packagesCRAN, packagesGithub)


### ########################################################## ###
#####           *COMPREHENSIVE LISTS OF PPT'S*                ####
### ########################################################## ###

# In Functional Data Analysis (FDA), data transformation is essential. Specifically, we must format the data into 
# data.frames where each participant (PPT) is represented by a single row, and each column corresponds to one of the 7505 
# PET values from a selected brain slice (79*95=7505). This approach is necessary for FDA as it allows us to analyze one 
# brain slice at a time. Additionally, we need to predefine lists of PPTs by their sex and diagnostic group 
# categories (e.g., AD, CN, male, female) to facilitate the creation of SCC matrices tailored to these groupings.

# Create lists of names for Alzheimer's Disease (AD) and Control (CN) groups using unique participant identifiers
names_AD <- unique(database$PPT[database$group == "AD"])
names_CN <- unique(database$PPT[database$group == "CN"])

# Create lists of names for males and females using unique participant identifiers
names_M <- unique(database$PPT[database$sex == "F"])
names_F <- unique(database$PPT[database$sex == "M"])

# Create combined lists for AD and CN by sex using unique participant identifiers
names_AD_M <- unique(database$PPT[database$group == "AD" & database$sex == "M"])
names_AD_F <- unique(database$PPT[database$group == "AD" & database$sex == "F"])
names_CN_M <- unique(database$PPT[database$group == "CN" & database$sex == "M"])
names_CN_F <- unique(database$PPT[database$group == "CN" & database$sex == "F"])


### ########################################################## ###
#####                  *GENERATE MATRICES*                    ####
### ########################################################## ###


# !!! debugear aqui


# Function to obtain dimensions from a DICOM file
getDimensions <- function(filename = NULL) {
  # If no filename is provided, find the first .img file in the current working directory
  if (is.null(filename)) {
    files <- list.files(path = getwd(), pattern = "\\.img$", full.names = TRUE, recursive = FALSE)
    if (length(files) == 0) {
      stop("No .img files found in the directory.")
    }
    filename <- files[1]  # Use the first file found
  }
  
  # Load the data using oro.nifti
  file <- oro.nifti::readNIfTI(fname = filename, verbose = FALSE, warn = -1, reorient = TRUE)
  
  # Get dimensions X and Y, and calculate 'dim'
  xDim <- file@dim_[2]
  yDim <- file@dim_[3]
  dim <- xDim * yDim
  
  # Return a list with X, Y, and 'dim' dimensions
  return(list(xDim = xDim, yDim = yDim, dim = dim))
}

# If 'dim' is not predefined, obtain it using:
dimensions <- getDimensions()
dim <- dimensions$dim

# Function to create matrices based on participant identifier (PPT) and sex
matrixCreator <- function(names, z_slice = 30) {
  # Ensure 'dim' is defined. Here we assume 'dim' is obtained from 'getDimensions' and globally available
  if (!exists("dim")) {
    stop("Dimension 'dim' is not defined. Please define it using getDimensions().")
  }
  
  # Initialize the SCC matrix with appropriate dimensions
  SCC_matrix <- matrix(ncol = dim, nrow = 0)
  
  for (i in seq_along(names)) {
    # Subset data for each participant at the specified Z slice and take the first 'dim' entries
    Y <- subset(database, database$ppt == names[i] & database$z == z_slice)
    Y <- Y[1:dim, "pet_normal"]  
    Y <- as.matrix(Y)
    Y <- t(Y)
    Y[is.na(Y)] <- 0  # Replace NA and NaN values with 0
    
    # Bind each participant's data to the SCC matrix
    SCC_matrix <- rbind(SCC_matrix, Y)
  }
  
  return(SCC_matrix)
}

# Create matrices for each group and sex using the matrixCreator function

matrices <- list()
matrices$AD <- matrixCreator(names_AD)
matrices$CN <- matrixCreator(names_CN)
matrices$Male <- matrixCreator(names_M)
matrices$Female <- matrixCreator(names_F)

matrices$AD_M <- matrixCreator(names_AD_M)
matrices$AD_F <- matrixCreator(names_AD_F)
matrices$CN_M <- matrixCreator(names_CN_M)
matrices$CN_F <- matrixCreator(names_CN_F)







### ########################################################## ###
#####           *FUNCTIONS FOR SCCs DATA.FRAME*               ####
### ########################################################## ###

## FUNCTIONS TO CREATE Y's FOR SCC's in the shape of a dataframe with only one row. 
## We then loop it in order to get one row p/PPT (that is the second code paragraph which accompanies every function)
## IN SUMMARY: FUNCTION TO GET A SINGLE PPT DATA, TRANSFORM IT TO FUNCTIONAL DATA STYLE, AND THEN LOOP THAT PROCESS FOR ALL THE PPT'S YOU NEED
## This is not the cleanest way to do this but I find it the most static for a shareable script. We could also define the Z level as other argument
## of our function and then create a new matrix with every new loop (...)


# CN function:

YcreatorCN <- function(i){
  Y <- subset(database, database$PPT==names_CN[i] & database$z==30) # Y = Choose by PPT and Z=30, group=change depending on the SCC you need
  Y <- Y[1:7505,8] # Make sure to keep only 7505 (problematic before) and column 8 (normalized PET responses)
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        # CN LOOP: now every CN-PPT will be one row, as ImageSCC package requires
        
      SCC_matrix_CN <- matrix(ncol = 7505,nrow=0)
      for (i in 1:length(names_CN)){
        temp<-YcreatorCN(i)
        SCC_matrix_CN<-rbind(SCC_matrix_CN,temp)    
      }

    
# AD function:

YcreatorAD <- function(i){
  Y <- subset(database, database$PPT==names_AD[i] & database$z==30) 
  Y <- Y[1:7505,8] 
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        # AD LOOP: 
        
        SCC_matrix_AD <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_AD)){
          temp<-YcreatorAD(i)
          SCC_matrix_AD<-rbind(SCC_matrix_AD,temp)    
        }

        
### THE FOLLOWING FUNCTIONS AND LOOPS ARE FOR COMPARATIONS ACCORDING TO SEX:
                
    
# AD + Female function:

YcreatorAD_F <- function(i){
  Y <- subset(database, database$PPT==names_AD_F[i] & database$z==30) 
  Y <- Y[1:7505,8]
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        # AD + Female LOOP:
        
        SCC_matrix_AD_F <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_AD_F)){
          temp<-YcreatorAD_F(i)
          SCC_matrix_AD_F<-rbind(SCC_matrix_AD_F,temp)    
        }


# AD + Male function:

YcreatorAD_M <- function(i){
  Y <- subset(database, database$PPT==names_AD_M[i] & database$z==30) 
  Y <- Y[1:7505,8] 
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        
        # AD + Male LOOP: 
        
        SCC_matrix_AD_M <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_AD_M)){
          temp<-YcreatorAD_M(i)
          SCC_matrix_AD_M<-rbind(SCC_matrix_AD_M,temp)    
        }


    
# CN + Female function:

YcreatorCN_F <- function(i){
  Y <- subset(database, database$PPT==names_CN_F[i] & database$z==30) 
  Y <- Y[1:7505,8] 
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        # CN + Female LOOP: 
        
        SCC_matrix_CN_F <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_CN_F)){
          temp<-YcreatorCN_F(i)
          SCC_matrix_CN_F<-rbind(SCC_matrix_CN_F,temp)    
        }



# CN + Male function:

YcreatorCN_M <- function(i){
  Y <- subset(database, database$PPT==names_CN_M[i] & database$z==30) 
  Y <- Y[1:7505,8] 
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}
    
        # CN + Male LOOP: 
        
        SCC_matrix_CN_M <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_CN_M)){
          temp<-YcreatorCN_M(i)
          SCC_matrix_CN_M<-rbind(SCC_matrix_CN_M,temp)    
        }

        
        
### THE FOLLOWING FUNCTIONS AND LOOPS ARE FOR COMPARATIONS ACCORDING TO AGE:
                

# CN <75 function:

YcreatorCN_less_75 <- function(i){
  Y <- subset(database, database$PPT==names_CN_less_75[i] & database$z==30) 
  Y <- Y[1:7505,8] 
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        # CN <75 loop:
        
        SCC_matrix_CN_less_75 <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_CN_less_75)){
          temp<-YcreatorCN_less_75(i)
          SCC_matrix_CN_less_75<-rbind(SCC_matrix_CN_less_75,temp)    
        }


# CN >75 function:

YcreatorCN_more_75 <- function(i){
  Y <- subset(database, database$PPT==names_CN_more_75[i] & database$z==30) 
  Y <- Y[1:7505,8] 
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        # CN >75 loop:
        
        SCC_matrix_CN_more_75 <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_CN_more_75)){
          temp<-YcreatorCN_more_75(i)
          SCC_matrix_CN_more_75<-rbind(SCC_matrix_CN_more_75,temp)    
        }


# AD <75 function:

YcreatorAD_less_75 <- function(i){
  Y <- subset(database, database$PPT==names_AD_less_75[i] & database$z==30)
  Y <- Y[1:7505,8] 
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        # AD <75 loop:
        
        SCC_matrix_AD_less_75 <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_AD_less_75)){
          temp<-YcreatorAD_less_75(i)
          SCC_matrix_AD_less_75<-rbind(SCC_matrix_AD_less_75,temp)    
        }


# AD >75 function:

YcreatorAD_more_75 <- function(i){
  Y <- subset(database, database$PPT==names_AD_more_75[i] & database$z==30) 
  Y <- Y[1:7505,8] 
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

        # AD >75 loop:
        
        SCC_matrix_AD_more_75 <- matrix(ncol = 7505,nrow=0)
        for (i in 1:length(names_AD_more_75)){
          temp<-YcreatorAD_more_75(i)
          SCC_matrix_AD_more_75<-rbind(SCC_matrix_AD_more_75,temp)    
        }


        
### ########################################################## ###
#####                *CONTOURS OF NEURO-DATA*                 ####
### ########################################################## ###
        
        
## Now that we have our data prepared in the structure we need it to be, we also need the triangulation parameters.
## However, before generating the Delaunay Triangulation grid, we need to know the contour/boundaries of our data,
## otherwise we won't be adjusting to the true shape of our data and the results will be critically affected.

        
# Z are the coordinates where data is measured:

x <- rep(1:79, each = 95, length.out = 7505) 
y <- rep(1:95,length.out = 7505)
Z <- cbind(as.matrix(x),as.matrix(y)); Z


# Triangulation Parameters: this code is made for a custom triangulation for the position of our PET data
# The package provides some simple templates but for optimal performance we need to find the boundaries of our data MANUALLY
# If we want to extract contour for other slice (not Z=30) we would have to go back, change corresponding function+loop and then with
# the new SCC matrix, come back here and extract boundaries.

dat <- cbind(Z,as.matrix(SCC_matrix_AD[1,])) # We take a sample to calculate boundaries
dat <- as.data.frame(dat)
dat[is.na(dat)] <- 0
sum(is.na(dat$pet)) # should be = 0
head(dat)

  ## memory.size(max = TRUE) # In my case this was necessary using a laptop.


rownames(dat) <- NULL # Remove row numbering, for some reason this can be problematic sometimes

df = getContourLines(dat[1:7504,], # select all rows but for the last one
                     levels = c(0)) # and search for the jump between PET=0 and other value, that will be the boundary

ggplot(df,aes(x,y,colour = z)) + geom_path() # Display of brain boundaries for Z=30

contour30 = df # Contour in z=30 -> contour30 (change if using a different Z)
head(contour30);str(contour30)


### Now for the different sets of coordinates as we will need external boundaries and internal holes in a list:

contour <- function(x){
  
  aa <- contour30[contour30$GID == x,] # We keep GID==x (0,1,2,3...)
  a <- aa[,5:6] # Then we keep just the coordinates 
  print(a) # and then print in order to loop and make a list
}

coord <- list()

for (i in 0:max(contour30$GID)) { #change contour30 to any other name previously assigned if necessary
  
  coord[[i + 1]] <- contour(i)
  rownames(coord[[i + 1]]) <- NULL
}


### Test the results are coherent:

head(coord[[1]],10); plot(coord[[1]])   # external boundaries
head(coord[[2]],10); points(coord[[2]]) # first hole
head(coord[[3]],10); points(coord[[3]]) # second hole


### ########################################################## ###
#####             *TRIANGULATION PARAMETERS*                  ####
### ########################################################## ###

# Now we create the triangulation grid using these coordinates and a N fineness degree value
# N = An integer parameter controlling the fineness of the triangulation and subsequent triangulation. As n increases the fineness increases. Usually, n = 8 seems to be a good choice.


VT8 = TriMesh(coord[[1]],8,list(as.matrix(coord[[2]]),as.matrix(coord[[3]])))
head(VT8$V,10);head(VT8$Tr,10) # Vertices' coordinates and nodes you have to link to create the triangulation

# N=8 is usually a good enough fineness degree although you can use finer triangulations which will take more time to compute

VT15 = TriMesh(coord[[1]],15,list(as.matrix(coord[[2]]),as.matrix(coord[[3]])))
head(VT15$V,10);head(VT15$Tr,10)

VT25 = TriMesh(coord[[1]],25,list(as.matrix(coord[[2]]),as.matrix(coord[[3]])))
head(VT25$V,10);head(VT25$Tr,10)

par(pin = c(5.25,6.33),
    mai = c(0.5,0.5,0.5,0.5))

TriPlot(VT8$V, VT8$Tr, col = 1, lwd = 1)
TriPlot(VT15$V, VT15$Tr, col = 1, lwd = 1)
TriPlot(VT25$V, VT25$Tr, col = 1, lwd = 1)

# In order to be consistent we use common names Brain.V and Brain.Tr. From here onwards most of the names follow the ones provided by Wang et al (2019)

Brain.V <- VT15[[1]]
Brain.Tr <- VT15[[2]]

head(Brain.V);head(Brain.Tr)

V.est = as.matrix(Brain.V) # Brain.v <- cbind(Brain.V[,2],Brain.V[,1]) # In case you need to transpose the data
Tr.est = as.matrix(Brain.Tr)

V.band=as.matrix(Brain.V)
Tr.band=as.matrix(Brain.Tr) 

# Response Variable:

Y_CN = SCC_matrix_CN
Y_AD = SCC_matrix_AD

Y_AD_F = SCC_matrix_AD_F
Y_AD_M = SCC_matrix_AD_M

Y_CN_F = SCC_matrix_CN_F
Y_CN_M = SCC_matrix_CN_M

Y_CN_L_75 = SCC_matrix_CN_less_75
Y_CN_M_75 = SCC_matrix_CN_more_75

Y_AD_L_75 = SCC_matrix_AD_less_75
Y_AD_M_75 = SCC_matrix_AD_more_75

### ########################################################## ###
#####        *OTHER PARAMETERS FOR SCC ESTIMATION*            ####
### ########################################################## ###

# Following Wang et al's recomendations:

d.est = 5 # degree of spline for mean function 
d.band = 2 # degree of spline for SCC
r = 1 # smoothing parameter
lambda = 10^{seq(-6,3,0.5)} # penalty parameters
alpha.grid = c(0.10,0.05,0.01) # vector of confidence levels

### ########################################################## ###
#####               *CONSTRUCTION OF SCC'S*                   ####
### ########################################################## ###

# Run one sample SCC construction:

x <- rep(1:91, each = 109, length.out = 9919) 
y <- rep(1:109, length.out = 9919)
Z <- cbind(as.matrix(x),as.matrix(y)); Z
 
SCC_CN_1 = scc.image(Ya = SCC_CN,Z = Z,d.est = d.est,d.band = d.band,r = r,
                     V.est.a = V.est,Tr.est.a = Tr.est,
                     V.band.a = V.band,Tr.band.a = Tr.band,
                     penalty = TRUE,lambda = lambda,alpha.grid = alpha.grid,adjust.sigma = TRUE)

  plot(SCC_CN_1,
       breaks=seq(from=0,to=2,length.out = 65),
       col=,
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Control Group",
       col.sub="black",
       family ="serif")
  
    # This line of code shows triangulations and then SCCs for different alpha levels in a single-group 


SCC_AD_1=scc.image(Ya=Y_AD,Z=Z,d.est=d.est,d.band=d.band,r=r,
                   V.est.a=V.est,Tr.est.a=Tr.est,
                   V.band.a=V.band,Tr.band.a=Tr.band,
                   penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)


  plot(SCC_AD_1,
       breaks=seq(from=0,to=2,length.out = 65),
       col=,
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Alzheimer Group",
       col.sub="black",
       family ="serif")

   # Same here for the other group (AD). These images are orientative, our main goal is obtaining SCCs for the difference
   # between groups' mean functions.
  
   ## IN CASE YOU WANT TO EXPORT PLOTS THIS IS EASIER AFTER (...) :  

  par(las=1,
          col.main="white",
          col.lab="white",
          col.sub="white")
      
      plot(SCC_AD_1,
           axes=FALSE,
           breaks=seq(from=0,to=1.6,length.out = 65), #1.6 because max(SCC_AD_1$scc) and max(SCC_CN_1$scc) are both below
           col=,
           family ="serif",
           horizontal=F)
      dev.off()



########################################################################################/
#   AND NOW FOR A COMPARATION BETWEEN CN AND AD WHICH IS THE IMPORTANT SETUP            /
########################################################################################/ 

      

SCC_COMP_1=scc.image(Ya=Y_AD,Yb=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                     V.est.a=V.est,Tr.est.a=Tr.est,
                     V.band.a=V.band,Tr.band.a=Tr.band,
                     penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


# With Yb included we now have a Two-group SCC for the difference between mean functions of Yb and Ya

## IMPORTANT NOTE: ESTIMATED MEAN FUNCTIONS APPEAR IN THE ORDER THEY APPEAR IN THE CODE (FIRST Ya THEN Yb)
## BUT SCC IS ACTUALLY CALCULATED FOR THE DIFFERENCE BETWEEN Yb MINUS Ya

    
    # FIRST WAY OF VISUALIZATION:

    plot(SCC_COMP_1,
         # breaks=c(0,2),
         # col="turquoise",
         breaks=seq(from=0,to=2,length.out = 65),
         xlab="Longitudinal (1-95)",
         ylab="Transversal (1-79)",
         sub="Difference between estimated mean functions: CNs - ADs",
         col.sub="red",
         family ="serif")

    # SECOND WAY: JUST ONE COLOR AND THEN WE OVERLAY A SERIES OF POINTS
    # RUN THIS NEXT PLOT CODE, STOP IN ONE OF THE ESTIMATED MEAN FUNCTIONS, THEN RUN "POINTS" TO OVERLAY THEM
    
    plot(SCC_COMP_1,
         breaks=c(0,2),
         col="turquoise",
         # breaks=seq(from=0,to=2,length.out = 65),
         xlab="Longitudinal (1-95)",
         ylab="Transversal (1-79)",
         sub="Difference between estimated mean functions: CNs - ADs",
         col.sub="red",
         family ="serif")


    
my_points <- function(aa){
  
  Z.band <- matrix(aa$Z.band,ncol=2) # Positions
  z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Separated positions
  n1 <- length(z1); n2 <- length(z2) # Lengths of those positions
  scc <- matrix(NA,n1*n2,2) # Matrix with value of that SCC in that position
  ind.inside.band <- aa$ind.inside.cover # Keep only regions inside triangulation
  scc[ind.inside.band,] <- aa$scc[,,2] # Assign SCC to those areas
  scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # LIMITS: minimum of inferior, maximum of superior
  
  scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
  scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
  scc.l.mtx[scc.l.mtx<0]=NA # The ones that work just fine are substituded by a NA as we don't want to represent them, only positive LowerSCCs are represented
  scc.u.mtx[scc.u.mtx>0]=NA # The ones that work just fine are substituded by a NA as we don't want to represent them, only positive UpperSCCs are represented
  
  points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Points where mean difference is positive (first image is stronger)
  points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Points where mean difference is negative (second image is stronger)
  
  pointers<-list(points.P,points.N)
  print(pointers)  
}


  # NOT NECESSARY TO RUN THESE TWO LINES
  image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regions where the difference between mean functions is positive (falls above 0). That is, activity in Image1 is stronger than in Image 2 in that area
  image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regions where the difference between mean functions is negative (falls below 0). That is, activity in Image2 is stronger than in Image 1 in that area 
  #

points_1 <- my_points(SCC_COMP_1) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)

## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC

  points(points_1[[1]],
         type="p",
         pch=15,
         col="navy")  
  
  points(points_1[[2]],
         type="p",
         pch=15,
         col="yellow") 

    # Type: p,l,b,c,o,s,S,h,n
    # pch = 0:18 =:46
    # col=
    # bg= background
    # lwd= line width


##############################//

  
SCC_COMP_2=scc.image(Ya=Y_CN,Yb=Y_AD,Z=Z,d.est=d.est,d.band=d.band,r=r,
                     V.est.a=V.est,Tr.est.a=Tr.est,
                     V.band.a=V.band,Tr.band.a=Tr.band,
                     penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


  plot(SCC_COMP_2, # Comment args 1&2 and uncoment arg 3 to change visualization
       breaks=c(0,2),  
       col="turquoise",
       #breaks=seq(from=0,to=2,length.out = 65),
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Difference between estimated mean functions: CNs - ADs",
       col.sub="red",
       family ="serif")


  points_2 <- my_points(SCC_COMP_2) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)
  
  ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
  
      points(points_1[[1]],
             type="p",
             pch=15,
             col="navy")  
      
      points(points_1[[2]],
             type="p",
             pch=15,
             col="yellow") 
    
      # Type: p,l,b,c,o,s,S,h,n
      # pch = 0:18 =:46
      # col=
      # bg= background
      # lwd= line width



####################################################/
#   COMPARISON BETWEEN AD FEMALES AND AD MALES
####################################################/ 


SCC_COMP_MF1=scc.image(Ya=Y_F,Yb=Y_M,Z=Z,d.est=d.est,d.band=d.band,r=r,
                       V.est.a=V.est,Tr.est.a=Tr.est,
                       V.band.a=V.band,Tr.band.a=Tr.band,
                       penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


  plot(SCC_COMP_MF1,
       breaks=c(0,2),
       col="turquoise",
       # breaks=seq(from=0,to=2,length.out = 65),
       # col=,
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Difference between AD estimated mean functions for: Male AD - Female AD",
       col.sub="red",
       family ="serif")


    points_MF1 <- my_points(SCC_COMP_MF1) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)
      
        ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
        
        points(points_MF1[[1]],
               type="p",
               pch=15,
               col="navy")  
        
        points(points_MF1[[2]],
               type="p",
               pch=15,
               col="yellow") 
        
        # Type: p,l,b,c,o,s,S,h,n
        # pch = 0:18 =:46
        # col=
        # bg= background
        # lwd= line width


##############################//

SCC_COMP_MF2=scc.image(Ya=Y_M,Yb=Y_F,Z=Z,d.est=d.est,d.band=d.band,r=r,
                       V.est.a=V.est,Tr.est.a=Tr.est,
                       V.band.a=V.band,Tr.band.a=Tr.band,
                       penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


    plot(SCC_COMP_MF2,
         breaks=c(0,2),
         col="turquoise",
         #breaks=seq(from=0,to=2,length.out = 65),
         xlab="Longitudinal (1-95)",
         ylab="Transversal (1-79)",
         sub="Difference between AD estimated mean functions for: Male AD - Female AD",
         col.sub="red",
         family ="serif")


     points_MF2 <- my_points(SCC_COMP_MF2) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)
    
          ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
          
          points(points_MF2[[1]],
                 type="p",
                 pch=15,
                 col="navy")  
          
          points(points_MF2[[2]],
                 type="p",
                 pch=15,
                 col="yellow") 
          
          # Type: p,l,b,c,o,s,S,h,n
          # pch = 0:18 =:46
          # col=
          # bg= background
          # lwd= line width
    
    
          
################################################################################/
#   AND NOW FOR A COMPARATION BETWEEN AD male AND CN male AND SAME FOR FEMALE
################################################################################/



SCC_COMP_M_ADCN=scc.image(Ya=Y_CN_M,Yb=Y_AD_M,Z=Z,d.est=d.est,d.band=d.band,r=r,
                          V.est.a=V.est,Tr.est.a=Tr.est,
                          V.band.a=V.band,Tr.band.a=Tr.band,
                          penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    

    
    plot(SCC_COMP_M_ADCN,
         breaks=c(0,2),
         col="turquoise",
         # breaks=seq(from=0,to=2,length.out = 65),
         xlab="Longitudinal (1-95)",
         ylab="Transversal (1-79)",
         sub="Difference between Control males and Alzheimer males.",
         col.sub="red",
         family ="serif")

    
    points_M_ADCN <- my_points(SCC_COMP_M_ADCN) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)

        ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
        
        points(points_M_ADCN[[1]],
               type="p",
               pch=15,
               col="navy")  
        
        points(points_M_ADCN[[2]],
               type="p",
               pch=15,
               col="yellow") 
        
        # Type: p,l,b,c,o,s,S,h,n
        # pch = 0:18 =:46
        # col=
        # bg= background
        # lwd= line width


#########/


SCC_COMP_F_ADCN=scc.image(Ya=Y_CN_F,Yb=Y_AD_F,Z=Z,d.est=d.est,d.band=d.band,r=r,
                          V.est.a=V.est,Tr.est.a=Tr.est,
                          V.band.a=V.band,Tr.band.a=Tr.band,
                          penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


  plot(SCC_COMP_F_ADCN,
       breaks=c(0,2),
       col="turquoise",
       # breaks=seq(from=0,to=2,length.out = 65),
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Difference between Control females and Alzheimer females",
       col.sub="red",
       family ="serif")

  
      points_F_ADCN <- my_points(SCC_COMP_F_ADCN) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)
      
      ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
      
      points(points_F_ADCN[[1]],
             type="p",
             pch=15,
             col="yellow")  
      
      points(points_F_ADCN[[2]],
             type="p",
             pch=15,
             col="navy") 
      
      # Type: p,l,b,c,o,s,S,h,n
      # pch = 0:18 =:46
      # col=
      # bg= background
      # lwd= line width
      
  
################################################################################/
#   AND NOW FOR A COMPARATION BETWEEN AD AND CN FILTERING BY AGE BLOCK
################################################################################/


# CN VS AD <75
      
SCC_COMP_LESS_75=scc.image(Ya=Y_AD_L_75,Yb=Y_CN_L_75,Z=Z,d.est=d.est,d.band=d.band,r=r,
                           V.est.a=V.est,Tr.est.a=Tr.est,
                           V.band.a=V.band,Tr.band.a=Tr.band,
                           penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


  plot(SCC_COMP_LESS_75,
       # breaks=c(0,2),
       # col="turquoise",
       breaks=seq(from=0,to=2,length.out = 65),
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Difference between Control and Alzheimer for <75 y.o. patients",
       col.sub="red",
       family ="serif")

  
      points_LESS_75 <- my_points(SCC_COMP_LESS_75) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)
      
      ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
      
      points(points_LESS_75[[1]],
             type="p",
             pch=15,
             col="navy")  
      
      points(points_LESS_75[[2]],
             type="p",
             pch=15,
             col="yellow") 
      
      # Type: p,l,b,c,o,s,S,h,n
      # pch = 0:18 =:46
      # col=
      # bg= background
      # lwd= line width


#########################/

# ONLY MALES <75

SCC_COMP_MORE_75=scc.image(Ya=Y_AD_M_75,Yb=Y_CN_M_75,Z=Z,d.est=d.est,d.band=d.band,r=r,
                           V.est.a=V.est,Tr.est.a=Tr.est,
                           V.band.a=V.band,Tr.band.a=Tr.band,
                           penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    



    plot(SCC_COMP_MORE_75,
         # breaks=c(0,2),
         # col="turquoise",
         breaks=seq(from=0,to=2,length.out = 65),
         xlab="Longitudinal (1-95)",
         ylab="Transversal (1-79)",
         sub="Difference between Control and Alzheimer for >75 y.o. patients",
         col.sub="red",
         family ="serif")

  
        points_MORE_75 <- my_points(SCC_COMP_MORE_75) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)
    
        ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
        
        points(points_MORE_75[[1]],
               type="p",
               pch=15,
               col="navy")  
        
        points(points_MORE_75[[2]],
               type="p",
               pch=15,
               col="yellow") 
        
        # Type: p,l,b,c,o,s,S,h,n
        # pch = 0:18 =:46
        # col=
        # bg= background
        # lwd= line width

        

#############################    /

        
# COMPARING AD'S IN DIFFERENT AGE BLOCKS        

SCC_COMP_ADS_AGE=scc.image(Ya=Y_AD_M_75,Yb=Y_AD_L_75,Z=Z,d.est=d.est,d.band=d.band,r=r,
                           V.est.a=V.est,Tr.est.a=Tr.est,
                           V.band.a=V.band,Tr.band.a=Tr.band,
                           penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    



  plot(SCC_COMP_ADS_AGE,
       breaks=c(0,2),
       col="turquoise",
       # breaks=seq(from=0,to=2,length.out = 65),
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Difference between Control and Alzheimer for >75 y.o. patients",
       col.sub="red",
       family ="serif")

  
      points_ADS_AGE <- my_points(SCC_COMP_ADS_AGE) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)
      
      ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
      
      points(points_ADS_AGE[[1]],
             type="p",
             pch=15,
             col="navy")  
      
      points(points_ADS_AGE[[2]],
             type="p",
             pch=15,
             col="yellow") 
      
      # Type: p,l,b,c,o,s,S,h,n
      # pch = 0:18 =:46
      # col=
      # bg= background
      # lwd= line width


##################/          


# COMPARING CN'S OF DIFFERENT AGE BLOCKS

SCC_COMP_CNS_AGE=scc.image(Ya=Y_CN_M_75,Yb=Y_CN_L_75,Z=Z,d.est=d.est,d.band=d.band,r=r,
                           V.est.a=V.est,Tr.est.a=Tr.est,
                           V.band.a=V.band,Tr.band.a=Tr.band,
                           penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    

  
  plot(SCC_COMP_CNS_AGE,
       # breaks=c(0,2),
       # col="turquoise",
       breaks=seq(from=0,to=2,length.out = 65),
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Difference between Control and Alzheimer for >75 y.o. patients",
       col.sub="red",
       family ="serif")

  
      points_CNS_AGE <- my_points(SCC_COMP_CNS_AGE) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)
      
      ## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC
      
      points(points_CNS_AGE[[1]],
             type="p",
             pch=15,
             col="navy")  
      
      points(points_CNS_AGE[[2]],
             type="p",
             pch=15,
             col="yellow") 
      
      # Type: p,l,b,c,o,s,S,h,n
      # pch = 0:18 =:46
      # col=
      # bg= background
      # lwd= line width


