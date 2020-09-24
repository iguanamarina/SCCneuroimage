#######################################################################\
#   
#   THE OBJECTIVES OF THIS MASTER-SCRIPT ARE:
#     
#     - TO IMPORT NIFTI FILES INTO 'R'
#     - TO CREATE AN ORDERED DATABASE FROM A SERIES OF PET IMAGES (.NIFTI)
#     - 
#     - TO GET A LIST OF THE PPT'S NAMES        
#     - CREATE A NEW DATABASE FOR SCCs
#     - INCLUDE NEW TRIANGULATION PARAMETERS
#     - TEST SCC'S WITH DEFINITIVE DATA
#
#
#######################################################################\


### ################################################## ###
#####                NEUROIMAGE DATA                  ####
### ################################################## ###


## PREVIOUS PROCESSING OF PET IMAGES (REALIGNEMENT, WRAPPING, CORREGISTRARION, NORMALIZATION, AND MASKING) DONE IN 'MATLAB'
## CHECK OUT: ../SCCneuroimage/Matlab Preprocessing Code/
## THESE ARE EXTENSIONS OF SPM MATLAB CODE WITH SLIGHT MODIFICATIONS IN ORDER TO LOOP THE PROCESS OVER FOR OUR 126 PPT'S
  

### ################################################## ###
#####                 *PREAMBLE*                      ####
####         Installation of necessary packgs.         ###
### ################################################## ###


install.packages(c("mgcv","gamair","oro.nifti","memsic"))
library(mgcv);library(gamair);library(oro.nifti);library(memisc)
memory.size(max = TRUE) # Not useful in Linux systems

### ################################################## ###
#####            *IMPORT NIFTI FILES*                 ####
####            f.clean  (PPT,z,x,y,pet)               ###
### ################################################## ###


## Set as working directory -> directory with .hdr/.img NIfTI files, mask, and Demographics .csv file

setwd("/Users/Juan A. Arias/Documents/GitHub/SCCneuroimage/PETimg_masked") # My Directory

demo <- read.csv2("Demographics.csv")

f.clean <- function(name) { #### f.clean is meant for CLEANING ONE SINGLE-PPT DATA, then we loop it
  
  # Read NIFTI image, transform it to dataframe, preserve slice Z and organize the table
  
  ## Load Data
  
  file <- readNIfTI(fname = name, verbose = FALSE, warn = -1, reorient = TRUE, call = NULL, read_data = TRUE)
  namex <- as.character(name)
  n = img_data(file)
  n = to.data.frame(n)
  
  ## Prepare data.frame base where further data from the loop will be integrated
  
  dataframe <- data.frame(z=integer(),x=integer(),y=integer(),pet=integer()) 
  
  # Loop for 79 slices of Z in the NiFtI image -> move to dataframe
  
  for (i in seq(1:79)) {
    
    n_lim = n[n$Var2==i,] # Select just one Z slice
    n_lim$Var1=NULL
    n_lim$Var2=NULL
    
    z <-rep(i,length.out=7505)
    x <-rep(1:79, each=95, length.out = 7505) 
    y <-rep(1:95,length.out = 7505)
    attach(n_lim)
    pet<-c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`,`13`,`14`,`15`,`16`,`17`,`18`,`19`,`20`,
           `21`,`22`,`23`,`24`,`25`,`26`,`27`,`28`,`29`,`30`,`31`,`32`,`33`,`34`,`35`,`36`,`37`,`38`,`39`,`40`,
           `41`,`42`,`43`,`44`,`45`,`46`,`47`,`48`,`49`,`50`,`51`,`52`,`53`,`54`,`55`,`56`,`57`,`58`,`59`,`60`,
           `61`,`62`,`63`,`64`,`65`,`66`,`67`,`68`,`69`,`70`,`71`,`72`,`73`,`74`,`75`,`76`,`77`,`78`,`79`)
    detach(n_lim)
    
    temp0 = data.frame(z,x,y,pet) # temporal dataframe
    temp1 <- print(temp0) # unsure whether this is necessary but, if things work, don't touch them
    dataframe <- rbind(dataframe,temp1) # sum new data with previous data
    
  }
  
  # Demographics: PPT, group (AD/CN), sex, age.
  
  demog <- demo[demo$PPT==namex,]
  
  PPT <- rep(demog$PPT,length.out=7505)
  group <-rep(demog$Group,length.out=7505)
  sex <-rep(demog$Sex,length.out=7505)
  age <-rep(demog$Age,length.out=7505)
  
  temp2 <- data.frame(PPT,group,sex,age)
  dataframe <- cbind(temp2,dataframe)
  
  print(dataframe) # Necessary for assigning an object name

}


# Example of conversion from NIFTI to R dataframe:

example = f.clean("003_S_1059")
head(example)  # Some values are Zeros due to the masking process


### ################################################## ###
#####             *CREATE   DATABASE*                 ####
### ################################################## ###


files <- list.files(path=getwd(), pattern="*.img", full.names=F, recursive=FALSE) # list of files
files <- gsub(files, pattern=".img$", replacement="") # remove file extension .img
files

database <- data.frame(PPT=integer(),group=integer(),sex=integer(),age=integer(),z=integer(),x=integer(),y=integer(),pet=integer())
#create data.frame to include data

for (i in 1:length(files)) { #loop to include every PPT in the dataframe
  
  temporal <- f.clean(files[i])
  database <- rbind(database,temporal)
  
}

nrow(database[database$pet<0,]) # There are negative values (ilogical)

database$pet[database$pet<0]<- NaN  # Convert negatives to NaN 
database$pet[database$pet==0] <- NaN # Convert zeros to NaN 


### ####################################################### ###
#####           *MEAN AVERAGE NORMALIZATION*               ####
### ####################################################### ###


mean.data <- data.frame(pet=integer(),pet_normal=integer())

for (i in 1:length(files)) {
  
  temp <- database[database$PPT==files[i],8]
  #temp <- temp[,8]
  mean <- mean(as.numeric(temp),na.rm=T)
  pet_normal <- as.data.frame(temp/mean)
  colnames(pet_normal)<-c("pet_normal")
  temp <- cbind(temp,pet_normal)
  show(temp)    
  mean.data <- rbind(mean.data,temp)
  
}

mean.data <- mean.data[,2]
database <- cbind(database,mean.data)
colnames(database)<-c("PPT","group","sex","age","z","x","y","pet","pet_normal")

# Export complete database (with non normalized and normalized values)
## Set as working directory -> directory where you want to save Database.Rdata

setwd("/Users/Juan A. Arias/Documents/GitHub/SCCneuroimage") 
write.table(database,file="Database.csv",sep=";",na="NA") #export as .csv IF DESIRABLE
save(database,file="database.Rda") #export as .Rda IF DESIRABLE


# now to continue with our workflow:

database<-database[,-8]   # Usually you will be working with normalized PET values but you can skip this step 
                          # if you want to keep the non-normalized values for something


### ########################################################## ###
#####          SIMULTANEOUS CONFIDENCE CORRIDORS              ####
### ########################################################## ###


# load("C:/Users/Juan A. Arias/Documents/GitHub/SCCneuroimage/database.Rda") # if not loaded already

# This part of the code needs a series of packages which to date (09/2020) are only available at GitHub

# Install packages from GitHub: #

install.packages(c("devtools","remotes","readr","imager","itsadug","ggplot2","contoureR"))
library(devtools);library(remotes);library(readr);library(imager);library(itsadug);library(ggplot2);library(contoureR)

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE) # Forces installation to finish even with warning messages (EXTREMELY NECESSARY)
remotes::install_github("funstatpackages/BPST")
remotes::install_github("funstatpackages/Triangulation")
remotes::install_github("funstatpackages/ImageSCC")

library(BPST);library(Triangulation);library(ImageSCC)



# Preliminary modifications and exploratory analysis:

Data<-database
attach(Data)
head(Data)
str(Data)


### ########################################################## ###
#####           *COMPREHENSIVE LISTS OF PPT'S*                ####
### ########################################################## ###


# Working in a Functional Data Analysis (FDA) framework requires some transformations of the data. Namely, we will need to
# create data.frames in which every PPT corresponds to one row and each column to one PET value out of 7505 in a slice. 
# This is the way of approaching the data in FDA and obviously we will be limited to one brain slice at a time. 
# In order to create the SCC matrices we will also need the PPTs names listed depending on their sex, group... that's
# why we do it here in advance for AD, CN, female & male, <75 & >75... and combinations.


names<-as.vector(Data[,1])
names<-names[duplicated(names)!=TRUE]
((length(names))*592895)==nrow(Data)  # Every PPT has 592.895 measures. PPT*measures=nrow(Data). Should return a TRUE


# LIST OF CONTROL PPT'S #

names_CN<-as.vector(Data[,1])[Data$group=="CN"]
names_CN<-names_CN[duplicated(names_CN)!=TRUE]; names_CN

# LIST OF ALZHEIMER DISEASE PPT'S #

names_AD<-as.vector(Data[,1])[Data$group=="AD"]
names_AD<-names_AD[duplicated(names_AD)!=TRUE]; names_AD

# LIST OF MALE CONTROL PPT'S #

names_CN_M<-as.vector(Data[,1])[Data$sex=="M" & Data$group=="CN"]
names_CN_M<-names_CN_M[duplicated(names_CN_M)!=TRUE]; names_CN_M

# LIST OF FEMALE CONTROL PPT'S #

names_CN_F<-as.vector(Data[,1])[Data$sex=="F" & Data$group=="CN"]
names_CN_F<-names_CN_F[duplicated(names_CN_F)!=TRUE]; names_CN_F

# LIST OF FEMALE ALZHEIMER PPT'S #

names_AD_F<-as.vector(Data[,1])[Data$sex=="F" & Data$group=="AD"]
names_AD_F<-names_AD_F[duplicated(names_AD_F)!=TRUE];names_AD_F

# LIST OF MALE ALZHEIMER PPT'S #

names_AD_M<-as.vector(Data[,1])[Data$sex=="M" & Data$group=="AD"]
names_AD_M<-names_AD_M[duplicated(names_AD_M)!=TRUE]; names_AD_M

# LIST OF ALZHEIMER PPT'S <=75 Y.O. #

names_AD_less_75<-as.vector(Data[,1])[Data$age<=75 & Data$group=="AD"]
names_AD_less_75<-names_AD_less_75[duplicated(names_AD_less_75)!=TRUE]; names_AD_less_75

# LIST OF ALZHEIMER PPT'S >75 Y.O. #

names_AD_more_75<-as.vector(Data[,1])[Data$age>75 & Data$group=="AD"]
names_AD_more_75<-names_AD_more_75[duplicated(names_AD_more_75)!=TRUE]; names_AD_more_75

# LIST OF CONTROL PPT'S <=75 Y.O. #

names_CN_less_75<-as.vector(Data[,1])[Data$age<=75 & Data$group=="CN"]
names_CN_less_75<-names_CN_less_75[duplicated(names_CN_less_75)!=TRUE]; names_CN_less_75

# LIST OF CONTROL PPT'S >75 Y.O. #

names_CN_more_75<-as.vector(Data[,1])[Data$age>75 & Data$group=="CN"]
names_CN_more_75<-names_CN_more_75[duplicated(names_CN_more_75)!=TRUE]; names_CN_more_75



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
  Y <- subset(Data, Data$PPT==names_CN[i] & Data$z==30) # Y = Choose by PPT and Z=30, group=change depending on the SCC you need
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
  Y <- subset(Data, Data$PPT==names_AD[i] & Data$z==30) 
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
  Y <- subset(Data, Data$PPT==names_AD_F[i] & Data$z==30) 
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
  Y <- subset(Data, Data$PPT==names_AD_M[i] & Data$z==30) 
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
  Y <- subset(Data, Data$PPT==names_CN_F[i] & Data$z==30) 
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
  Y <- subset(Data, Data$PPT==names_CN_M[i] & Data$z==30) 
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
  Y <- subset(Data, Data$PPT==names_CN_less_75[i] & Data$z==30) 
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
  Y <- subset(Data, Data$PPT==names_CN_more_75[i] & Data$z==30) 
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
  Y <- subset(Data, Data$PPT==names_AD_less_75[i] & Data$z==30)
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
  Y <- subset(Data, Data$PPT==names_AD_more_75[i] & Data$z==30) 
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

x <-rep(1:79, each=95, length.out = 7505) 
y <-rep(1:95,length.out = 7505)
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
                     levels=c(0)) # and search for the jump between PET=0 and other value, that will be the boundary

ggplot(df,aes(x,y,colour=z)) + geom_path() # Display of brain boundaries for Z=30

contour30=df # Contour in z=30 -> contour30 (change if using a different Z)
head(contour30);str(contour30)


### Now for the different sets of coordinates as we will need external boundaries and internal holes in a list:

contour<-function(x){
  
  aa<-contour30[contour30$GID==x,] # We keep GID==x (0,1,2,3...)
  a<-aa[,5:6] # Then we keep just the coordinates 
  print(a) # and then print in order to loop and make a list
}

coord<-list()

for (i in 0:max(contour30$GID)){ #change contour30 to any other name previously assigned if necessary
  
  coord[[i+1]]<-contour(i)
  rownames(coord[[i+1]]) <- NULL
}


### Test the results are coherent:

head(coord[[1]],10); plot(coord[[1]])   # external boundaries
head(coord[[2]],10); points(coord[[2]]) # first hole
head(coord[[3]],10); points(coord[[3]]) # second hole


### ########################################################## ###
#####             *TRIANGULATION PARAMETERS*                  ####
### ########################################################## ###


# An integer parameter controlling the fineness of the triangulation and subsequent triangulation. As n increases the fineness increases. Usually, n = 8 seems to be a good choice.

library(Triangulation)

VT8=TriMesh(coord[[1]],8,list(as.matrix(coord[[2]]),as.matrix(coord[[3]])))
head(VT8$V,10);head(VT8$Tr,10)

VT15=TriMesh(coord[[1]],15,list(as.matrix(coord[[2]]),as.matrix(coord[[3]])))
head(VT15$V,10);head(VT15$Tr,10)

VT25=TriMesh(coord[[1]],25,list(as.matrix(coord[[2]]),as.matrix(coord[[3]])))
head(VT25$V,10);head(VT25$Tr,10)


# In order to be consistent with the rest of the code:

#(!!)PLAY AROUND HERE
Brain.V <- VT15[[1]]
Brain.Tr <- VT15[[2]]

head(Brain.V);head(Brain.Tr)


V.est=as.matrix(Brain.V)
# Brain.v<-cbind(Brain.V[,2],Brain.V[,1])
Tr.est=as.matrix(Brain.Tr)

V.band=as.matrix(Brain.V)
Tr.band=as.matrix(Brain.Tr) 


# Response Variable (multiple):

Y_CN=SCC_matrix_CN
Y_AD=SCC_matrix_AD

Y_AD_F=SCC_matrix_AD_F
Y_AD_M=SCC_matrix_AD_M

Y_CN_F=SCC_matrix_CN_F
Y_CN_M=SCC_matrix_CN_M

Y_CN_L_75=SCC_matrix_CN_less_75
Y_CN_M_75=SCC_matrix_CN_more_75

Y_AD_L_75=SCC_matrix_AD_less_75
Y_AD_M_75=SCC_matrix_AD_more_75

# Parameters for SCC estimation: 

#(!!) play around with these

d.est=5 # degree of spline for mean function 
d.band=2 # degree of spline for SCC
r=1 # smoothing parameter
lambda=10^{seq(-6,3,0.5)} # penalty parameters
alpha.grid=c(0.10,0.05,0.01) # vector of confidence levels


# Run one sample SCC construction one time:

library(ImageSCC)
library(fields)

SCC_CN_1=scc.image(Ya=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                   V.est.a=V.est,Tr.est.a=Tr.est,
                   V.band.a=V.band,Tr.band.a=Tr.band,
                   penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)

plot(SCC_CN_1,
     breaks=seq(from=0,to=2,length.out = 65),
     col=,
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Control Group",
     col.sub="black",
     family ="serif")


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

par()
par(mfrow=c(1,1),
    col.main="white",
    col.lab="white",
    col.sub="white",
    fg="black",
    mfrow=c(1,1),
    mar=c(2,4,4,2),
    pin=c(9,7),
    las=3)

par(las=1)

plot(SCC_AD_1,
     axes=FALSE,
     breaks=seq(from=0,to=2,length.out = 65),
     col=,
     family ="serif",
     horizontal=F)

dev.off()



####################################################/
#   AND NOW FOR A COMPARATION BETWEEN CN AND AD
####################################################/ 


SCC_COMP_1=scc.image(Ya=Y_AD,Yb=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                     V.est.a=V.est,Tr.est.a=Tr.est,
                     V.band.a=V.band,Tr.band.a=Tr.band,
                     penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    
# With Yb included: Two-group SCC is constructed for the difference between the mean functions of Yb and Ya

## IMPORTANTISIMO: LAS E.M.F APARECEN EN EL ORDEN EN QUE SE PONEN EN SCC.IMAGE (PRIMERO Ya Y LUEGO Yb)
## PERO LOS SCC SE HACEN PARA LA DIFERENCIA ENTRE Yb Y Ya (!!!!!!!!!)

plot(SCC_COMP_1,
     breaks=c(0,2),
     col="turquoise",
     #breaks=seq(from=0,to=2,length.out = 65),
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Difference between estimated mean functions: CNs - ADs",
     col.sub="red",
     family ="serif")


aa=SCC_COMP_1

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan
image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=".",
       col="navy",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="orange",
       cex=9) 

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width



SCC_COMP_2=scc.image(Ya=Y_CN,Yb=Y_AD,Z=Z,d.est=d.est,d.band=d.band,r=r,
                     V.est.a=V.est,Tr.est.a=Tr.est,
                     V.band.a=V.band,Tr.band.a=Tr.band,
                     penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


plot(SCC_COMP_2,
     breaks=c(0,2),
     col="turquoise",
     #breaks=seq(from=0,to=2,length.out = 65),
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Difference between estimated mean functions: CNs - ADs",
     col.sub="red",
     family ="serif")

plot(SCC_COMP_2,
     breaks=c(0,2),
     col="turquoise",
     #breaks=seq(from=0,to=2,length.out = 65),
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Difference between estimated mean functions: CNs - ADs",
     col.sub="red",
     family ="serif")  

# Si haces la diferencia de medias y no hay diferencias significativas, 
# el cero (0) tiene que estar incluido en los intervalos propuestos, indicando 
# que la diferencia es pequeña y que bien podr�???a no ser ninguna diferencia en absoluto. 
# En cambio, si en el intervalo de confianza inferior hay valores superiores a 0, esto 
# indica que en esa ubicación (en ese pixel/voxel) el intervalo de confianza está por 
# encima del 0 (0 cae debajo) y que forzosamente en ese punto el valor será superior al 
# esperable por una diferencia de medias de imagenes 'iguales', o sea que el valor en la 
# primera imagen era más alto que en la segunda. De la misma forma, en el intervalo de 
# confianza superior pueden aparecer valores inferiores a cero (0) que indican que el 
# cero va a caer encima del intervalo de confianza, indicando que el valor real en ese 
# pixel/voxel es inferior al esperable por la diferencia de dos imagenes 'iguales' y que 
# en la diferencia entre medias, la segunda imagen era más potente o con valores más altos 
# en esa posición.


aa=SCC_COMP_2

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan
image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=2,
       col="yellow")  

points(points.N,
       type="p",
       pch=2,
       col="navy") 

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


aa=SCC_COMP_MF1

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)


points(points.P,
       type="p",
       pch=".",
       col="navy",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="orange",
       cex=9)

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width


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


aa=SCC_COMP_MF2

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=2,
       col="yellow")  

points(points.N,
       type="p",
       pch=2,
       col="orange") 

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width

################################################################################/
#   AND NOW FOR A COMPARATION BETWEEN AD male AND CN male AND SAME FOR female
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

aa=SCC_COMP_M_ADCN

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=".",
       col="orange",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="navy",
       cex=9)

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width


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

aa=SCC_COMP_F_ADCN

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=".",
       col="orange",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="navy",
       cex=9) 

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width

#############################/


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

aa=SCC_COMP_LESS_75

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=".",
       col="navy",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="orange",
       cex=9) 

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width


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

aa=SCC_COMP_MORE_75

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=".",
       col="navy",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="orange",
       cex=9) 

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width



#############################    /


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

aa=SCC_COMP_ADS_AGE

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=".",
       col="navy",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="orange",
       cex=9) 

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width


##################/          



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

aa=SCC_COMP_CNS_AGE

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=".",
       col="navy",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="orange",
       cex=9) 

# Type: p,l,b,c,o,s,S,h,n
# pch = 0:18 =:46
# col=
# bg= background
# lwd= line width

############################/




## prueba con un unico paciente contra todos los controles

singlePPT<-YcreatorAD(13)

# Response Variable (multiple):

Y_CN=SCC_matrix_CN
Y_AD=rbind(singlePPT,singlePPT)

# Parameters for SCC estimation:

d.est=5 # degree of spline for mean function 
d.band=2 # degree of spline for SCC
r=1 # smoothing parameter
lambda=10^{seq(-6,3,0.5)} # penalty parameters
alpha.grid=c(0.01,0.005,0.001) # vector of confidence levels


SCC_COMP_single=scc.image(Ya=Y_AD,Yb=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                          V.est.a=V.est,Tr.est.a=Tr.est,
                          V.band.a=V.band,Tr.band.a=Tr.band,
                          penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


plot(SCC_COMP_single,
     breaks=c(0,2),
     col="turquoise",
     # breaks=seq(from=0,to=2,length.out = 65),
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Difference between 1 random AD and CN group",
     col.sub="red")


aa=SCC_COMP_single

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan

image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points(points.P,
       type="p",
       pch=".",
       col="orange",
       cex=9)  

points(points.N,
       type="p",
       pch=".",
       col="navy",
       cex=9) 


##################################################/


# Percentage of locations where the true difference between mean functions value fall within the SCC

mfit$Yhat <- rbind(SCC_COMP_1$mu.hat.a,SCC_COMP_1$mu.hat.b)

xx <- Z[,1]
uu <- unique(xx)
n1 <- length(uu)
yy <- Z[,2]
vv <- unique(yy)
n2 <- length(vv)
npix <- (as.numeric(Z[nrow(Z),1]))*(as.numeric(Z[nrow(Z),2]));

# Generate true mean function
if(mu.func==1){
  beta.true <- 20*((xx-0.5)^2+(yy-0.5)^2)
}else if(mu.func==2){
  beta.true <- 5*exp(-15*((xx-0.5)^2+(yy-0.5)^2))+0.5
}else if(mu.func==3){
  beta.true <- 3.2*(-xx^3+yy^3)+2.4
}else if(mu.func==4){
  # beta.true <- -10*sin(7*pi*(xx+0.09))+10*sin(7*pi*(yy-0.14))+2.5
  beta.true <- -10*sin(5*pi*(xx+0.22))+10*sin(5*pi*(yy-0.18))+2.8
}
ind.outside <- setdiff(1:npix,ind.inside)
beta.true[ind.outside] <- NA
beta.true <- as.matrix(beta.true)
beta.diff=beta.true[,2]-beta.true[,1]


apply(SCC_COMP_2$scc,3,FUN=function(scc){
  sum((scc[,1]<beta.diff[ind.inside]) & (scc[,2]>beta.diff[ind.inside]))/length(ind.inside)
})


# Check if 0 fall within the SCC everywhere
apply(SCC_COMP_2$scc,3,FUN=function(scc){
  any((scc[,1]>0) | (scc[,2]<0))
})


# Simulation 1000 times to obtain empirical coverage rate and empirical power
nsim=100
coverage=lapply(1:nsim,FUN=function(iter){
  set.seed(iter)
  cat('Iteration No.',iter,'\t')
  ptm0=Sys.time()
  dat=data2g.image(na=na,nb=nb,Z=Z,ind.inside=ind.inside,mu1.func=mu1.func,
                   noise.type='Func',lam1=lam1,lam2=lam2,iter=iter,delta=delta)
  Ya=dat$Ya
  Yb=dat$Yb
  beta.diff=dat$beta.true[,2]-dat$beta.true[,1]
  
  out=scc.image(Ya=Ya,Yb=Yb,Z=Z,V.est.a=V.est.a,Tr.est.a=Tr.est.a,V.band.a=V.band.a,Tr.band.a=Tr.band.a,
                V.est.b=V.est.b,Tr.est.b=Tr.est.b,V.band.b=V.band.b,Tr.band.b=Tr.band.b,
                d.est=d.est,d.band=d.band,r=r,
                penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)
  cover.true=apply(out$scc,3,FUN=function(scc){
    sum((scc[,1]<beta.diff[ind.inside]) & (scc[,2]>beta.diff[ind.inside]))/length(ind.inside)
  })
  reject.null=apply(out$scc,3,FUN=function(scc){
    any((scc[,1]>0) | (scc[,2]<0))
  })
  ptm1=Sys.time()
  cat('Time: ',ptm1-ptm0,'\n')
  return(list(bw=out$bw,cover.true=cover.true,reject.null=reject.null))
})

# Empirical test power:
reject.all=sapply(coverage,'[[',3)
apply(reject.all,1,mean)

# Empirical coverage rate
cover.all=sapply(coverage,'[[',2)
apply(cover.all==1,1,mean)

# Mean bandwidth
bw.all=sapply(coverage,'[[',1)
apply(bw.all,1,mean)