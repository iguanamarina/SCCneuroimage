#######################################################################\
#   THE GOALS OF THIS SCRIPT ARE:
#     - CREATE DATABASE
#     - 
#     - TO GET A LIST OF THE PPT'S NAMES        
#     - CREATE A NEW DATABASE FOR SCCs
#     - INCLUDE NEW TRIANGULATION PARAMETERS
#     - TEST SCC'S WITH DEFINITIVE DATA
#######################################################################\


## REALIGNEMENT, WRAPPING, CORREGISTRARION, NORMALIZATION, AND MASKING DONE IN 'MATLAB':

## INSTALL NEURO-PACKAGES: ##

install.packages(C("mgcv","gamair","oro.nifti","memsic"))

library(mgcv);library(gamair);library(oro.nifti);library(memisc)


##########################################################\
###             f.clean  (PPT,z,x,y,pet)               ###
##########################################################\

## Set as working directory -> directory with .hdr/.img NIfTI files and Demographics

setwd("~/MEGA/PhD/4. Wood Example (GAM's)/Brain Imaging Example (Simon Wood)/PETimg_masked")

f.clean <- function(name) { #### f.clean is meant for CLEANING ONE SINGLE PPT DATA
  
  # Read the NIFTI image, transform it to dataframe, preserve slice Z and organize the table
  
  ## Load Data
  
  file <- readNIfTI(fname = name, verbose = FALSE, warn = -1, reorient = TRUE, call = NULL, read_data = TRUE)
  namex <- as.character(name)
  n = img_data(file)
  n = to.data.frame(n)
  
  ## Prepare data.frame base where surther data from the loop will be integrated
  
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
    
    temp0 = data.frame(z,x,y,pet) # tmeporal dataframe
    temp1 <- print(temp0) # is this necessary?
    dataframe <- rbind(dataframe,temp1) # sum new data with previous data
    
  }
  
  #Demographics: PPT, group (AD/CN), sex, age.
  
  demo <- read.csv2("Demographics.csv")
  demo <- demo[demo$PPT==namex,]
  
  PPT <- rep(demo$PPT,length.out=7505)
  group <-rep(demo$Group,length.out=7505)
  sex <-rep(demo$Sex,length.out=7505)
  age <-rep(demo$Age,length.out=7505)
  
  temp2 <- data.frame(PPT,group,sex,age)
  dataframe <- cbind(temp2,dataframe)
  
  print(dataframe) # Necessary for asigning an object name
  
  
}

#Example(s):
"003_S_1059" = f.clean("003_S_1059")
head(`003_S_1059`)  # Some values are Zeros due to the masking process



##########################################################
###               COMPLETE   DATABASE                  ###
##########################################################

files <- list.files(path="~/MEGA/PhD/4. Wood Example (GAM's)/Brain Imaging Example (Simon Wood)/PETimg_masked", 
                    pattern="*.img", full.names=F, recursive=FALSE) # list of files

files <- gsub(files, pattern=".img$", replacement="") # remove file extension .img

database <- data.frame(PPT=integer(),group=integer(),sex=integer(),age=integer(),z=integer(),x=integer(),y=integer(),pet=integer())
#create data.frame to include data

for (i in 1:length(files)) { #loop to include every PPT in the dataframe
  
  temporal <- f.clean(files[i])
  database <- rbind(database,temporal)
  
}

nrow(database[database$pet<0,]) # Existen valores negativos

database$pet[database$pet<0]<- NaN  # Se convierten tanto los negativos como los zeros en NaN !!! CAREFUL HERE por que no despues del mean normalization=?

database$pet[database$pet==0] <- NaN 

write.csv2(database,file="Database.csv",sep=";",na="NA") #export as .csv IF DESIRABLE



##########################################################
###            MEAN AVERAGE NORMALIZATION              ###
##########################################################


mean.data <- data.frame(pet=integer(),pet_normal=integer())
# Data is the data.framecreated for masked and mean averaged data

memory.size(max = TRUE)

for (i in 1:length(files)) {
  
  temp <- database[database$PPT==files[i],]
  temp <- temp[,8]
  mean <- mean(as.numeric(temp),na.rm=T)
  pet_normal <- as.data.frame(temp/mean)
  colnames(pet_normal)<-c("pet_normal")
  temp <- cbind(temp,pet_normal)
  show(temp)    
  mean.data <- rbind(mean.data,temp)
  
}

mean.data <- mean.data[,2]
colnames(mean.data)<-c("pet")
database <- cbind(database,mean.data)
database$pet<-database$mean.data
database<-database[,-9]

write.table(database,file="Database.csv",sep=";",na="NA") #export as .csv IF DESIRABLE

save(database,file="data.Rda")


# Now we have the database with our complete processed stuff and we can move on with this towards SCC's


##########################################################
###            HERE STARTS PROPERLY SCC CODE           ###
##########################################################

# load("F:/MEGAsync/PhD/data.Rda")

# Install packages from GitHub: #

install.packages("devtools");library(devtools) 
install.packages("remotes");library(remotes) 
install.packages("readr")
install.packages("imager")
install.packages("itsadug") 
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE) # Forces installation to finish even with warning messages (EXTREMELY NECESSARY)
remotes::install_github("funstatpackages/BPST")
remotes::install_github("funstatpackages/Triangulation")
remotes::install_github("funstatpackages/ImageSCC")


library(BPST);library(Triangulation);library(ImageSCC);library(readr);library(imager);library(itsadug);library(fields)


#########################################################################//

# Preliminary modifications and visualization:

Data<-database
attach(Data)
head(Data)
str(Data)


#########################################################################//


# LIST PPT's:

names<-as.vector(Data[,1])
names<-names[duplicated(names)!=TRUE]
((length(names))*592895)==nrow(Data)  # Every PPT has 592.895 measures. PPT*measures=nrow(Data). Should return a TRUE

names_CN<-as.vector(Data[,1])[Data$group=="CN"]
names_CN<-names_CN[duplicated(names_CN)!=TRUE]
nCN=length(names_CN)

names_AD<-as.vector(Data[,1])[Data$group=="AD"]
names_AD<-names_AD[duplicated(names_AD)!=TRUE]
nAD=length(names_AD)


## FUNCTION TO CREATE Y's FOR SCC's in the shape of a dataframe with only one row

YcreatorCN <- function(i){
  Y <- subset(Data, Data$PPT==names_CN[i] & Data$z==30) # Y = Choose by PPT and Z always 30, group=change depending on the SCC you need
  Y <- Y[1:7505,8] # Make sure to keep only 7505 (problematic before) and column 9 (responses)
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

YcreatorAD <- function(i){
  Y <- subset(Data, Data$PPT==names_AD[i] & Data$z==30) # Y = Choose by PPT and Z always 30, group=change depending on the SCC you need
  Y <- Y[1:7505,8] # Make sure to keep only 7505 (problematic before) and column 9 (responses)
  Y <- as.matrix(Y)
  Y=t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}


# CN LOOP: now every CN-PPT will be one row, as ImageSCC requires

SCC_matrix_CN <- matrix(ncol = 7505,nrow=0)
for (i in 1:nCN){
  temp<-YcreatorCN(i)
  SCC_matrix_CN<-rbind(SCC_matrix_CN,temp)    
}


# AD LOOP: now every AD-PPT will be one row, as ImageSCC requires

SCC_matrix_AD <- matrix(ncol = 7505,nrow=0)
for (i in 1:nAD){
  temp<-YcreatorAD(i)
  SCC_matrix_AD<-rbind(SCC_matrix_AD,temp)    
}



#########################################################################
#     NOW THAT WE ALREADY HAVE THIS LIST AMD THE SCRIPT TO GET MORE     #
#     WE TRY TO REPRODUCE SCC CORRIDORS FOR MULTIPLE PPT'S              #
#########################################################################


# I'm using my previous file "sample" because it is the fastest way to get the coords. I need,
# but any other way is valid

sample <- read_csv("sample.csv", col_types = cols(pet = col_number(), x = col_integer(), y = col_integer()))
sample <- as.matrix(sample)
colnames(sample) <- c("X","Y","pet")
Z <- sample[,1:2]
Z


# Triangulation Parameters: Still not the good ones

Brain.V3  # This is just for this test, mi Brain images do not correspond with the shapes provided by the package so I'll develop a new one with 'Triangulation' package.
plot(Brain.V3,main="Vertices to Use") # The orientation of the vertices does not correspond with my data (see attached image)

Brain<-cbind(Brain.V3[,2],Brain.V3[,1]) # in order to adjust vertices I do a trasposition
plot(Brain,main="New Brain Vertices")

V.est=as.matrix(Brain)
V.est=cbind((V.est[,1]*79),(V.est[,2]*95)) 
# in order to keep the original proportions of V.est while adapting it to my 79*95 dimensional data
plot(V.est,main = "New Proportional Brain")  # Now oriented to fit my data
Tr.est=as.matrix(Brain.Tr3)  

V.band=as.matrix(Brain)
V.band=cbind((V.band[,1]*79),(V.band[,2]*95)) # Applying the same logic
Tr.band=as.matrix(Brain.Tr3)  


# Response Variable (multiple):

Y_CN=SCC_matrix_CN
Y_AD=SCC_matrix_AD

# Parameters for SCC estimation:

d.est=5 # degree of spline for mean function 
d.band=2 # degree of spline for SCC
r=1 # smoothing parameter
lambda=10^{seq(-6,3,0.5)} # penalty parameters
alpha.grid=c(0.10,0.05,0.01) # vector of confidence levels


# Run one sample SCC construction one time:


SCC_CN_1=scc.image(Ya=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                   V.est.a=V.est,Tr.est.a=Tr.est,
                   V.band.a=V.band,Tr.band.a=Tr.band,
                   penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)

plot(SCC_CN_1,breaks=seq(from=0,to=2,length.out = 65),
     col=,
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Control Group",
     col.sub="red")



SCC_AD_1=scc.image(Ya=Y_AD,Z=Z,d.est=d.est,d.band=d.band,r=r,
                   V.est.a=V.est,Tr.est.a=Tr.est,
                   V.band.a=V.band,Tr.band.a=Tr.band,
                   penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)

plot(SCC_AD_1,breaks=seq(from=0,to=2,length.out = 65),
     col=,
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Alzheimer Group",
     col.sub="red")


####################################################
#   AND NOW FOR A COMPARATION BETWEEN CN AND AD
#################################################### 


SCC_COMP_1=scc.image(Ya=Y_AD,Yb=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                     V.est.a=V.est,Tr.est.a=Tr.est,
                     V.band.a=V.band,Tr.band.a=Tr.band,
                     penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    
# With Yb included: Two-group SCC is constructed for the difference between the mean functions of Yb and Ya

plot(SCC_COMP_1,breaks=seq(from=0,to=2,length.out = 65),
     col=,
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Difference between the estimated mean functions between Controls and Alzheimer's",
     col.sub="red")




######TEST#####


SCC_COMP_2=scc.image(Ya=Y_AD,Yb=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                     V.est.a=V.est,Tr.est.a=Tr.est,
                     V.band.a=V.band,Tr.band.a=Tr.band,
                     penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    
# With Yb included: Two-group SCC is constructed for the difference between the mean functions of Yb and Ya

plot(SCC_sample,breaks=seq(from=0,to=2,length.out = 65),
     col=,
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Difference between the estimated mean functions between Controls and Alzheimer's",
     col.sub="red")

Z.band <- matrix(SCC_COMP_1$Z.band,ncol=2)
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]);
n1 <- length(z1); n2 <- length(z2)
scc <- matrix(NA,n1*n2,2)
ind.inside.band <- SCC_COMP_1$ind.inside.cover
scc[ind.inside.band,] <- SCC_COMP_1$scc[,,2]
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE))
scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1)
scc.l.mtx[scc.l.mtx<0]=NA
image.plot(z2,z1,scc.l.mtx, zlim = scc.limit)

points<-which(scc.l.mtx>0,arr.ind=TRUE)

points(points)


#############################    


## prueba con un unico paciente contra todos los controles

singlePPT<-YcreatorAD(45)

# Response Variable (multiple):

Y_CN=SCC_matrix_CN
Y_AD=rbind(singlePPT,singlePPT)

# Parameters for SCC estimation:

d.est=5 # degree of spline for mean function 
d.band=2 # degree of spline for SCC
r=1 # smoothing parameter
lambda=10^{seq(-6,3,0.5)} # penalty parameters
alpha.grid=c(0.10,0.05,0.01) # vector of confidence levels


SCC_COMP_single=scc.image(Ya=Y_AD,Yb=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                          V.est.a=V.est,Tr.est.a=Tr.est,
                          V.band.a=V.band,Tr.band.a=Tr.band,
                          penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)    


plot(SCC_COMP_single,breaks=seq(from=0,to=2,length.out = 65),
     col=,
     xlab="Longitudinal (1-95)",
     ylab="Transversal (1-79)",
     sub="Difference between 1 random AD and CN group",
     col.sub="red")



Z.band <- matrix(SCC_COMP_single$Z.band,ncol=2)
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]);
n1 <- length(z1); n2 <- length(z2)
scc <- matrix(NA,n1*n2,2)
ind.inside.band <- SCC_COMP_single$ind.inside.cover
scc[ind.inside.band,] <- SCC_COMP_single$scc[,,2]
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE))
scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1)
scc.l.mtx[scc.l.mtx<0]=NA
image.plot(z2,z1,scc.l.mtx, zlim = scc.limit)

points<-which(scc.l.mtx>0,arr.ind=TRUE)

points(points, pch = 3)


##################################################






# Percentage of locations where the true difference between mean functions value fall within the SCC

mfit$Yhat <- rbind(out$mu.hat.a,out$mu.hat.b)

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