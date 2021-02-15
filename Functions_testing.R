new.f.clean <- function(name) { #### f.clean is meant for CLEANING ONE SINGLE-PPT DATA, then we loop it
  
  # Read NIFTI image, transform it to dataframe, preserve slice Z and organize the table
  
  ## Load Data
  
  file <- readNIfTI(fname = name, verbose = FALSE, warn = -1, reorient = TRUE, call = NULL, read_data = TRUE)
  namex <- as.character(name)
  n = to.data.frame(img_data(file))
  
  ## Prepare data.frame base where further data from the loop will be integrated
  
  dataframe <- data.frame(z=integer(),x=integer(),y=integer(),pet=integer()) 
  
  # Loop for 79 slices of Z in the NiFtI image -> move to dataframe
  
  for (i in 1:79) {
    
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

example2 = new.f.clean("003_S_1059")
head(example2)  # Some values are Zeros due to the masking process