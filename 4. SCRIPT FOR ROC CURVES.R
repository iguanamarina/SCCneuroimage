############################## ################### ############## ### 
##
## Script name: SCRIPT FOR ROC CURVES
##
## Purpose of script: To compute sort-of-ROC curves (sensibility vs especificity) 
## for SCC and SPM and get AUC in z=30, roiAD, hypo=0'4. This might not be correct, 
## anyway here is our attempt. 
##
## Date Created: 2021-11-26
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://messy-dataset.xyz
##
## Notes: This script is self-contained and nothing should be necessary but for
## pressing enter again and again (after hyperparameter setup). However, results
## may differ with slight changes and, still, Carmen didn't give an OK to this. 
##   
############################## ################### ############## ### 


####
# PREAMBLE: ----
####


#* Set working directory: ----

setwd("~/GitHub/SCCneuroimage") 

#* Tune Options: ----
options(scipen = 6, digits = 6) # View outputs in non-scientific notation
memory.limit(30000000)     # This is needed on some PCs to increase memory allowance

#* Load up packages: ---- 

library(gamair);library(oro.nifti);library(memisc);library(devtools);library(remotes);library(readr);library(imager);library(itsadug);library(ggplot2);library(contoureR);library(fields);library(BPST);library(Triangulation);library(ImageSCC);library(tidyr); library(dplyr);library(stringr); library(threadr);library(memisc); library(refreg); library(MESS); library(sfsmisc); library(Bolstad2)

#* Load up functions: ----

load("~/GitHub/SCCneuroimage/Functions/f.clean.RData") # Function for NiFTi -> df
load("~/GitHub/SCCneuroimage/Functions/my_points.RData") # Function for getting relevant points from SCC
load(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z), "/", "contour", as.numeric(param.z), ".RData")) # Contour at z
load(paste0("~/GitHub/SCCneuroimage/z", as.numeric(param.z), "/", "SCC_CN.RData")) # Data for controls
load(paste0("~/GitHub/SCCneuroimage/z30/pseudoROC/SCC_MULTIPLE_ALPHA_roiAD_4.RData")) # Data for pathological

#* Hyper-parameters: ----

param.z = 30
region <- "roiAD"
roi <- 4
number <- paste0("C", 1:25)


####  
# PART 1: COMPUTE SCC ----------
####


#* Parameters for SCC estimation ----

Brain.V <- VT[[1]]; Brain.Tr <- VT[[2]]
V.est = as.matrix(Brain.V); Tr.est = as.matrix(Brain.Tr)
V.band = as.matrix(Brain.V); Tr.band = as.matrix(Brain.Tr)

#* Mean Average Normalization: ----
      
  for (k in 1:nrow(SCC_CN)) {
    
    temp <- SCC_CN[k, ]
    mean <- mean(as.numeric(temp), na.rm = T)
    SCC_CN[k, ] <- (temp/mean)
  }
  
  for (k in 1:nrow(SCC_matrix)) {
    
    temp <- SCC_matrix[k, ]
    mean <- mean(as.numeric(temp), na.rm = T)
    SCC_matrix[k, ] <- (temp/mean)
  }

#* Wang et al's parameters for SCCs: ----
    
  d.est = 5 # degree of spline for mean function  5
  d.band = 2 # degree of spline for SCC  2
  r = 1 # smoothing parameter  1
  lambda = 10^{seq(-6,3,0.5)} # penalty parameters
  alpha.grid = c(0, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 0.999, 0.9999, 1) 
  # vector of confidence levels HARDCORE
  x <- rep(1:91, each = 109, length.out = 9919) 
  y <- rep(1:109, length.out = 9919)
  Z <- cbind(as.matrix(x),as.matrix(y)); head(Z)


#* Build SCC: ----

setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z30/SCC30results_multiple_alpha")

SCC_MULTIPLE_ALPHA = scc.image(Ya = SCC_matrix, Yb = SCC_CN, Z = Z, 
                               d.est = d.est, d.band = d.band, r = r,
                               V.est.a = V.est, Tr.est.a = Tr.est,
                               V.band.a = V.band, Tr.band.a = Tr.band,
                               penalty = TRUE, lambda = lambda, 
                               alpha.grid = alpha.grid, 
                               adjust.sigma = TRUE)    

#* Save/Load SCC: ----
    
#save(SCC_MULTIPLE_ALPHA, file = paste0("SCC_MULTIPLE_ALPHA_", region, "_", roi,".RData"))
load(paste0("~/GitHub/SCCneuroimage/z30/pseudoROC/SCC_MULTIPLE_ALPHA_roiAD_4.RData")) # Data for pathological
     
  
####  
# PART 2: THEORETICAL POINTS----------
####  


region <- "wroiAD" # Because of nomenclature reasons
ROI_data <- list() # Cargar las databases

for (i in 1:length(region)) {
  for (k in 1:length(number)) {
    
  setwd("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/roisNormalizadas/tables")
  ROI_data[[paste0(as.character(region[i]), "_", as.character(number[k]))]] <- readRDS(
    paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/roisNormalizadas/tables/ROItable_", 
           region[i], "_", number[k], ".RDS"))
  }
}  

#* T_points: Theoretical points at Z ----

T_points <- list() 

for (i in 1:length(region)) {
  
  for (j in 1:length(number)) {
    
    T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]] <- ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]][ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]]$z == param.z & ROI_data[[paste0(as.character(region[i]), "_", as.character(number[j]))]]$pet == 1, 3:4]
    
    T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]] <- unite(as.data.frame(T_points[[paste0(as.character(region[i]), "_", as.character(number[j]))]]), newcol, c(y,x), remove = T)
        
    }
  }
      
rm(ROI_data) # No longer necessary


####  
# PART 3: SENS/ESP FOR SCC----------
####  


#* H_points: Hypotetical points for SCCs ----

# setwd(paste0("C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/SCC matrix/z", as.numeric(param.z) ,"/SCC", as.numeric(param.z), "results", "_multiple_alpha"))
region <- "roiAD" # again, nomenclature reasons in the neuroimage dataset
alpha <- SCC_MULTIPLE_ALPHA[["alpha"]]
  
H_points <- list()
  
for (i in 1:length(alpha)) {

  H_points[[as.character(alpha[i])]] <- my_points(SCC_MULTIPLE_ALPHA, i)[[1]] 
  H_points[[as.character(alpha[i])]] <- unite(as.data.frame(H_points[[i]]), newcol, c(row, col), remove = T)
  
}

#* Total Coords: ----

x <- rep(1:91, each = 109, length.out = 9919); y <- rep(1:109, length.out = 9919)
total.coords <- data.frame(y, x) ###!!!!!
total.coords <- unite(as.data.frame(total.coords), newcol, c(y, x), remove = T)


#* SCC_sens_esp Loop: ----

SCC_sens_esp <- data.frame(alpha = integer(), group = integer(), sens = integer(), esp = integer()) # df final
region <- "wroiAD" # again, nomenclature
  
for (i in 1:length(alpha)) {
  
  for (j in 1:length(number)) {
  
    inters <- inner_join(H_points[[as.character(alpha[i])]],
                         T_points[[paste0(as.character(region), "_", as.character(number[j]))]])
    
    sensibilitySCC <- nrow(inters)/nrow(T_points[[paste0(as.character(region), "_", as.character(number[j]))]])*100 
    # p(correctly identify hypoactive pixel)
    
    true_neg <- anti_join(total.coords, T_points[[paste0(as.character(region), "_", as.character(number[j]))]])
    hypo_neg <- anti_join(total.coords, H_points[[as.character(alpha[i])]])
    anti_inters <- inner_join(true_neg, hypo_neg)
    
    specificitySCC <- nrow(anti_inters)/nrow(true_neg)*100
    # p(correctly identify healthy pixel)
  
    temp <- data.frame(alpha = alpha[i], group = number[j], sens = sensibilitySCC, esp = specificitySCC)  
    SCC_sens_esp <- rbind(SCC_sens_esp, temp)
    
  }
}

setwd("~/GitHub/SCCneuroimage/z30/pseudoROC")
saveRDS(SCC_sens_esp, file = "SCC_sens_esp.RDS")
tableSCC <- readRDS(paste0("~/GitHub/SCCneuroimage/z30/pseudoROC", "/SCC_sens_esp.RDS"))

tableSCC$alpha <- factor(tableSCC$alpha,      
                         levels = alpha)
tableSCC$inv_esp <- 100 - tableSCC$esp

View(tableSCC)


####  
# PART 4: SENS/ESP FOR SPM----------
####  

#* SPM_sens_esp Loop: ----

SPM_sens_esp <- data.frame(alpha = integer(), group = integer(), sens = integer(), esp = integer())
region <- "wroiAD"
  
for (i in 1:length(alpha)) {

  if (i == 2) { # En i=2 el alpha es 1e-04 y eso da problemas con setwd, de ahí este IF-ELSE
  setwd(paste0("~/GitHub/SCCneuroimage/z30/pseudoROC", "/SPM/", formatC(alpha[2], format = "f", digits = 4)))
  } else {
  setwd(paste0("~/GitHub/SCCneuroimage/z30/pseudoROC", "/SPM/", as.character(alpha[i])))
  }
  
  binary <- f.clean("binary.nii")
  H_points_SPM <- binary[binary$z == as.numeric(param.z) & binary$pet == 1, 2:3]
  H_points_SPM <- unite(as.data.frame(H_points_SPM), newcol, c(y, x), remove = T)
  
  sens_esp <- data.frame(region = integer(), group = integer(), sens = integer(), esp = integer())
  
  for (j in 1:length(number)) {
  
    inters <- inner_join(H_points_SPM,
                         T_points[[paste0(as.character(region), "_", as.character(number[j]))]])
  
    sensibilitySPM <- nrow(inters)/nrow(T_points[[paste0(as.character(region), "_", as.character(number[j]))]])*100 
    # p(correctly identify hypoactive pixel)
    
    true_neg <- anti_join(total.coords, T_points[[paste0(as.character(region), "_", as.character(number[j]))]])
    hypo_neg <- anti_join(total.coords, H_points_SPM)
    anti_inters <- inner_join(true_neg, hypo_neg)
    
    specificitySPM <- nrow(anti_inters)/nrow(true_neg)*100
    # p(correctly identify healthy pixel)
    
    temp <- data.frame(alpha = alpha[i], group = number[j], sens = sensibilitySPM, esp = specificitySPM)  
    sens_esp <- rbind(sens_esp, temp)

  }

  SPM_sens_esp <- rbind(SPM_sens_esp, sens_esp)

}

setwd("~/GitHub/SCCneuroimage/z30/pseudoROC")
saveRDS(SPM_sens_esp, file = "SPM_sens_esp.RDS")
tableSPM <- readRDS(paste0("~/GitHub/SCCneuroimage/z30/pseudoROC/", "SPM_sens_esp.RDS"))

tableSPM$alpha <- factor(tableSPM$alpha,      
                         levels = alpha)
tableSPM$inv_esp <- 100 - tableSPM$esp

View(tableSPM)


####  
# PART 5: ROC CURVES AND AUC CALCULATION----------
####  

#* For SCC: ----

dat = tableSCC
dat$sens = dat$sens/100
dat$inv_esp = dat$inv_esp/100

sen_i = c(0, dat$sens, 1) # Add for extrapolation
esp_i = c(0, dat$inv_esp, 1) # Add for extrapolation

dat_i = data.frame(sen_i, esp_i)

  ## Infer data out of range:
  
  m1 = ACE(y = "sen_i", predictor = "~s(esp_i)",
       restriction = "correlation",
       data = dat_i)
  
  nw1 = data.frame(esp_i = seq(0,1,0.01))
  pred = tanh(predict(m1$fit, newdata = nw1))


#* For SPM: ----

dat2 = tableSPM
dat2$sens = dat2$sens/100
dat2$inv_esp = dat2$inv_esp/100

sen_i = c(0, dat2$sens, 1) # Add for extrapolation
esp_i = c(0, dat2$inv_esp, 1) # Add for extrapolation

dat_i = data.frame(sen_i, esp_i)

  ## Infer data out of range:
  
  m2 = ACE(y = "sen_i", predictor = "~s(esp_i)",
       restriction = "correlation",
       data = dat_i)
  
  nw2 = data.frame(esp_i = seq(0,1,0.01))
  pred2 = tanh(predict(m2$fit, newdata = nw2))


#* ROC Plots: ----

  
#** For SCC:  ----
  
plot(nw1$esp_i, pred, type = "o", xlab = "1-Especificity", ylab = "Sensibility", main = "Estimated Sens vs Esp Graph for SCC")
  points(dat$sens ~ dat$inv_esp, col = as.factor(dat$alpha))
  abline(h = 1, v = 0)

#** For SPM:  ----
  
plot(nw2$esp_i, pred2, type = "o", xlab = "1-Especificity", ylab = "Sensibility", main = "Estimated Sens vs Esp Graph for SPM")
  points(dat2$sens ~ dat2$inv_esp, col = as.factor(dat$alpha))
  abline(h = 1, v = 0) 

#** Both plotted together:----

plot(nw2$esp_i, pred2, type = "l", xlab = "1-Especificity", ylab = "Sensibility", main = "Estimated Sens vs Esp Graph for SPM/SCC")
  lines(nw1$esp_i, pred, col = "black")  
  points(dat2$sens ~ dat2$inv_esp, col = "darkblue")
  points(dat$sens ~ dat$inv_esp, col = "red")
  abline(h = 1, v = 0)
  
  
#* AUCs: ----

df_aucSCC = data.frame(x = nw1$esp_i,
                      y = pred)
df_aucSPM = data.frame(x = nw2$esp_i,
                      y = pred2)

#** for SCC: ----

AUCscc <- data.frame(AUC = auc(df_aucSCC$x, df_aucSCC$y, type = 'spline'), # 0.920059
                     integrate.xy = integrate.xy(df_aucSCC$x, df_aucSCC$y), # 0.920058
                     sintegral = sintegral(df_aucSCC$x, df_aucSCC$y)$int) # 0.919923

#** for SPM: ----

AUCspm <- data.frame(AUC = auc(df_aucSPM$x, df_aucSPM$y, type = 'spline'), # 0.896913
                     integrate.xy = integrate.xy(df_aucSPM$x, df_aucSPM$y), # 0.896914
                     sintegral = sintegral(df_aucSPM$x, df_aucSPM$y)$int) # 0.89687

AUC <- rbind(AUCscc, AUCspm)
rownames(AUC) <- c("SCC", "SPM")
View(AUC)

