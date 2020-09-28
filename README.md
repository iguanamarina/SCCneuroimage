# SCCneuroimage


## SUMMARY

In this GitHub repository you will find everything you need to replicate results obtained in: "Arias, J. A., Cadarso Suarez, C., & Fernandez Aguiar, P. (**Under review**). *Simultaneous Confidence Corridors in neuroimage data analysis: applications for Alzheimer's Disease diagnosis*" including:


* Positron Emission Tomography data extracted from the Alzheimer's Disease Neuroimaging Initiative database.

* Matlab code for the pre-processing of these images.

* R scripts for the ellaboration of SCCs for neuroimaging data.


## KEY STEPS:


1. PET images pre-processing
  + Realignement
  + Unwrapping
  + Corregistration
  + Normalization
  + Masking
  
2. Import PET images into R
  + Import one (1) individual
  + Loop for all participants
  + Create a complete and clean database (vars: PPT, group, sex, age, z, x, y, PET)
  + Mean average normalization
  
3. Calculate Simultaneous Confidence Corridors (SCC)
  + List of PPTs according to their group/sex/age or a combination
  + Convert data to a Functional Data setup (SCC matrices)
  + Extract contours for PET data
  + Compute triangulations over these contours
  + Construct SCCs for one-sample case
  + Construct SCCs for the difference between estimated mean functions
  
4. Visualization of regions suffering AD-induced neural loss
  + Visualize SCC for the difference between two groups
  + Overlay points falling above or below estimated confidence corridors
  + Exploration of results
