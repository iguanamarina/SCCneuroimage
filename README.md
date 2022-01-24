# SCCneuroimage

## SUMMARY

In this GitHub repository you will find everything you need to replicate results obtained in: "Arias, J. A., Cadarso Suarez, C., & Fernandez Aguiar, P. (**Under review**). *Simultaneous Confidence Corridors in neuroimage data analysis: applications for Alzheimer's Disease diagnosis*" including:

-   Positron Emission Tomography (PET) data extracted from the Alzheimer's Disease Neuroimaging Initiative database.

-   PET data and ROIs for simulation study.

-   Matlab code for the pre-processing of these images.

-   R scripts for the ellaboration of SCCs for neuroimaging data.

-   R scripts for the evaluation of SCCs compared to classical SPM for these neuroimaging datasets.

## KEY STEPS:

### PET images pre-processing

-   Realignement
-   Unwrapping
-   Corregistration
-   Normalization
-   Masking

### Import PET images into R

-   Import one (1) individual
-   Loop for all participants
-   Create a complete and clean database (vars: PPT, group, sex, age, z, x, y, PET)
-   Mean average normalization

### Calculate Simultaneous Confidence Corridors (SCC)

-   List of PPTs according to their group/sex/age or a combination
-   Convert data to a Functional Data setup (SCC matrices)
-   Extract contours for PET data
-   Compute triangulations over these contours
-   Construct SCCs for one-sample case
-   Construct SCCs for the difference between estimated mean functions

### Visualization of regions suffering AD-induced neural loss

-   Visualize SCC for the difference between two groups
-   Overlay points falling above or below estimated confidence corridors
-   Exploration of results

### Evaluation of obtained predictions

## CHEATSHEET:

1.  Open "`MASTERSCRIPT (for simulations).R`"

2.  Press `Enter` in everything, un-comment lines 36 & 37 if necessary, modify `param.z` if necessary.

3.  Create folders as requested (yet to be improved) and copy/paste SCC results if you have them already (as it is with my case).

4.  Also create SPM folder with `binary.nii` files. This has to be carried out manually.

5.  
