for i = 1:length(aa)

% aa=dir('*.nii')

% DATOS PET:

PET_matrix = niftiread(strcat('C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Niftis Simulacion/',aa(i).name)); 
% Guarda datos de la imagen PET en matriz "PET_matrix" (91x109x91); formato 'single'

% DATOS MÁSCARA:

mask_matrix = niftiread('w_mask_ICV.nii'); 
% Guarda datos de la  máscara en una matriz "mask_matrix" (91x109x91) y formato 'uint8'
mask_matrix_single = single(mask_matrix); 
% Casting para pasar los datos de mask_matrix a fomato 'single'

% DATOS HEADER IMPLICITOS:

info_PET = niftiinfo(strcat('C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/Niftis Simulacion/',aa(i).name)); 
% Read header info
info_PET_masked=info_PET;
% Así forzaremos a las nuevas imágenes enmascaradas a tener las propiedades que nos interesan
info_PET_masked.Filename=strcat('C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/masked/','masked_',aa(i).name); 
% Modifico del struct el nombre del archivo, todo lo demás lo quiero igual

% MASKING:

PET_matrix_masked = PET_matrix.*mask_matrix_single;
% Hago el enmascaramiento per se

% WRITE FILES:

niftiwrite(PET_matrix_masked, strcat('C:/Users/Juan A. Arias/Desktop/Simulaciones PET CHUS/masked/','masked_',aa(i).name),info_PET);
% Guardo mi PET enmascarado

disp(i)

end
