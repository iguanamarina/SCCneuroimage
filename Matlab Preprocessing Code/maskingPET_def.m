for i = 1:length(aa)

% aa=dir('*.img')

% Input/Output files

path_PET = fopen(strcat('C:/Users/Juan A. Arias/Desktop/PETimg/',aa(i).name),'r');
path_MASK = fopen('C:/Users/Juan A. Arias/Desktop/mask.img','r');
path_EXPORT = fopen(strcat('C:/Users/Juan A. Arias/Desktop/PETmasked/',aa(i).name),'w');

% Input vectors
PET = fread(path_PET,'float','ieee-be');
MASK = fread(path_MASK,'char');

% Masking
MPET = MASK.*PET/100;

% Save output file
fwrite(path_EXPORT,MPET,'float','ieee-be');

disp(i)

end

