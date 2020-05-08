for i = 1:53

ppath_hdr = strcat('/Usuario/Desktop/Alzheimer Disease vs Control/1.AD/MRI/',aa(i+2).name,'/**/*.hdr');
ppath_img = strcat('/Usuario/Desktop/Alzheimer Disease vs Control/1.AD/MRI/',aa(i+2).name,'/**/*.img');

ppath_hdr = dir(ppath_hdr);
ppath_img = dir(ppath_img);

mkdir(strcat('/Usuario/Desktop/Alzheimer Disease vs Control/4.Analysis/AD/MRI/',aa(i+2).name));

if length(ppath_hdr) == 1
    
copyfile(strcat(ppath_hdr.folder,'/',ppath_hdr.name),strcat('/Usuario/Desktop/Alzheimer Disease vs Control/4.Analysis/AD/MRI/',aa(i+2).name,'/R1.hdr'));
copyfile(strcat(ppath_img.folder,'/',ppath_img.name),strcat('/Usuario/Desktop/Alzheimer Disease vs Control/4.Analysis/AD/MRI/',aa(i+2).name,'/R1.img'));

end

end
