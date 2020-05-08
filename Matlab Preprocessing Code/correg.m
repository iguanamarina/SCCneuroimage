function a = correg(path_MRI,path_PET)

path_MRI_1 = strcat(path_MRI,',1');

path_PET_1 = strcat(path_PET,',1');

spm('defaults','fmri');
spm_jobman('initcfg');


matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {path_MRI_1};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {path_PET_1};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

a = 0;

spm_jobman('run',matlabbatch);