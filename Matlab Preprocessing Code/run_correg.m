for i = 1:length(aa)-2
    
    
    ppath = strcat('./CN_MRI_Base/',aa(i+2).name,'/**/*img');
    
    a_MRI = dir(ppath);
    
    path_MRI = strcat(a_MRI.folder,'/',a_MRI.name);
    
    
    
    ppath = strcat('./CN_FDG-PET_Base/',aa(i+2).name,'/**/*img');
    
    a_PET = dir(ppath);
    
    path_PET = strcat(a_PET.folder,'/',a_PET.name);
    
    
    
    if length(a_PET) == 1
    
        a = correg(path_MRI,path_PET);
        
    end
    
    disp(i)
    
    
end
