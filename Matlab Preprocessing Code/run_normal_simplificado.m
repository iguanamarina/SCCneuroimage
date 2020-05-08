for i = 1:length(aa)-2
    
    
    ppath1 = strcat('/home/alexis/Escritorio/ADNI-Juan/MRI_bien/',aa(i+2).name,'/T1.img');
    
        path_MRI = ppath1;
    
        
    ppath2 = strcat('/home/alexis/Escritorio/ADNI-Juan/PET_bien/',aa(i+2).name,'/P1.img');
    
        path_rPET = ppath2;
        
        a_PET = dir(ppath2);
    
   
    ppathw = strcat('/home/alexis/Escritorio/ADNI-Juan/PET_bien/',aa(i+2).name,'/wP1.img');
        
        a_wPET = dir(ppathw);
    
        
    
    if length(a_PET) == 1 && isempty(a_wPET)
    
       
        a = normal(path_MRI,path_rPET);
    
    end
    
    
    
    disp(i)    
     
    
end
