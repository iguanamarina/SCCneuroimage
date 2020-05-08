for i = 1:length(aa)-2
    
    
    path_MRI = strcat('./MRI_Bien/',aa(i+2).name,'/T1.img');
    
   
    
    path_rPET = strcat('./PET_Bien/',aa(i+2).name,'/P1.img');
    
  
    
    a_rPET = dir(ppath);
    
    path_rPET = strcat(a_rPET.folder,'/',a_rPET(2).name); 
    %always choose second file: coregistered PET
    
    
    if length(a_rPET) > 1
        
    a = normal(path_MRI,path_rPET);
       
    %Si hay más de 1 file es que hay corregister: run "normal"
    
    end
    
    disp(i)
    
    
end

%%%%%%%%%%% NO USAR ESTE %%%%%%%%%%%%%%%