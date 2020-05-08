for i = 1:length(aa)-2
    
    
    path_MRI = strcat('./MRI_Bien/',aa(i+2).name,'/**/*img');
    
   
    path_rPET = strcat('./PET_Bien/',aa(i+2).name,'/**/*img');
    
  
    a = normal(path_MRI,path_rPET);
       
    disp(i)
    
    end
    
     
    
