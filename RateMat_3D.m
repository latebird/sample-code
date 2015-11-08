function [Rate_Mtx,Os_idx] = RateMat_3D(transitionRate,DecayRate,...
                             length, width, height,num_Os)
n             = length*width*height+2;
idx           = zeros(7,3);
Rate_Mtx      = zeros (n,n);
Os_idx        = zeros(1,num_Os);
flag          = 0;


for a=1:length
     for b=1:width
        for c=1:height
           idx(1,:)=[a,b,c]  ;
           idx(2,:)=[a-1,b,c];
           idx(3,:)=[a+1,b,c];
           idx(4,:)=[a,b-1,c];
           idx(5,:)=[a,b+1,c];
           idx(6,:)=[a,b,c-1];
           idx(7,:)=[a,b,c+1];
         
          % Modify the neighbor index base on perodic bandary condition    
                
           site    = (a-1)*width*height+(b-1)*height+c;
           counter = 0; 
           for i=2:7
               
               if 0.5<idx(i,1) && idx(i,1)<length+0.5&&0.5<idx(i,2)&&...
                  idx(i,2)<width+0.5&& 0.5<idx(i,3)&&idx(i,3)<height+0.5
                      neighor = (idx(i,1)-1)*width*height+(idx(i,2)-1)*height+idx(i,3);
                      Rate_Mtx(site,neighor) = transitionRate(1,1);
                      counter                = counter + 1;
               end
               
           end
           Rate_Mtx(site, site) = -counter*transitionRate(1,1)-DecayRate(1);
           Rate_Mtx(n-1,site)   = DecayRate(1);
        end
     end
end

  %-------------------
  % Select the Os sites    |
  %-------------------
  temp        = randperm(n);
  Os_idx     = temp(1:num_Os);
  clear temp;
  Sitename  = zeros(n,1);
  
  for i=1: num_Os
        Sitename(Os_idx(i)) = 1;
  end
 %---------------------------------------------------------
 %Loop used to modify the rate matrix due to the add of Os| 
 %---------------------------------------------------------

 for Os= 1 : num_Os    
            c = rem(Os_idx(Os),height);
            if c ==0
                c =height;
            end
            b = rem((Os_idx(Os)-c)/height, width)+1;
            a  = (Os_idx(Os)-(b-1)*height-c)/(height*width)+1;
          
           idx(1,:)=[a,b,c]  ;
           idx(2,:)=[a-1,b,c];
           idx(3,:)=[a+1,b,c];
           idx(4,:)=[a,b-1,c];
           idx(5,:)=[a,b+1,c];
           idx(6,:)=[a,b,c-1];
           idx(7,:)=[a,b,c+1];
         
          % Modify the neighbor index base on perodic bandary condition    
           
         %----------------------------------------------
         % Store the Os site index in  the array Os_idx|
         %----------------------------------------------
           
         site = Os_idx(Os);
         neighor_Os_counter=0;
         counter           =0;
                 
         for i=2:7
                
                flag=0;
                %-----------------------------------
                %Find the index of nearest neigbhor|
                %-----------------------------------
                
                if 0.5<idx(i,1) && idx(i,1)<length+0.5&&0.5<idx(i,2)&&...
                    idx(i,2)<width+0.5&& 0.5<idx(i,3)&&idx(i,3)<height+0.5
                                     
                        neighor=(idx(i,1)-1)*width*height+(idx(i,2)-1)*height+idx(i,3);
                %-------------------------------------
                % determine nearest neigbor atom type|
                %-------------------------------------
                 
                         flag = Sitename(neighor);
                %------------------------------------------------------------
                % Counter how many Os are next to the site we are concerning|
                %------------------------------------------------------------

                        neighor_Os_counter = neighor_Os_counter+flag;
                        counter            = counter + 1;
                %--------------------------------------------------------------------
                %Modify the off diagnal element value of the rate matrix accordingly|
                %--------------------------------------------------------------------

                         if flag == 0
                              Rate_Mtx(site,neighor) = transitionRate(1,2);
                              Rate_Mtx(neighor,site) = transitionRate(2,1);
                         else
                              Rate_Mtx(site,neighor) = transitionRate(2,2);
                              Rate_Mtx(neighor,site) = transitionRate(2,2);
                        end
                end
         end
                 
          Rate_Mtx(site,site)=-counter*transitionRate(2,1)-DecayRate(2)-neighor_Os_counter*transitionRate(2,2);
          Rate_Mtx(n,site) = DecayRate(2); 
 end
end
  
         
