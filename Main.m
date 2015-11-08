% This is the main program for 3D MOF


mytime=clock; mytimefmt=datestr(mytime);
fprintf(1,'Program Start at %s \n', mytimefmt);

% Read input files to get the value for basic rate constants and 
% MOF stcutrue information.
    
    kRuhop          = 1./50;
    kRuDecay        = 1./250.;
    kOshop          = 1./50;
    kOsDecay        = 1./243.;
    kOsTrap         = 1./50;
    kOsUntrap       = kOsTrap*exp(-12);
    Ospctg          = 0.0116;
    tsteps          = 2000;
    NumMOF          = 2^22;

    OsExciteOpt     = 0.0;

    ckeckstart      = 20;
    ckeckstep       = 12;
    outstep         = 6;
    checkpoint      = ckeckstart + ckeckstep ;
    tol_00          = 1.0d-03;
    tol_01          = tol_00 / 1.0;
    tol_02          = tol_01 / 1.0;

    tol             = tol_00;
    tol_flag        = 1;
  
% Stucture information

    length = 20;
    width  = 20;
    height = 20;
    n      = length*width*height+2;
    
% Initialization    
    
    OsPopu_av         = zeros(1,tsteps); 
    OsPopu_av_tmp     = zeros(1,tsteps);
    OsPopu             = zeros(1,tsteps);
    RuOssum         = zeros(1,tsteps);
    RuPopu_av       = zeros(1,tsteps);
    Num_Os          = floor( double(n-2)*Ospctg);
   
    
    % Store the rate constant in matrix form 
    
    [transitionRate,DecayRate] = ...
    initialization(kRuhop,kRuDecay,kOshop,kOsDecay,kOsTrap,kOsUntrap);
     


     RuPopu_Plus_RuEmi_av = zeros(1,tsteps);
     OsEmi_av             = zeros(1,tsteps);
     OsPopu_unNormal_av   = zeros(1,tsteps);
     
     
%Main Loop to gain tha average
for iternum=1:NumMOF
% Generate the rate matrix
 
    [Rate_Mtx,Osidx] = ...
    RateMat_3D(transitionRate,DecayRate,...
    length, width, height,Num_Os);
% % Calculate the Exp of Rate_Mtx 
   
    EXP_Rate_Mtx=expm(Rate_Mtx);
% Monitor the time propogate process.
    


    RuPopu_Plus_RuEmi = zeros(1,tsteps);
    OsEmi             = zeros(1,tsteps);
    OsPopu_unNormal   = zeros(1,tsteps);
    [OsPopu,RuOssum,RuPopu_Plus_RuEmi,OsEmi,OsPopu_unNormal]=  Time_Propgte(n, tsteps, Num_Os, EXP_Rate_Mtx, Osidx, OsExciteOpt);
% Add  to accumulate to gain the avarage value 

    RuPopu_Plus_RuEmi_av=RuPopu_Plus_RuEmi_av+RuPopu_Plus_RuEmi;
    OsEmi_av            = OsEmi_av + OsEmi;
    OsPopu_unNormal_av  = OsPopu_unNormal_av+OsPopu_unNormal;
     
    
    
    OsPopu_En_total(iternum,1:tsteps)     =  OsPopu;  
    RuOssum_all(iternum,1:tsteps) =  RuOssum;
    RuPopu_av=RuPopu_av +RuOssum;
    OsPopu_av = OsPopu_av + OsPopu;
% Check for converge
   
   if rem(iternum,outstep)==0
      fprintf(1,' Current iteration :  %d \n', iternum); 
   end
  
   if iternum==ckeckstart
      OsPopu_av_tmp=OsPopu_av;
   end
   flag = 0;
 
    % (2) Check if converged 
      if iternum == checkpoint 
         flag        = Cverg_check(iternum,ckeckstep,OsPopu_av,OsPopu_av_tmp,tol);
         OsPopu_av_tmp = OsPopu_av;
         fid=fopen('OsPopu_av.dat','W');
         fprintf(fid,'Iteration Number is %d \n', iternum );
         fprintf(fid, '%f\n', OsPopu_av/max(OsPopu_av));
         fclose(fid);
          
         fid=fopen('RuPopu_av.dat','W');
         fprintf(fid,'Iteration Number is %d \n', iternum );
         fprintf(fid, '%f\n', RuPopu_av/max(RuPopu_av));
         fclose(fid);
         
         
         RutoOs=zeros(1,tsteps-1);
         OsEmi_temp =zeros(1,tsteps-1);
         R=zeros(1,tsteps-1);
         
         for i=1:tsteps-1
              RutoOs(i)= RuPopu_Plus_RuEmi_av(i)-RuPopu_Plus_RuEmi_av(i+1);
              OsEmi_temp(i) = OsEmi_av(i)-OsEmi_av(i+1);
              R(i)= RutoOs(i)/OsEmi_temp(i);
         end
            
         fid=fopen('RUtoOs_VS_OsEmi.dat','W');
         fprintf(fid,'Iteration Number is %d \n', iternum );
         fprintf(fid, '%f \n', R);
         fclose(fid);
         
    
         
         fid=fopen('OsPopu_En_total.dat','W');
         fprintf(fid,'Iteration Number is %d \n', iternum );
         for outputloop=1:iternum
             fprintf(fid, '%f \n', OsPopu_En_total(outputloop,:));
         end
         fclose(fid);
           
         fid=fopen('RuOssum_all.dat','W');
         fprintf(fid,'Iteration Number is %d \n', iternum );
         for outputloop=1:iternum
             fprintf(fid, '%f \n', RuOssum_all(outputloop,:));
         end
         fclose(fid);

         checkpoint = checkpoint + ckeckstep;
      end
         
     % (3) If converged stop the program,or increase the number of loop by delta_n.
         if flag == 1

            if tol_flag == 1

                fprintf(1,'Calculation with tolerence 1 converged \n');

                OsPopu_av_tmp_0 = OsPopu_av / max(OsPopu_av);

                fid=fopen('OsPopu_av_00.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                fprintf(fid, '%f \n', OsPopu_av_tmp_0);
                fclose(fid);
                
                fid=fopen('RuPopu_av_00.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                fprintf(fid, '%f\n', RuPopu_av/max(RuPopu_av));
                fclose(fid);
         
          
                RutoOs=zeros(1,tsteps-1);
                OsEmi_temp =zeros(1,tsteps-1);
                R=zeros(1,tsteps-1);
          
               for i=1:tsteps-1
                    RutoOs(i)= RuPopu_Plus_RuEmi_av(i)-RuPopu_Plus_RuEmi_av(i+1);
                    OsEmi_temp(i) = OsEmi_av(i)-OsEmi_av(i+1);
                    R(i)= RutoOs(i)/OsEmi_temp(i);
               end
            
               fid=fopen('RutoOs_VS_OsEmi_00.dat','W');
               fprintf(fid,'Iteration Number is %d \n', iternum );
               fprintf(fid, '%f \n', R);
               fclose(fid);
             
                                   
              trap_efficency = (OsEmi_av+ OsPopu_unNormal_av)/(iternum*(n-2-Num_Os));
              fprintf(1,'EnT efficency in 3D  MOF = %8.3f \n',trap_efficency(tsteps));
              trap_efficency=0;
         
                fid=fopen('BackupOsPopu_En_total_00.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                for outputloop=1:iternum
                    fprintf(fid, '%f \n', OsPopu_En_total(outputloop,:));
                end
                fclose(fid);
                
                fid=fopen('BackupRuOssum_all_00.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                for outputloop=1:iternum
                    fprintf(fid, '%f \n', RuOssum_all(outputloop,:));
                end
                fclose(fid);
                
                tol = tol_01;
                tol_flag = tol_flag + 1;

             elseif tol_flag == 2
              
                fprintf(1,'Calculation with tolerence 2 converged \n'); 
                OsPopu_av_tmp_0 = OsPopu_av / max(OsPopu_av);
            
                fid=fopen('OsPopu_av_01.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                fprintf(fid, '%f \n', OsPopu_av_tmp_0);
                fclose(fid);

                
                fid=fopen('RuPopu_av_01.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                fprintf(fid, '%f\n', RuPopu_av/max(RuPopu_av));
                fclose(fid);
         
          
                RutoOs=zeros(1,tsteps-1);
                OsEmi_temp =zeros(1,tsteps-1);
                R=zeros(1,tsteps-1);
          
               for i=1:tsteps-1
                    RutoOs(i)= RuPopu_Plus_RuEmi_av(i)-RuPopu_Plus_RuEmi_av(i+1);
                    OsEmi_temp(i) = OsEmi_av(i)-OsEmi_av(i+1);
                    R(i)= RutoOs(i)/OsEmi_temp(i);
               end
            
               fid=fopen('RutoOs_VS_OsEmi_01.dat','W');
               fprintf(fid,'Iteration Number is %d \n', iternum );
               fprintf(fid, '%f \n', R);
               fclose(fid);
               
               trap_efficency = (OsEmi_av+ OsPopu_unNormal_av)/(iternum*(n-2-Num_Os));
               fprintf(1,'EnT efficency in 3D  MOF = %8.3f \n',trap_efficency(tsteps));
               trap_efficency=0;
               
               
               fid=fopen('BackupOsPopu_En_total_01.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                for outputloop=1:iternum
                    fprintf(fid, '%f \n', OsPopu_En_total(outputloop,:));
                end
                fclose(fid);
 
 
                fid=fopen('BackupRuOssum_all_01.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                for outputloop=1:iternum
                    fprintf(fid, '%f \n', RuOssum_all(outputloop,:));
                end
                fclose(fid);

                tol = tol_02;
                tol_flag = tol_flag + 1;

             elseif tol_flag == 3

                fprintf(1,'Calculation with tolerence 3 converged \n');
                OsPopu_av_tmp_0 = OsPopu_av / max(OsPopu_av);

                fid=fopen('OsPopu_av_02.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                fprintf(fid, '%f \n', OsPopu_av_tmp_0);
                fclose(fid);

                fid=fopen('RuPopu_av_02.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                fprintf(fid, '%f\n', RuPopu_av/max(RuPopu_av));
                fclose(fid);
         
          
                RutoOs=zeros(1,tsteps-1);
                OsEmi_temp =zeros(1,tsteps-1);
                R=zeros(1,tsteps-1);
          
               for i=1:tsteps-1
                    RutoOs(i)= RuPopu_Plus_RuEmi_av(i)-RuPopu_Plus_RuEmi_av(i+1);
                    OsEmi_temp(i) = OsEmi_av(i)-OsEmi_av(i+1);
                    R(i)= RutoOs(i)/OsEmi_temp(i);
               end
            
               fid=fopen('RutoOs_VS_OsEmi_02.dat','W');
               fprintf(fid,'Iteration Number is %d \n', iternum );
               fprintf(fid, '%f \n', R);
               fclose(fid);
              
               trap_efficency = (OsEmi_av+ OsPopu_unNormal_av)/(iternum*(n-2-Num_Os));
               fprintf(1,'EnT efficency in 3D  MOF = %8.3f \n',trap_efficency(tsteps));
               trap_efficency=0;        
                
                fid=fopen('BackupOsPopu_En_total_02.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                for outputloop=1:iternum
                    fprintf(fid, '%f \n', OsPopu_En_total(outputloop,:));
                end
                fclose(fid);
                
                fid=fopen('BackupRuOssum_all_02.dat','W');
                fprintf(fid,'Iteration Number is %d \n', iternum );
                for outputloop=1:iternum
                    fprintf(fid, '%f \n', RuOssum_all(outputloop,:));
                end
                fclose(fid);

                break;

             end

         end


         if iternum == NumMOF
            fprintf(1,'Max iteration number reached. Calculation does not converged\n');
            fid=fopen('OsPopu_av.dat','W');
            fprintf(fid,'Iteration Number is %d \n', iternum );
            fprintf(fid, '%f\n', OsPopu_av);
            fclose(fid);

            fid=fopen('OsPopu_En_total.dat','W');
            fprintf(fid,'Iteration Number is %d \n', iternum );
            for outputloop=1:iternum
                fprintf(fid, '%f \n', OsPopu_En_total(outputloop,:)); 
            end
            fclose(fid);

            fid=fopen('RuOssum_all.dat','W');
            fprintf(fid,'Iteration Number is %d \n', iternum );
            for outputloop=1:iternum
               fprintf(fid, '%f \n', RuOssum_all(outputloop,:));
            end
            fclose(fid);

         end 
         
end

            
      quit;
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
