% This is the function simualte the excitation propogation with a defined 
% rate matrix. The return of this function is the a time-dependent Ru excited 
% state populaiton.

function [F,RuOssum,RuPopu_Plus_RuEmi,OsEmi,OsPopu_unNormal] = Time_Propgte(lengthmof,tsteps,num_Os,U,m, OsExciteOpt)

     
     RuPopu_Plus_RuEmi = zeros(1,tsteps);
     OsEmi             = zeros(1,tsteps);
     
     % define the intial excited state 
     C0 =ones(1,lengthmof);
     C0(lengthmof-1)= 0;
     C0(lengthmof)  = 0;
     
     F  = zeros(1,tsteps);
     P  = zeros(lengthmof,tsteps);
     OsPopu_unNormal= zeros(1,tsteps);
     RuOssum =zeros(1,tsteps);
      
        for i=1:num_Os
           C0(m(i))=1.0*OsExciteOpt;
        end 
        
     P(1:lengthmof,1)=C0;
     x=sum(P(:,1));
     
     % Doing the time propagation
     for k=2:tsteps
         P(:,k) =  U*P(:,k-1);
         x0=sum(P(:,k))/x;
         P(:,k) = P(:,k)/x0;        
     end 
     
     % Computing the excited sates population on both Ru and Os
     for k=1:num_Os
         F=P(m(k), 1:tsteps)+F;
     end 
     
     for k=1:lengthmof-2
         RuOssum =P(k, 1:tsteps)+RuOssum;
     end 
     Fmax =  max(F);
     RuOssum=RuOssum-F;
     RuPopu_Plus_RuEmi=P(lengthmof-1,:)+RuOssum;
     OsEmi=P(lengthmof,:);
     RuOssum=RuOssum/max(RuOssum);
     
     OsPopu_unNormal=F;
     F = F/Fmax;
     
     
end
