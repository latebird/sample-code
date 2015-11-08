function [flag]=Cverg_check(n,delta_n,Fi,Fi_pre,tol)

   Fi    = Fi/n;

   Fi_pre= Fi_pre/(n-delta_n);

   Fi    = Fi-Fi_pre;
   Fi    = abs(Fi);
   temp  = max(Fi);
   flag  = 0;

   if temp < tol
      flag=1;
   end 
   
end
   
