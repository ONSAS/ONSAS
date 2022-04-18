function windVel = windVelNonLinearStatic(x,t)
  windy = 5*t     ;
  windz =  0  ;
  windVel = [0 windy windz]'; 
end