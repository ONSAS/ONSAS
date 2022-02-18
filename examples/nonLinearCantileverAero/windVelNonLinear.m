function windVel = windVelNonLinear(x,t)
  windy = 40*t     ;
  windz =  0  ;
  windVel = [0 windy windz]'; 
end