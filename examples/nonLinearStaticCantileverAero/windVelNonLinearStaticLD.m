function windVel = windVelNonLinearStaticLD(x,t)
  windy = 100*t;
  windz =  0  ;
  windVel = [0 windy windz]'; 
end