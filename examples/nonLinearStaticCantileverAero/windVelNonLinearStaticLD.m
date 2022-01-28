function windVel = windVelNonLinearStaticLD(x,t)
  windy = 25*t;
  windz =  0  ;
  windVel = [0 windy windz]'; 
end