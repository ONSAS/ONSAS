function windVel = windVel (x, t)
  windy = 0;
  windz = 10*t ;
  
  windVel = [0 windy windz]'; 
end