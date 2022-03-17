function windVel = windVelNonLinearDynamic(x,t)
  finalTime = 5 ;
  constWindTime = finalTime / 10 ;
  vwindMax = 80 ;
  windy = t * vwindMax / constWindTime * (t <= constWindTime) + vwindMax * (t > constWindTime) ;
  windz =  0  ;

  windVel = [0 windy windz]'; 
end