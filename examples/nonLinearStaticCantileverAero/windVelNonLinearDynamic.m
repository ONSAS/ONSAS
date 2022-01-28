function windVel = windVelNonLinearDynamic(x,t)
  finalTime = 10 ;
  constWindTime = finalTime/2 ;
  vwindMax = 25 ;
  windy = t * vwindMax / constWindTime * (t <= constWindTime) + vwindMax * (t > constWindTime) ;
  windz =  0  ;

  windVel = [0 windy windz]'; 
end