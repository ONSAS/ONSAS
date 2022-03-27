function windVel = windVelNonLinearDynamic2D(x,t)
  finalTime = 200 ;
  constWindTime = finalTime / 30 ;
  vwindMax = 5 ;
  windy = t * vwindMax / constWindTime * (t <= constWindTime) + vwindMax * (t > constWindTime) ;
  windz =  0  ;

  windVel = [0 windy windz]'; 
end