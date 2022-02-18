function windVel = windVelDynamic(x,t)
  % finalTime = 150  ;
  % vwindMax  = 3 ;
  % thetaMax  = 40  ;
  % constantWindTime = 20 ;
  % thetaVel  = deg2rad(t / constantWindTime * thetaMax * (t <= constantWindTime) + thetaMax * (t > constantWindTime) ) ;
  % windy =  t * vwindMax * cos( thetaVel ) / constantWindTime * (t <= constantWindTime) + vwindMax * cos( deg2rad(thetaMax) ) * (t > constantWindTime) ;
  % windz = -t * vwindMax * sin( thetaVel ) / constantWindTime * (t <= constantWindTime) - vwindMax * sin( deg2rad(thetaMax) ) * (t > constantWindTime) ;
  % windVel = [ 0 windy windz ]' ; %Must be a column vector
  windVel = [10 0 0]' ;
end 