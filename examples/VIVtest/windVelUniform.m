function windVel = windVelUniform(x,t)
  global vwindMax;
%   finalTime = 0.2 ;
%   constWindTime = finalTime / 10 ;
  % constant profile
  windx =  0.035 ; % Ur = 5.6
  % ramp profile
  %windx = t * vwindMax / constWindTime * (t <= constWindTime) + vwindMax * (t > constWindTime) ;
  windVel = [-windx 0 0]'; 
end