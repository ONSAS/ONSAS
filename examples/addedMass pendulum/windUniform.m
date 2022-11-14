function windVel = windUniform(x,t)
global finalTime; global vwindMax;
  constWindTime  = 0;%finalTime/10 ;
  %windx =  vwindMax / constWindTime * (t <= constWindTime) + vwindMax * (t > constWindTime) ;
  f = 0.2;
  Amp = 0.01;
  windx =  Amp*sin(2*pi*f*t);
  %windx =  Amp *t* (t <= constWindTime) - Amp *t*(t > constWindTime) ;
  windVel = [-windx 0 0]';
end