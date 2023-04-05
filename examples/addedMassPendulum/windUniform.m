function windVel = windUniform(x,t)
global finalTime; global vwindMax;
  f = 0.2;
  Amp = 0.01;
  windx =  Amp*sin(2*pi*f*t);
  windVel = [-windx 0 0]';
end
