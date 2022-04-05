function windVel = windVelNonLinearDynamic3D(x,t)
  % time when the fluid velocity is null
  zeroWindTime = 30 ;
  % period of sinusoidal wind
  T = zeroWindTime / 1 ; 
  %The angular frequency is
  Omega1 = 2*pi / T  ;
  %And the amplitude
  amplitude = 10     ;
  % High frecuency factor
  highFrec = 2.0716 ;
  omega2 = 2 * pi * highFrec ;
  % amplitude factor
  amplitudeFactorY = 0.2 ;
  amplitudeFactorZ = 0.5 ;
  % the axis y fluid velocity
  windy = ( amplitude * sin( Omega1 * t) + amplitude * amplitudeFactorY * sin( highFrec * t) ) .* (t < zeroWindTime) ;
  % the axis z fluid velocity
  windz = amplitude * amplitudeFactorZ  * sin( highFrec * t)  .* (t < zeroWindTime)  ;
    %yhe period of the sinusoidal component
  windVel = [0 windy windz]'; 
end