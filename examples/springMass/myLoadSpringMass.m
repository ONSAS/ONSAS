function f = myLoadSpringMass( t)

%Force data
k        = 39.47 ; % spring constant
m        = 1     ; % mass of the system
p0       = 40    ; % amplitude of applied load

%md The free vibration motion parameters are:
omegaN       = sqrt( k / m )           ;
%md The frequency of the sinusoidal external force is set:
omegaBar     = 4*omegaN ;

f=zeros(12,1);
f(7) = p0 *sin( omegaBar*t) ; 
