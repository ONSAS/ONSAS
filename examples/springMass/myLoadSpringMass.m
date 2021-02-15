function f = myLoadSpringMass( t)

%Force data
omegaBar = 2   ;
p0       = 0.1 ;

f=zeros(12,1);
f(7) = p0 *sin( omegaBar*t) ; 
