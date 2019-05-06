function f = myInfLinLoadFunc(t);

Fz = -1 ;
l = 2 ;

f = [[0 0 0 0  Fz * (l-t)^2  * [ 3*t + (l-t) ] / l^3 * (t<=l) 0]' ; ...
   [0 0 0 0  Fz * t^2 * [ t + 3*(l-t) ] / l^3 * (t<l) + Fz * (2*l-t)^2 * [ 3*(t-l) + (2*l-t) ] / (l)^3 * (t>=l) 0]' ; ...
   [0 0 0 0  Fz * (t-l)^2      * [ t-l + 3*(2*l-t) ] / (l)^3 * (t>=l) 0]'  ] ; 
