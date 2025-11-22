
// Dimensions
Lx = 10 ;
Ly = 1.0 ;

ms = 0.5 ;
f = 1 ;

Point(1) = {0,0,0, ms} ;
Point(2) = {Lx,0,0, f*ms} ;
Point(3) = {Lx,Ly,0, f*ms} ;
Point(4) = {0,Ly,0, ms} ;

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Curve ("00_01_01") = {4};
Physical Curve ("00_01_02") = {2};
Physical Surface ("01_02_00") = {1};
