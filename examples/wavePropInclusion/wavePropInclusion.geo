lx = 1.0    ;
ly = 0.25*lx ;

ms = 0.03*lx ;

r  = 0.4*ly ;
x0 = 0.5*lx ;
y0 = 0.5*ly ;

Point(1) = {0,       0, 0, ms};
Point(2) = {lx,      0, 0, ms};
Point(3) = {lx,     ly, 0, ms};
Point(4) = {0,      ly, 0, ms};

Point(5) = {x0-r , y0, 0, 0.5*ms};
Point(6) = {x0   , y0, 0, 0.5*ms};
Point(7) = {x0+r , y0, 0, 0.5*ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Circle(5) = {5, 6, 7};
Circle(6) = {7, 6, 5};

Line Loop(1) = {1, 2, 3, 4,-5,-6};
Line Loop(2) = {5, 6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Point ("00_01_01_00") = {1};
Physical Point ("00_01_02_00") = {2};
Physical Line  ("00_01_03_00") = {2};

Physical Surface("01_03_00_00") = {1};
Physical Surface("02_03_00_00") = {2};
