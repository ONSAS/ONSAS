
// Dimensions
Lx = 2.0 ;
Ly = 1.0 ;
Lz = 1.0 ;

ms = 0.5 ;

Point(1) = {0,0,0, ms} ;
Point(2) = {0,0,Lz, ms} ;
Point(3) = {0,Ly,Lz, ms} ;
Point(4) = {0,Ly,0, ms} ;

Point(5) = {Lx,0,0, ms} ;
Point(6) = {Lx,0,Lz, ms} ;
Point(7) = {Lx,Ly,Lz, ms} ;
Point(8) = {Lx,Ly,0, ms} ;

Line(1) = {4, 3};
Line(2) = {3, 7};
Line(3) = {7, 8};
Line(4) = {8, 4};

Line(5) = {1, 2};
Line(6) = {2, 6};
Line(7) = {6, 5};
Line(8) = {5, 1};

Line(9) = {1, 4};
Line(10) = {3, 2};
Line(11) = {5, 8};
Line(12) = {7, 6};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {9, -4, -11, 8};
Plane Surface(2) = {2};

Curve Loop(3) = {11, -3, 12, 7};
Plane Surface(3) = {3};

Curve Loop(4) = {10, 6, -12, -2};
Plane Surface(4) = {4};

Curve Loop(5) = {5, -10, -1, -9};
Plane Surface(5) = {5};

Curve Loop(6) = {-8, -7, -6, -5};
Plane Surface(6) = {6};

Surface Loop(1) = {5, 6, 2, 1, 4, 3};
Volume(1) = {1};

//+
Physical Surface ("00_01_02_00") = {5};
Physical Surface ("00_01_03_00") = {6};
Physical Surface ("00_01_04_00") = {2};
Physical Surface ("00_01_01_00") = {3};
//+
Physical Volume  ("01_02_00_00") = {1};
