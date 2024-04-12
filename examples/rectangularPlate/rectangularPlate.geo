

// sizes for another case: rInt = 0.05 ; rExt = 0.15 ;
Lx = 2 ; Ly = 1 ;

ms = .08 ; //

Point(1) = {  0, 0,  0, ms};
Point(2) = { Lx, 0,  0, ms};
Point(3) = { Lx, Ly, 0, ms};
Point(4) = {  0, Ly, 0, ms};

Line(1) = {1, 2}; //
Line(2) = {2, 3}; //
Line(3) = {3, 4}; //
Line(4) = {4, 1}; //

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Line("00_01_01_00") = {1,3};
Physical Line("00_01_00_00") = {2,4};

Physical Surface("01_02_02_00") = {1};
