
l = 3200 ;

ms = 0.05*l ;
//ms = 0.5*l ;

Point(1) = {0, 0,  0,     ms};
Point(2) = {l, 0,  0,     ms};
Point(3) = {l, l,  0, 0.5*ms};
Point(4) = {0, l,  0,     ms};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(1) = {5, 6, 7, 8};

Plane Surface(1) = {1};

Physical Point   ("00_01_01_00")  = {1}; // fixed x and y
Physical Line    ("00_02_02_00")  = {6}; // supported in x
Physical Point   ("00_01_03_00")  = {3}; // loaded in y
Physical Surface ("01_03_00_00") = {1};
