
l = 3200 ;

ms = 0.05*l ;

Point(1) = {    0, 0,  0,     ms};
Point(2) = {    l, 0,  0,     ms};
Point(3) = {    l, l,  0, 0.5*ms};
Point(4) = {0.8*l, l,  0, 0.5*ms};
Point(5) = {    0, l,  0,     ms};

Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 4};
Line(14) = {4, 5};
Line(15) = {5, 1};

Line Loop(1) = {11, 12, 13, 14, 15};

Plane Surface(1) = {1};

Physical Point   ("00_01_01_00") = {1}; // fixed x and y
Physical Line    ("00_02_02_00") = {12}; // supported in x
Physical Point   ("00_01_03_00") = {3}; // loaded in y
Physical Surface ("01_03_00_00") = {1};
