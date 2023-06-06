//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 500, 0, 2*Pi};
//+
Circle(2) = {0, 0, 42.16, 250, 0, 2*Pi};
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/3} {
  Duplicata { Point{2}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/3} {
  Duplicata { Point{3}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/3} {
  Duplicata { Point{4}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/3} {
  Duplicata { Point{5}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/3} {
  Duplicata { Point{6}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{1}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{8}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{9}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{10}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{11}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{12}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{13}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{14}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{15}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{16}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
  Duplicata { Point{17}; }
}
//+
Point(19) = {0, 0, 72.16, 1};
//+
Line(3) = {19, 2};
//+
Line(4) = {19, 3};
//+
Line(5) = {19, 4};
//+
Line(6) = {19, 5};
//+
Line(7) = {19, 6};
//+
Line(8) = {19, 7};
//+
Line(9) = {2, 3};
//+
Line(10) = {3, 4};
//+
Line(11) = {4, 5};
//+
Line(12) = {5, 6};
//+
Line(13) = {6, 7};
//+
Line(14) = {7, 2};
//+
Line(15) = {2, 8};
//+
Line(16) = {8, 3};
//+
Line(17) = {3, 10};
//+
Line(18) = {10, 4};
//+
Line(19) = {4, 12};
//+
Line(20) = {12, 5};
//+
Line(21) = {5, 14};
//+
Line(22) = {14, 6};
//+
Line(23) = {6, 16};
//+
Line(24) = {16, 7};
//+
Line(25) = {7, 18};
//+
Line(26) = {18, 2};
//+
Recursive Delete {
  Point{17}; Point{15}; Point{13}; Point{11}; Point{9}; Curve{1}; Curve{2}; 
}

Physical Point("00_01_01_00") = {8,10,12,14,16,18};

Physical Point("00_01_02_00") = {19};

Physical Line ("01_02_00_00") = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
