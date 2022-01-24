Mesh.Algorithm = 6;
lc     =  .50 ;
lc2    =  .9 ;

radio  =  .70 ;         // radius inclusion
radio2 = 4.00 ;         // radius semi-sphere
alt1   = 3.00*radio ;
alt2   = 0 ;

angCarg = 0.1*Pi ;

// --- points inclusion ---
Point(1) = {0.0    ,        +alt1 ,    0.0  , lc};
Point(2) = {radio  ,        +alt1 ,    0.0  , lc};
Point(3) = {0      ,  radio +alt1 ,    0.0  , lc};
Point(4) = {-radio ,        +alt1 ,    0.0  , lc};
Point(5) = {0      , -radio +alt1 ,    0.0  , lc};
Point(6) = {0      , 0      +alt1 ,  -radio , lc};
Point(7) = {0      , 0      +alt1 ,   radio , lc};
// -----------------------------------------------


// --- points semi-sphere ---
Point(11) = {0.0      , 0.0+alt2 , 0.0,lc2};
Point(12) = {  radio2 , 0.0+alt2 , 0.0 , lc2};
Point(13) = {0        , radio2+alt2 , 0.0 , lc2};
Point(14) = { -radio2 , 0+alt2 , 0.0 , lc2};
Point(16) = {0        , 0+alt2,-radio2,lc2};
Point(17) = {0        , 0+alt2,radio2,lc2};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {3,1,6};
Circle(6) = {6,1,5};
Circle(7) = {5,1,7};
Circle(8) = {7,1,3};
Circle(9) = {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};

Circle(21) = {12,11,13};
Circle(22) = {13,11,14};
Circle(25) = {13,11,16};
Circle(28) = {17,11,13};
Circle(29) = {12,11,17};
Circle(210) = {17,11,14};
Circle(211) = {14,11,16};
Circle(212) = {16,11,12};

Line Loop(13) = {2,8,-10};
Ruled Surface(14) = {13};
Line Loop(15) = {10,3,7};
Ruled Surface(16) = {15};
Line Loop(17) = {-8,-9,1};
Ruled Surface(18) = {17};
Line Loop(19) = {-11,-2,5};
Ruled Surface(20) = {19};
Line Loop(21) = {-5,-12,-1};
Ruled Surface(22) = {21};
Line Loop(23) = {-3,11,6};
Ruled Surface(24) = {23};
Line Loop(25) = {-7,4,9};
Ruled Surface(26) = {25};
Line Loop(27) = {-4,12,-6};
Ruled Surface(28) = {27};

Line Loop(113) = {22,28,-210};
Ruled Surface(114) = {113};

Line Loop(117) = {-28,-29,21};
Ruled Surface(118) = {117};

Line Loop(119) = {-211,-22,25};
Ruled Surface(120) = {119};

Line Loop(121) = {-25,-212,-21};
Ruled Surface(122) = {121};

Line Loop(200) = {29,210,211,212};
Plane Surface(123) = {200};

Surface Loop(29) = {28,26,16,14,20,24,22,18};
Volume(30) = {29};

Surface Loop(31) = {-28,-26,-16,-14,-20,-24,-22,-18,-114,-118,-120,-122,-123};
Volume(32) = {31};
//+

Physical Surface("00_01_01_00") = {123};
Physical Surface("00_01_02_00") = {120};

// inclusion
Physical Volume("01_02_00_00") = {30};

// semi-sphere
Physical Volume("02_02_00_00") = {32};
