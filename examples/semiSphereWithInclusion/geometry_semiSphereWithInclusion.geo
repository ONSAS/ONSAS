Mesh.Algorithm = 6;
lc     =  .0125 ;
lc2    =  .025 ;

radio  = 0.03 ;         // radius inclusion
radio2 = 0.15 ;         // radius semi-sphere
radio3 = 0.02  ;         // radius pressure
alt1   = radio2*.5 ;
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
Point(13) = {0        , 0+alt2,-radio2,lc2};
Point(14) = { -radio2 , 0+alt2 , 0.0 , lc2};
Point(15) = {0        , 0+alt2,radio2,lc2};

alt3 = Sqrt( radio2*radio2 - radio3*radio3  ) ;

Point(16) = {  radio3 , alt3 , 0.0 , lc2};
Point(17) = {0        , alt3,-radio3,lc2};
Point(18) = { -radio3 , alt3 , 0.0 , lc2};
Point(19) = {0        , alt3,radio3,lc2};

Point(20) = {0        , alt3,0,lc2};

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

Circle(21) = {12,11,16};
Circle(22) = {13,11,17};
Circle(25) = {14,11,18};
Circle(28) = {15,11,19};

Circle(29) = {12,11,13};
Circle(210) = {13,11,14};
Circle(211) = {14,11,15};
Circle(212) = {15,11,12};

Circle(31) = {16,20,17};
Circle(32) = {17,20,18};
Circle(33) = {18,20,19};
Circle(34) = {19,20,16};

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

Line Loop(111) = {21,-34,-28,212};
Line Loop(112) = {22,-31,-21,29};
Line Loop(113) = {25,-32,-22,210};
Line Loop(114) = {28,-33,-25,211};

// pressure
Line Loop(115) = {31,32,33,34};

// bottom
Line Loop(116) = {-212,-211,-210,-29};

Ruled Surface(121) = {111};
Ruled Surface(122) = {112};
Ruled Surface(123) = {113};
Ruled Surface(124) = {114};
Ruled Surface(125) = {115};
Ruled Surface(126) = {116};

Surface Loop(29) = {28,26,16,14,20,24,22,18};
Volume(30) = {29};

//Surface Loop(31) = {-28,-26,-16,-14,-20,-24,-22,-18,-114,-118,-120,-122,-123};
Surface Loop(31) = {-28,-26,-16,-14,-20,-24,-22,-18,121,122,123,124,125,126};
Volume(32) = {31};
//+

Physical Surface("00_01_01_00") = {125};
Physical Surface("00_01_02_00") = {126};

// inclusion
Physical Volume("01_02_00_00") = {30};

// semi-sphere
Physical Volume("02_02_00_00") = {32};
