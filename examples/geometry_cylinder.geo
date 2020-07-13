ms = 0.01 ;
//ms = 0.009 ;
msi = 0.5*ms ;

largo_canio = 0.01 ;
radio_int   = 0.20 ;
radio_ext   = 0.24 ;

Point(1)  = { 0.0, 0, 0,  ms };
Point(2)  = { radio_int, 0, 0,  msi };
Point(3)  = { radio_ext, 0, 0,  ms  };
Point(4)  = { 0, radio_int, 0,  msi };
Point(5)  = { 0, radio_ext, 0,  ms  };

Point(6)  = {   0, 0, largo_canio , ms };
Point(7)  = { radio_int, 0, largo_canio , msi };
Point(8)  = { radio_ext, 0, largo_canio , ms };
Point(9)  = { 0, radio_int, largo_canio , msi };
Point(10) = { 0, radio_ext, largo_canio , ms };

Circle(1) 	= {2, 1, 4};
Line(2) 	= {4, 5} ;
Circle(3) 	= {5, 1, 3} ;
Line(4) 	= {3, 2} ;

Circle(1+4) = {2+5, 1+5, 4+5};
Line(2+4) 	= {4+5, 5+5} ;
Circle(3+4) = {5+5, 1+5, 3+5} ;
Line(4+4) 	= {3+5, 2+5} ;

Line(9) = {2,2+5} ;
Line(9+1) = {2+1,2+5+1};
Line(9+1+1) = {2+1+1,2+5+1+1};
Line(9+1+1+1) = {2+1+1+1,2+5+1+1+1};

Line Loop(1) = {1,2,3,4} ;
Line Loop(1+1) = {-5,-8,-7,-6} ;
Line Loop(1+1+1) = {-2,11,6,-12} ;
Line Loop(4) = {-4,10,8,-9} ;

Line Loop(5) = {-3, 12, 7, -10} ;
Line Loop(6) = {9, 5,-11,-1 } ;

Plane Surface(1) = {1};
Plane Surface(2) = {2}; 
Plane Surface(3) = {3}; 
Plane Surface(4) = {4};

Ruled Surface(5) = {5};
Ruled Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Physical Surface("05x00x00x01x00") = {3} ;
Physical Surface("05x00x00x02x00") = {4} ;
Physical Surface("05x00x00x03x00") = {1,2} ;

Physical Surface("05x00x00x00x00") = {5} ;
Physical Surface("05x00x00x00x01") = {6} ;

Physical Volume("03x01x01x00x00") = {1} ;

