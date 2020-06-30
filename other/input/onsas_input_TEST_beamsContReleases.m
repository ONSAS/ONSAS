% ------------------------------------
% Test problem: continous beam with internal releases.
% ------------------------------------


inputONSASversion = '0.1.9';
problemName = 'beamsContReleases' ;

% Constitutive properties
E  = 2.1e5 ; nu = 0.3 ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% Geometrical properties
l = 2   ; b = 0.1 ; % width given by the local axis y
h  = b ; % fixed-ratio heigth given by local z
A  = b*h ; Iy = (b*h^3) / 12 ; Iz = (h*b^3) / 12 ; It = 1.0 ;
secGeomProps = [ A Iy Iz It ] ;

% Nodes
Nodes = [ -1.5*l 0 0  ; ...
          -0.5*l 0 0  ; ... 
               0 0 0  ; ... 
           0.5*l 0 0  ;
           1.5*l 0 0  ] ;

% Boundary conditions
nodalSprings = [ 1  inf  inf  inf  inf  inf  inf  ; ...
                 5  inf  inf  inf  inf  inf  inf  ] ;

% Conectivty matrix
Conec = [ 1 2 0 0 1 1 2 ; ...
          2 3 0 0 1 1 2 ; ...
          3 4 0 0 1 1 2 ; ...
          4 5 0 0 1 1 2 ] ;

% Releases
Releases = [ 1 0 1 0 1 ; ...
             4 1 0 1 0 ] ;

% Nodal loads
Fload =  1 ; Fy = Fload*1 ; Fz = Fload*1 ;
nodalConstantLoads = [ 3 0 0 Fy 0 Fz 0 ] ;

% Analysis parameters
nonLinearAnalysisBoolean = 0 ; printflag = 2 ; tablesBoolean = 1 ;
plotParamsVector = [2] ;
linearDeformedScaleFactor = 1;

% Analytical solution 
analyticSolFlag = 4 ; G  = E / ( 2*(1+nu) ) ; Uanalyt = zeros(5*6,1) ;

% Fuerza y
deltaymens =  (0.5*Fy)*l^3/( 3*E*Iz) ;
thetazmens =  (0.5*Fy)*l^2/( 2*E*Iz) ;

deltaysimp =  (    Fy)*l^3/(48*E*Iz) ;
thetazsimp =  (    Fy)*l^2/(16*E*Iz) ;

% Fuerza z
deltazmens =  (0.5*Fz)*l^3/( 3*E*Iy) ;
thetaymens =  (0.5*Fz)*l^2/( 2*E*Iy) ;

deltazsimp =  (    Fz)*l^3/(48*E*Iy) ;
thetaysimp =  (    Fz)*l^2/(16*E*Iy) ;

analytSol = zeros( 4,2*6) ;
analytSol(1,6+(3:6) ) = [ deltaymens -thetaymens deltazmens  thetazmens ] ;

analytSol(2,0+(3:6) ) = [ deltaymens -thetaysimp deltazmens  thetazsimp ] ;
analytSol(2,6+(3:6) ) = [ deltaymens+deltaysimp   0   deltazmens+deltazsimp 0 ] ;

analytSol(3,0+(3:6) ) = [ deltaymens+deltaysimp   0   deltazmens+deltazsimp 0 ] ;
analytSol(3,6+(3:6) ) = [ deltaymens +thetaysimp  deltazmens -thetaysimp ] ;

analytSol(4,0+(3:6) ) = [ deltaymens +thetaymens  deltazmens -thetaymens ] ;

analyticSolDofs = (1:(5*6))' ;
analyticCheckTolerance = 1e-4 ;
