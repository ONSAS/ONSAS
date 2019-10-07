% -------------------
% Example 3.1 from journal article: Li and Khandelwal, 'Topology
%  optimization of geometrically nonlinear trusses with spurious eigenmodes
%  control', Engineering Structures, 131 (2017) 324-344, 2017.
% -------------------

inputONSASversion = '0.1.9';
problemName = 'VonMisesTrussNR' ;

Es = 210e9 ;
nu = 0 ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu ] ;

% each row shows the properties of each section: A, Iy Iz and J
A = 2.5e-4 ;  secGeomProps = [ A 2 2 4 ] ;

% in global system of coordinates
nodalSprings = [ 1  inf  0  inf  0  inf 0 ; ...
                 2    0  0  inf  0    0 0 ; ...
                 3  inf  0  inf  0  inf 0   ...
               ];

auxx = cos(65*pi/180) * 2 ; auxy = sin(65*pi/180) * 2 ;

Nodes = [      0  0     0  ; ...
            auxx  0  auxy  ; ...
          2*auxx  0     0  ] ;

Conec = [ 1 2 0 0 1 1 1 ;
          2 3 0 0 1 1 1 ] ;

nodalVariableLoads   = [ 2  0  0  0  0 -1  0 ];

% [ node nodaldof scalefactor(positive or negative) ]
controlDofInfo = [ 2 5 -1 ] ;

% analysis parameters
nonLinearAnalysisBoolean = 1 ;  dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ; 

stopTolIts     = 30     ;
stopTolDeltau  = 1.0e-5 ;  stopTolForces  = 1.0e-8 ;
targetLoadFactr = 3.0e7 ;  nLoadSteps     = 50    ;

plotParamsVector = [ 2  4 ]; plotsViewAxis = [ 0 -1 0] ;
printflag = 1 ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps ] ; 

%~ stopTolIts       = 30     ;
%~ stopTolDeltau    = 1.0e-10 ;
%~ stopTolForces    = 1.0e-6 ;
%~ targetLoadFactr  = 5e7    ;
%~ nLoadSteps       = 50    ;
%~ incremArcLen     = .1    ;

%~ numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            %~ targetLoadFactr nLoadSteps incremArcLen ] ; 

% analytical solution using engineering strain
analyticSolFlag = 2 ; analyticCheckTolerance = 1e-4 ;
l0 = sqrt(auxx^2 + auxy^2) ;
analyticFunc = @(w) -2 * Es*A* ( (  (auxy+(-w)).^2 + auxx^2 - l0^2 ) ...
  ./ (l0 * ( l0 + sqrt((auxy+(-w)).^2 + auxx^2) )) ) ...
  .* ( auxy + (-w) ) ./ ( sqrt((auxy+(-w)).^2 + auxx^2) ) ; 
