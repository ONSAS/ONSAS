% -------------------
% Example 3.1 from journal article: Li and Khandelwal, 'Topology
%  optimization of geometrically nonlinear trusses with spurious eigenmodes
%  control', Engineering Structures, 131 (2017) 324-344, 2017.
% -------------------

inputONSASversion = '0.1.8';  problemName = 'ej2b7' ;

Es = 2.1e5 ;
nu = 0 ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu ] ;

% each row shows the properties of each section: A, Iy Iz and J
A = 10000 ;  secGeomProps = [ A 1 1 1 ] ;

kres = 2100  ;

% in global system of coordinates
nodalSprings = [ 1  inf  0  inf  0  inf 0 ; ...
                 2    0  0  inf  0  kres 0 ; ...
                 3    0  0  inf  0  kres 0 ; ...
                 4    0  0  inf  0  inf 0   ...
               ];
auxx = 1000

Nodes = [      0  0     0  ; ...
            auxx  0     0  ; ...
          2*auxx  0     auxx*.0001  ; ...
          3*auxx  0     0  ] ;

Conec = [ 1 2 0 0 1 1 1 ; ...
          2 3 0 0 1 1 1 ; ...
          3 4 0 0 1 1 1 ] ;

nodalVariableLoads   = [ 4  -1  0  0  0  0  0 ];

% [ node nodaldof scalefactor(positive or negative) ]
controlDofInfo = [ 3 5 -1 ] ;

% analysis parameters
nonLinearAnalysisBoolean = 1 ;  dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ; 

stopTolIts     = 30     ;
stopTolDeltau  = 1.0e-5 ;  stopTolForces  = 1.0e-8 ;
%~ targetLoadFactr = 8e4 ;  nLoadSteps     = 200    ;
targetLoadFactr = kres*auxx/3*1.1 ;  nLoadSteps     = 200    ;
incremArcLen     = .01    ;


%~ plotParamsVector = [ 2  4 ];
plotParamsVector = [ 3 ];
plotsViewAxis = [ 0 -1 0] ;
printflag = 0 ;

numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps  incremArcLen] ; 
%~ numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps  ] ; 
