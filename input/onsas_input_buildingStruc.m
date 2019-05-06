
inputONSASversion = '0.1.8';
problemName = 'buildingStruc' ;

ndofpnode = 6 ;

% Material properties
E = 8500*(25+8)^(1/3)*100 ;
nu  = 0.3 ;

hyperElasParams = cell(1,1) ;  
hyperElasParams{1} = [ 1 E nu ] ;

% Geometric properties

% Flexible beam
a1 = 0.5 ;
b1 = 0.5 ;

%  rigid beams
Ab = a1*b1 ;
Izb = a1*b1^3/12 ;
Iyb = b1*a1^3/12 ;
beta = 0.141 ;
Itb = beta * b1 * a1^3 ;

secGeomProps = [ Ab   Iyb   Izb   Itb   ] ;

h = 3; nlevs = 6 ;

Nodes = [] ;

for i=1:nlevs
  Nodes = [ Nodes       ; ...
            0 0 (i-1)*h ; ... 
            4 0 (i-1)*h ; ...
            4 4 (i-1)*h ; ...
            0 4 (i-1)*h ] ;
end

nnodes = size(Nodes,1) ;

% Conectivity matrix 
Conec = [ ] ; 

for i=1:(nlevs-1)
  Conec = [ Conec; 
            (i-1)*4+1  (i  )*4+1 0 0 1 1 2 ; ...
            (i-1)*4+2  (i  )*4+2 0 0 1 1 2  ; ...
            (i-1)*4+3  (i  )*4+3 0 0 1 1 2  ; ...
            (i-1)*4+4  (i  )*4+4 0 0 1 1 2  ; ...
            (i  )*4+1  (i  )*4+2 0 0 1 1 2  ; ...
            (i  )*4+2  (i  )*4+3 0 0 1 1 2  ; ...
            (i  )*4+3  (i  )*4+4 0 0 1 1 2  ; ...
            (i  )*4+4  (i  )*4+1 0 0 1 1 2  ] ;
end
            
nelems = size(Conec,1) ; 

% Boundary conditions
nodalSprings = [ 1 inf inf inf inf inf inf ; ...
                 2 inf inf inf inf inf inf ; ...
                 3 inf inf inf inf inf inf ; ...
                 4 inf inf inf inf inf inf ] ;

%~ nodalConstantLoads   = [ 17  0 0 1 0 0 0 ] ;
nodalVariableLoads   = [ (nlevs-1)*4+1  0 0 1e-3 0 -1 0 ; ...
                         (nlevs-1)*4+2  0 0 1e-3 0 -1 0 ; ...
                         (nlevs-1)*4+3  0 0 0 0 -1 0 ; ...
                         (nlevs-1)*4+4  0 0 0 0 -1 0 ] ;


sectPar = [ 12 .5 .5 ];

% Analysis parameters
nonLinearAnalysisBoolean = 1 ; 
printflag = 2 ;
linearDeformedScaleFactor = 1.0 ;
controlDofInfo = [ nlevs*4 3 1 ] ;

% Plot options
%~ plotParamsVector = [ 2 4] ;
plotParamsVector = [ 3] ;
printflag = 0 ;

% Analytic solution 
analyticSolFlag = 0 ;

stopTolIts     = 30     ;
stopTolDeltau  = 1.0e-5 ;  stopTolForces  = 1.0e-8 ;
%~ targetLoadFactr = 1.0e3 ;  nLoadSteps     = 2    ;
%~ numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps ] ; 
targetLoadFactr = 1.0e5 ;  nLoadSteps     = 150    ;
numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps .05 ] ; 
