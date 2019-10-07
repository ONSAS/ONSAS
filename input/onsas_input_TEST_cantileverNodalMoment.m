% ---------------------------------------------as 
% ONSAS input file: cantilever with nodal moment example
% ---------------------------------------------

E = 200e9 ;  nu = 0.3 ;  l = 10   ; 

problemName = 'cantileverNodalMoment' ;
inputONSASversion = '0.1.9';

hyperElasParams = cell(1,1) ; hyperElasParams{1} = [1 E nu] ;

% --- cross section ---
b = .1 ;  h = .2 ;  Nelem = 20 ;

A  = b*h ;
Iy = b*h^3/12 ;
Iz = h*b^3/12 ;
It = 1 ;

sectPar = [ b h ]; 

secGeomProps = [ A Iy Iz It ] ;

nodalSprings = [ 1  inf  inf  inf  inf  inf  inf ] ;

Nodes = [ (0:(Nelem))'*l/Nelem zeros(Nelem+1,2) ] ;

Conec = [ (1:(Nelem))' (2:(Nelem+1))'  zeros(Nelem,2) (ones(Nelem,1)*[ 1 1 2]) ] ;

nodalVariableLoads   = [ Nelem+1  0 0 0 -1 0 0 ] ;

controlDofInfo = [ Nelem+1  4  -1 ] ;

% --- analysis parameters ---
nonLinearAnalysisBoolean = 1 ;
dynamicAnalysisBoolean   = 0 ; 

stopTolIts     = 30     ;
stopTolDeltau  = 1.0e-4 ;
stopTolForces  = 1.0e-4 ;
targetLoadFactr = 2*pi*E*Iy / l ;
nLoadSteps      = 20 ;

plotParamsVector = [ 2  2 ] ;    plotsViewAxis = [ 0 -1 0 ] ;
%~ plotParamsVector = [ 3 ] ;
printflag = 0 ;

numericalMethodParams = [ 1 ...
 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps ] ; 

analyticSolFlag = 1 ; analyticCheckTolerance = 1e-4 ;
analyticFunc = @(w) w*l / ( E * Iy )  ; 
