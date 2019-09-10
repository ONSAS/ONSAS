
inputONSASversion = '0.1.8';

problemName = 'validatAngleLoadsEx1' ;

E   = 200e9 ; 
nu  = 0 ;
rho = 0 ;

Iy  = (.1^4)/12 ;

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu rho] ;

A = 1   ;
P = 1e4 ;

secGeomProps = [ A 1 1 1 ] ;

%~ nelem = 2 ;
 %~ nelem = 4 ;
%~ nelem = 8 ;
%~ nelem = 16 ;
%~ nelem = 32 ;
nelem = 64 ;

L = 2 ;
l1 = L / nelem ;

nnodes = nelem +1;

if mod(nelem,2)~=0, stop, end

middleNode = 1+nelem/2 ;

Nodes = [ (0:(nnodes-1))'*l1 zeros(nnodes,2) ] ;

kz2 = 0;

nodalSprings = [ 1       inf  inf  inf  inf  inf  inf ; ...
                 nnodes  inf  inf  inf  inf  inf  inf ] ;

Conec = [ (1:nelem)' (2:nnodes)' zeros(nelem,2) ones(nelem,3) ] ;

nodalVariableLoads   = [ middleNode    0 0 0 0 -P 0 ] ;
%~ nodalConstantLoads   = [ nelem+1    1e-3 0 0 0  0 0 ] ;

EI = E*Iy ;

bendStiff = [ 0 ones(1,nnodes-2)*EI 0 ] ;

% analysis parameters
nonLinearAnalysisBoolean = 1 ; 
dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ;

controlDofInfo = [ middleNode  5  -1 ] ;

nonHomogeneousInitialCondU0 = [ 2 3 -.001 ] ;

printflag     = 0 ;
tablesBoolean = 1;

plotParamsVector = [ 2 2 ];

plotsViewAxis = [ 0 -1 0] ;

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-10 ;
stopTolForces    = 1.0e-6 ;
targetLoadFactr  = 1    ;
nLoadSteps       = 1    ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

% analytical solution 
analyticSolFlag = 0 ;
analytSol = [ -P*L^3 / (48*E*Iy) ] ;
