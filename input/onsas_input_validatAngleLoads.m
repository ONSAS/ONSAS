
inputONSASversion = '0.1.8';

problemName = 'validatAngleLoads' ;

E   = 800 ; 
nu  = 0 ;
rho = 0 ;

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu rho] ;

A = 1;
P = 1 ;

secGeomProps = [ A 1 1 1 ] ;

nelem = 10 ;

%Pl3/48EI

l1 = 2/nelem ;

nnodes = nelem +1;

Nodes = [ (0:(nnodes-1))'*l1 zeros(nnodes,2) ] ;

kz2 = 0;

nodalSprings = [ 1       inf  inf  inf  inf  inf  inf ; ...
                 %~ nnodes    0  inf  inf  inf  inf  inf ] ;
                 nnodes   inf  inf  inf  inf  inf  inf ] ;

Conec = [ (1:nelem)' (2:nnodes)' zeros(nelem,2) ones(nelem,3) ] ;

nodalVariableLoads   = [ 1+nelem/2 0 0 0 0 -1 0 ] ;

I = 
EI = 800/48 ;

bendStiff = [ 0 ones(1,nnodes-2)*EI 0 ] ;


% analysis parameters
nonLinearAnalysisBoolean = 1 ; 
dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ;

controlDofInfo = [ 1+nelem/2 5  -1 ] ;

nonHomogeneousInitialCondU0 = [ 2 5 -.001 ] ;

printflag     = 2 ;
tablesBoolean = 1;

plotParamsVector = [2 ];

plotsViewAxis = [ 0 -1 0] ;

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-10 ;
stopTolForces    = 1.0e-6 ;
targetLoadFactr  = 50    ;
nLoadSteps       = 20    ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

% analytical solution 
analyticSolFlag = 0 ;
