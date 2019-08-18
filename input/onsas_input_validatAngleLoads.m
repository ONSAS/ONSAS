
inputONSASversion = '0.1.8';

problemName = 'validatAngleLoads' ;

E   = 1 ; 
nu  = 0 ;
rho = 0 ;

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu rho] ;

A = 1;
P = 1 ;

secGeomProps = [ A 1 1 1 ] ;

l1 = 1 ;

Nodes = [   0  0   0 ; ...
           l1  0   0 ; ... 
         2*l1  0   0 ; ... 
         3*l1  0   0 ] ;

kz2 = 0;

nodalSprings = [ 1  inf  inf  inf  inf  inf  inf ; ...
                 2    0  inf  inf  inf  kz2  inf ; ...
                 3    0  inf  inf  inf  kz2  inf ; ...
                 4    0  inf  inf  inf  inf  inf ] ;

Conec = [ 1 2 0 0 1 1 1 ; ...
          2 3 0 0 1 1 1 ; ...
          3 4 0 0 1 1 1 ] ;

nodalVariableLoads   = [ 2  0 0 0 0 -P 0 ] ;
nodalConstantLoads   = [ 4  .1 0 0 0 0 0 ] ;


bendStiff = [ 0 1 1 0 ] ;

% analysis parameters
nonLinearAnalysisBoolean = 1 ; 
dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ;

controlDofInfo = [ 2 5 -1 ] ;


printflag = 2 ;
tablesBoolean = 1;

plotParamsVector = [2];

plotsViewAxis = [ 0 -1 0] ;

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-10 ;
stopTolForces    = 1.0e-6 ;
targetLoadFactr  = 1e-2   ;
nLoadSteps       = 4    ;

numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 

% analytical solution 
analyticSolFlag = 0 ;
