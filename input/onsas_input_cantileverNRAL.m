% --------------------------------------------------------------------------------------------------
% Test problem: rectangular cross-section cantilever beam submitted to nodal external loads at the free end.
% --------------------------------------------------------------------------------------------------

E  = 210e9   ;
nu = 0.3     ;
l  = 2       ; % reference length
b1 = .1      ;
b2 = .05     ;

inputONSASversion = '0.1.8';

problemName = 'cantileverNRAL' ;

% constitutive properties

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% geometrical properties

secGeomProps = [ b1*b2 b1*b2^3/12 b2*b1^3/12 1 ] ;

Nodes = [   0  0  0 ; ...
            0  0  l ] ;

% in global system of coordinates
nodalSprings = [ 1  inf  inf  inf  inf  inf  inf ] ;

Conec = [ 1 2 0 0 1 1 2 ] ;

nodalVariableLoads   = [ 2  0 0 0 1e-5 -1 0 ] ;

% analysis parameters
nonLinearAnalysisBoolean = 1 ; 
dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ;

% [ node nodaldof scalefactor(positive or negative) ]
controlDofInfo = [ 2 1 1 ] ;

printflag = 0 ;

plotParamsVector = [3 50];

sectPar = [ 12 b1 b2 ];

% analytical solution 
analyticSolFlag = 0 ;

Pcr = E * pi^2 * b1*b2^3/12 / ( 2*l)^2 ;

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-8 ;
stopTolForces    = 1.0e-8 ;
targetLoadFactr  = Pcr*2    ;
nLoadSteps       = 10000    ;
incremArcLen     = 1e-4    ;

numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps incremArcLen ] ; 


% --------------------------------------------------------------------------------------------------
