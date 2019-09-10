% --------------------------------------------------------------------------------------------------
% Test problem: rectangular cross-section cantilever beam submitted to nodal external loads at the free end.
% --------------------------------------------------------------------------------------------------

E   = 210e9   ;
nu  = 0.3     ;
l   = 2       ; % reference length
byL = 0.05    ;
bzL = 0.05    ;

inputONSASversion = '0.1.8';

problemName = 'cantilever' ;

% constitutive properties

hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu] ;

% geometrical properties
Iyy = byL * bzL^3 / 12 ;
Izz = bzL * byL^3 / 12 ;

secGeomProps = [ byL*bzL   Iyy   Izz   1 ] ;

nelems = 5 ;

Nodes = [  linspace(0,l,nelems+1)' zeros(nelems+1,2) ] ;

% in global system of coordinates
nodalSprings = [ 1  inf  inf  inf  inf  inf  inf ] ;

Conec = [ (1:(nelems))' (2:(nelems+1))' zeros(nelems,2) ones(nelems,2) 2*ones(nelems,1) ] ;

nodalVariableLoads   = [ nelems+1  -1 0   0 0   0 1e-3 ] ;

% analysis parameters
nonLinearAnalysisBoolean = 1 ; 
dynamicAnalysisBoolean   = 0 ; 

% [ node nodaldof scalefactor(positive or negative) ]
controlDofInfo = [ nelems+1 3 1 ] ;

printflag = 2 ;

plotParamsVector = [2 4  ];
%~ plotParamsVector = [2 50];

sectPar = [ 12 byL bzL ];

% analytical solution 
analyticSolFlag = 0 ;

reportBoolean = 1 ;

Pcr = min( [ E * pi^2 * Iyy / ( 2*l)^2, E * pi^2 * Izz / ( 2*l)^2 ] ) ;

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-8 ;
stopTolForces    = 1.0e-8 ;
targetLoadFactr  = Pcr*1.2    ;
nLoadSteps       = 360    ;
%~ targetLoadFactr  = Pcr*0.5    ;
%~ nLoadSteps       = 20    ;
incremArcLen     = 0.2e-3    ;

%~ numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            %~ targetLoadFactr nLoadSteps incremArcLen ] ; 
numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps ] ; 
% --------------------------------------------------------------------------------------------------
