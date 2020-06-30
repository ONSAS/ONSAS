% TEST example springmass
% ------------------------------------

% auxiliar numerical data
Es = 200e9 ;
b = .1 ;
A  = b*b ;
l0 = 2 ;

rhoprob =  8    ; % kg/m3
nu      =  0    ;

% method
timeIncr   =  0.01    ;
finalTime  =  5       ;
nLoadSteps = finalTime/timeIncr ;

DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;

% tolerances
stopTolDeltau = 1e-10           ; 
stopTolForces = 1e-10           ;
stopTolIts    = 30            ;
% ------------------------------------


% --- general data ---
inputONSASversion = '0.1.9';

problemName = 'rigidPendulum' ;
% ------------------------------------

% --- structural properties ---
rho = rhoprob ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu rho] ;


secGeomProps = [ A 0 0 0 ] ;

nodalSprings = [ 1  inf  0  inf  0  inf 0   ...
               ];

Nodes = [    0  0  0 ; ...
            l0  0  0 ] ;


Conec = [ 1 2 0 0 1 1 1 ] ; 

%~ loadFactorsFunc = @(t) p0 *sin( omegaBar*t) ; 

% -------------------
%~ nodalVariableLoads   = [ 2  1  0  0  0  0  0 ];
% or
%~ nodalVariableLoads   = [ 2  0  0  0  0  0  0 ];
nodalConstantLoads   = [ 2  0  0  0  0  -rhoprob*l0*A*.5*9.81  0 ];
% -------------------

controlDofInfo = [ 2 5 -1 ] ;
% ------------------------------

% ------------------------------
% analysis parameters
dynamicAnalysisBoolean   = 1 ; 

% initial conditions
u0   = 0;
udot0= 0;

%~ nonHomogeneousInitialCondU0 = [ 2 1 u0 ] ;

%~ nonHomogeneousInitialCondUdot0 =[ 2 1 udot0 ] ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

plotParamsVector = [3 101];
printflag = 2 ;

sectPar = [ 12 b b ] ;

analyticSolFlag = 0 ;
% ------------------------------------------------
