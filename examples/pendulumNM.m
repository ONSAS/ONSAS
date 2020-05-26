% ------------------------------------
% TEST example pendulum
% ------------------------------------

clear all, close all

% --- general data ---
inputONSASversion = '0.1.10';

dirOnsas = [ pwd '/..' ] ;
problemName = 'pendulumNM' ;
% ------------------------------------


Es = 10e11 ;
nu      =  0    ;

A  = 0.1 ;
l0 = 3.0443     ;

% 2*10 = rho*A*l0 ;

rho    = 2*10 / ( A * l0 ) ;

nodalDamping = 0.000 ;
%~ nodalMass    = 

% method
timeIncr   =  0.025    ;
finalTime  =  8                ;
%~ finalTime  =  1                ;
nLoadSteps = finalTime/timeIncr ;
DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;

% tolerances
stopTolDeltau = 1e-12           ; 
stopTolForces = 1e-12           ;
stopTolIts    = 30              ;
% ------------------------------------



% --- structural properties ---
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu rho] ;

secGeomProps = [ A 0 0 0 ] ;

nodalSprings = [ 1  inf  0  inf  0  inf 0  ...
               ];

Nodes = [    0  0  0   ; ...
             0  0  -l0 ] ;

Conec = [ 1 2 0 0 1 1 1 ] ; 


% -------------------
nodalConstantLoads   = [ 2  0  0  0  0  -rho*A*l0*0.5*9.8  0 ];
% or
%~ nodalVariableLoads   = [ 2  0  0  0  0  0  0 ];
%~ userLoadsFilename = 'myLoadSpringMass' ;

%~ nodalDamping = cres ;
% -------------------

controlDofInfo = [ 2 1 1 ] ;
% ------------------------------

% ------------------------------
% analysis parameters
dynamicAnalysisBoolean   = 1 ; 

% initial conditions
%~ u0    = 0 ;
udot0 = 7.725   ; % m/sec

%~ nonHomogeneousInitialCondU0    = [ 2 1 u0    ] ;
nonHomogeneousInitialCondUdot0 = [ 2 1 udot0 ] ;

numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

%~ plotParamsVector = [2 5 ];
plotParamsVector = [ 3 20 ];
printflag = 0 ;

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;
