% ------------------------------------
% TEST example pendulum
% ------------------------------------

clear all, close all

% --- general data ---
inputONSASversion = '0.1.10';

dirOnsas = [ pwd '/..' ] ;
problemName = 'rightAngleCantileverBeam' ;
% ------------------------------------

Es = 1 ;
nu = 0    ;
A  = 1 ;
l0 = 1     ;

% 2*10 = rho*A*l0 ;
mconc = 10;
rho    = 2*mconc / ( A * l0 ) ;

omega = sqrt(Es*A/l0/mconc) 
nodalDamping = 0 ;

perio = 1/ ( omega/(2*pi) )

%~ stop
% method
timeIncr   =  5    ;
finalTime  = 20*perio     
nLoadSteps = finalTime/timeIncr ;

%~ alphaHHT = 0 ;
alphaHHT = -0.05 ;


% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-10        ;
stopTolIts    = 30          ;
% ------------------------------------



% --- structural properties ---
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu rho] ;

secGeomProps = [ A 0 0 0 ] ;

nodalSprings = [  1  inf  0  inf  0  inf 0  ; ...
2  0  0  inf  0  inf 0  ...
               ];

Nodes = [    0  0  0   ; ...
             l0 0  0 ] ;
             %~ 0  0  -l0 ] ;

Conec = [ 1 2 0 0 1 1 1 ] ; 


% -------------------
nodalConstantLoads   = [ 2  0  0  0  0  -rho*A*l0*0.5*9.8  0 ];
% or
%~ nodalVariableLoads   = [ 2  0  0  0  0  0  0 ];
%~ userLoadsFilename = 'myLoadSpringMass' ;

%~ nodalDamping = cres ;
% -------------------

controlDofInfo = [ 2 1 +1 ] ;
% ------------------------------

% ------------------------------
% analysis parameters
dynamicAnalysisBoolean   = 1 ; 

% initial conditions
u0    = l0*0.1 ;
%~ udot0 = 7.725   ; % m/sec

nonHomogeneousInitialCondU0    = [ 2 1 u0    ] ;
%~ nonHomogeneousInitialCondUdot0 = [ 2 1 udot0 ] ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;

plotParamsVector = [0 ];
%~ plotParamsVector = [ 3 30 ];
printflag = 0 ;

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;


%~ angs = asin( (l0+controlDisps) ./ l0 ) * 180 / pi ;


lw  = 2   ; ms  = 5.5 ;
lw2 = 3.2 ; ms2 = 23 ;
plotfontsize = 22 ;

figure
plot(timesVec, controlDisps,'b--o','linewidth',lw,'markersize',ms);
