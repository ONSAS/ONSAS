% ------------------------------------------------------------------------------
% example Right-angle cantilever
% ------------------------------------------------------------------------------

clear all, close all
dirOnsas = [ pwd '/..' ] ;   problemName = 'rightAngleCantileverBeam' ;

% ------------------------------------
E   =  1e6  ;
nu  = -0.5  ;
A   =  1    ;
I   =  1e-3 ;
L   =  10   ;
rho =  1    ;

% inconsistent
J   = I ;
global Jrho
Jrho = diag( [ 20; 10; 10 ] ) ;

materialsParams = {[ rho 1 E nu ]} ;

crossSecsParams = [ A I I J ] ;

% method
timeIncr   =  0.025    ;
finalTime  = 5    
nLoadSteps = finalTime/timeIncr ;

alphaHHT = 0 ;
%~ alphaHHT = -0.05 ;

% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-10        ;
stopTolIts    = 30          ;
% ------------------------------------

nodalSprings = [  1  inf  0  inf  0  inf 0  ; ...
2  0  0  inf  0  inf 0  ...
               ];

nElemsPerBeam = 4 ;
Nodes = [ zeros(nElemsPerBeam+1,1)       linspace(0,L,nElemsPerBeam+1)'  zeros(nElemsPerBeam+1,1) ; ...
          -linspace(0,L,nElemsPerBeam+1)(2:end)' zeros(nElemsPerBeam,1)       zeros(nElemsPerBeam,1) ] ;

Conec = [ (1:(2*nElemsPerBeam))' (2:(2*nElemsPerBeam+1))' zeros(2*nElemsPerBeam,2) ones(2*nElemsPerBeam,2) 2*ones(2*nElemsPerBeam,1) ] ; 

% -------------------
nodalVariableLoads   = [ nElemsPerBeam+1  0  0  0  0  1  0 ];

controlDofs = [ nElemsPerBeam 5 +1 ] ;

loadFactorsFunc = @(t) t ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;

plotParamsVector = [0 ];
%~ plotParamsVector = [ 3 30 ];
printFlag = 0 ;

acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;

%~ lw  = 2   ; ms  = 5.5 ;
%~ lw2 = 3.2 ; ms2 = 23 ;
%~ plotfontsize = 22 ;

%~ figure
%~ plot(timesVec, controlDisps,'b--o','linewidth',lw,'markersize',ms);
