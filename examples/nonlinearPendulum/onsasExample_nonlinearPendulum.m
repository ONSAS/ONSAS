
clear all, close all

dirOnsas = [ pwd '/../../src' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path

problemName = 'nonlinearPendulum' ;
% ------------------------------------

% -- scalar params -----
Es  = 10e11  ; nu  = 0      ;
A   = 0.1    ; l0  = 3.0443 ;
m   = 10     ; g   = 9.81   ;

nodalDispDamping = 0 ;
rho = 2*m / ( A * l0 ) ;

% ----- geometry -----------------
Nodes = [     0 0  l0   ; ...
             l0 0  l0 ] ;

Conec = { [ 0 1 0 0 1  1   ] ; ...
          [ 0 1 1 0 2  2   ] ; ...
          [ 1 2 0 1 0  1 2 ] } ;


%  --- MELCS params ---

materialsParams = { [ rho 3 Es nu ] } ;

elementsParams = { 1 ; [ 2 0 ] } ;

loadsParams  = {[ 1 0  0  0  0  0  -rho*A*l0*0.5*g  0 ]};

crossSecsParams = { [2 sqrt(A) sqrt(A) ] } ;

springsParams = {[  inf  0  inf  0  inf 0 ] ; ...
                 [  0    0  inf  0  0   0 ] };

% -------------------

% method
timeIncr   =  0.05    ;
finalTime  =  4.13*8 ;

alphaHHT = -0.05 ;
%~ alphaHHT = 0 ;

% tolerances
stopTolDeltau = 0    ;
stopTolForces = 1e-6 ;
stopTolIts    = 30   ;
% ------------------------------------

controlDofs = [ 2 1 1 ] ;

% analysis parameters
numericalMethodParams = [ 4 ...
  timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;

plotParamsVector = [ 3 ];

% run ONSAS
ONSAS ;

controlDispsA = controlDisps ;

% ------------------------------------

problemName = 'simplePendulumTrussHHTSelfWeightDamp' ;

Conec = { [ 0 1 0 0 1  1   ] ; ...
          [ 0 1 0 0 2  2   ] ; ...
          [ 1 2 0 1 0  1 2 ] } ;

nodalDispDamping = .5 ;

loadsParams  = {};
booleanSelfWeightZ = 1 ;

% run ONSAS
ONSAS ;

controlDispsB = controlDisps ;

% ------------------------------------
%uNum = PenduloNL_HHT( l0, A, Es, m, nodalDispDamping, g, timeIncr, -alphaHHT, finalTime, stopTolForces, stopTolIts, 1, 0 ) ;

figure, hold on, grid on, MS = 10; LW = 1.5 ;
plot(timesVec, controlDispsA, 'b-o','markersize',MS,'linewidth',LW)
plot(timesVec, controlDispsB, 'g-x','markersize',MS,'linewidth',LW)
plot(timesVec, uNum(1,:)-l0 ,'r--','markersize',MS,'linewidth',LW)
xlabel('time (s)'), ylabel('control displacement')
legend('onsasA', 'onsasB', 'semi-analytic')
print('output.png','-dpng')
