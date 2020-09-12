% ------------------------------------
% springmass example
%
% Notation and analytical based on chapter 3 from
% Ray W. Clough and Joseph Penzien, Dynamics of Structures, Third Edition, 2003
% ------------------------------------

function onsasExample_springMass( dirOnsas, scalarParams )

close all

if nargin == 0
  dirOnsas = [ pwd '/..' ] ;
end
addpath( dirOnsas ) ;

problemName = 'springMass' ;

% scalar parameters

% spring mass system 
if nargin <= 1
  k        = 39.47 ;
  c        = 0.1   ;
  m        = 1     ;
  omegaBar = 4*sqrt(k/m) ;
  p0       = 40     ;
  u0       = 0.0   ; % initial displacement
end

% parameters for truss model
l   = 1   ;
A   = 0.1 ;
rho = m * 2 / ( A * l ) ;
E   = k * l /   A       ;
nodalDispDamping = c ;

omegaN = sqrt( k / m );
xi     = c / m  / ( 2 * omegaN ) ;
nodalDamping = c ;

freq   = omegaN / (2*pi)      ;
TN     = 2*pi / omegaN        ;
dtCrit = TN / pi              

% numerical method params
timeIncr      =  0.01 ;
finalTime     = 2*2*pi/omegaN                ;
stopTolDeltau = 1e-10           ; 
stopTolForces = 1e-10           ;
stopTolIts    = 30              ;
alphaHHT      = 0;
% ------------------------------------


% --- structural properties ---
materialsParams = {[ rho 1 E 0 ]} ;

elementsParams = { 1 ; [2 0]} ;

loadsParams    = { [ 1 1    1 0 0 0 0 0 ] } ;
loadFactorsFunc = @(t) p0 *sin( omegaBar*t ) ; 

crossSecsParams = {[ 2 sqrt(A) sqrt(A) ]} ;

springsParams  = {[ inf  0  inf  0  inf 0 ] ; ...
                  [ 0    0  inf  0  inf 0 ] } ;

Nodes = [    0  0  0 ; ...
            l  0  0 ] ;

Conec = { [ 0 1 0 0 1  1 ] ;
          [ 0 1 1 0 2  2 ] ;
          [ 1 2 0 1 0  1 2 ] } ; 

controlDofs = [ 2 1 +1 ] ;
% ------------------------------

% initial conditions
nonHomogeneousInitialCondU0    = [ 2 1 u0    ] ;
%~ nonHomogeneousInitialCondUdot0 = [ 2 1 udot0 ] ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;

plotParamsVector = [3 ];

if c == 0 && p0 == 0 % free undamped

  analyticSolFlag = 1 ;
  analyticFunc = @(t)   (   u0 * cos( omegaN * t )  ) ;
  analyticCheckTolerance = 2e-1 ;

else
  beta   = omegaBar / omegaN ;
  omegaD = omegaN * sqrt( 1-xi^2 ) ;

  G1 = (p0/k) * ( -2 * xi * beta   ) / ( ( 1 - beta^2 )^2 + ( 2 * xi * beta )^2 ) ;
  G2 = (p0/k) * (  1      - beta^2 ) / ( ( 1 - beta^2 )^2 + ( 2 * xi * beta )^2 ) ;
  if u0 < l
    A  = u0 - G1 ;
    B  =  (xi*omegaN*A - omegaBar*G2 ) / (omegaD);
  else
    error('this analytical solution is not valid for this u0 and l0');
  end
  
  analyticSolFlag = 1 ;
  analyticFunc = @(t) ...
    ( A * cos( omegaD * t ) + B * sin( omegaD * t ) ) .* exp( -xi * omegaN * t ) ...
    + G1 * cos( omegaBar * t ) + G2 * sin( omegaBar * t ) ;
    analyticCheckTolerance = 5e-2 ;
end
% ------------------------------------------------


global flagOutputMatrices
flagOutputMatrices = 1 ;

storeBoolean = 1 ;

ONSAS

load ../../output.mat 

save -mat ../../output.mat KT dampingMat massMat matUs loadFactors

