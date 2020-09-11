% ------------------------------------
% TEST example springmass
%
% From chapter 2 Ray W. Clough and Joseph Penzien, Dynamics of Structures, Third Edition, 2003
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
  c        = 0   ;
  m        = 1     ;
  omegaBar = 2*pi  ;
  p0       = 0     ;
  u0       = 0.2   ; % initial displacement
end

% parameters for truss model
l   = 1   ;
A   = 0.1 ;
rho = m * 2 / ( A * l ) ;
E   = k * l /   A       ;

omegaN = sqrt( k / m );
nodalDamping = c ;

freq   = omegaN / (2*pi)      ;
TN     = 2*pi / omegaN        ;
dtCrit = TN / pi              ;

% numerical method params
timeIncr      =  0.01 ;
finalTime     = 2*pi/omegaN                ;
stopTolDeltau = 1e-15           ; 
stopTolForces = 1e-12           ;
stopTolIts    = 30              ;
alphaHHT      = 0;
% ------------------------------------


% --- structural properties ---
materialsParams = {[ rho 1 E 0 ]} ;

elementsParams = { 1 ; [2 0]} ;

loadsParams    = { [ 1 1    1 0 0 0 0 0 ] } ;
loadFactorsFunc = @(t) p0 *sin( omegaBar*t ) ; 

crossSecsParams = {[ 3 sqrt(A*4/pi) ]} ;

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

end
%~ case 2 % forced
  %~ if u0 < l0
    %~ omegaReal = omegaN * sqrt( 1-ceda^2 ) ;
    %~ beta      = omegaBar/omegaN ;
    %~ G1        = (p0/kres) * ( -2 * ceda * beta / ( ( 1 - beta^2 )^2 + ( 2 * ceda * beta )^2 ) ) ;
    %~ G2        = (p0/kres) * ( ( 1 - beta^2 )   / ( ( 1 - beta^2 )^2 + ( 2 * ceda * beta )^2 ) ) ;
    %~ A         = u0 - G1 ;
    %~ B         =  (ceda*omegaN*A - omegaBar*G2 ) / (omegaReal);
  
    %~ analyticSolFlag = 1 ;
    %~ analyticFunc = @(t) ...
      %~ ( A*cos(omegaReal*t)+B*sin(omegaReal*t)).* exp( -ceda * omegaN * t ) ...
      %~ + G1 * cos( omegaBar * t ) + G2 * sin( omegaBar * t ) ;
    %~ analyticCheckTolerance = 5e-2 ;
  %~ else
    %~ error('this analytical solution is not valid for this u0 and l0');
  %~ end
%~ end
% ------------------------------------------------

ONSAS


figure
