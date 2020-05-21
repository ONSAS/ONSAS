% ------------------------------------
% TEST example springmass
% ------------------------------------

clear all, close all

% auxiliar numerical data
Es = .5 ;
A  = 1 ;
l0 = 4 ;

rhoprob =  1    ;
ceda    =  0.02 ;
nu      =  0    ;

% sForce data
omegaBar = 2   ;
p0       = 0.1 ;

%~ %Parameters
%~ targetLoadFactr=1.5;

rho    = rhoprob ;
kres   = (Es*A/l0)            ;
mres   = (rho * A*l0 /2)      ;
omegaN = sqrt( kres / mres );
cres   = 2*ceda*omegaN*mres   ;

nodalDamping = cres;

freq   = omegaN / (2*pi)      ;
TN     = 2*pi / omegaN        ;
dtCrit = TN / pi              ;

% method
timeIncr   =  0.001 * dtCrit    ;
finalTime  =  14                ;
nLoadSteps = finalTime/timeIncr ;
DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;

% tolerances
stopTolDeltau = 1e-12           ; 
stopTolForces = 1e-12           ;
stopTolIts    = 1000            ;
% ------------------------------------


% --- general data ---
inputONSASversion = '0.1.9';

dirOnsas = [ pwd '/..' ] ;
problemName = 'springMass_NM' ;
% ------------------------------------

% --- structural properties ---
rho = rhoprob ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu rho] ;


secGeomProps = [ A 0 0 0 ] ;

nodalSprings = [ 1  inf  0  inf  0  inf 0 ; ...
                 2    0  0  inf  0  inf 0   ...
               ];

Nodes = [    0  0  0 ; ...
            l0  0  0 ] ;

auxelemtype = 1 ;
Conec = [ 1 2 0 0 1 1 auxelemtype ] ; 

loadFactorsFunc = @(t) p0 *sin( omegaBar*t) ; 

% -------------------
%~ nodalVariableLoads   = [ 2  1  0  0  0  0  0 ];
% or
nodalVariableLoads   = [ 2  0  0  0  0  0  0 ];
userLoadsFilename = 'myLoadSpringMass' ;
% -------------------

controlDofInfo = [ 2 1 +1 ] ;
% ------------------------------

% ------------------------------
% analysis parameters
dynamicAnalysisBoolean   = 1 ; 

% initial conditions
u0    = 0.2 ;
udot0 = 0   ;

nonHomogeneousInitialCondU0    = [ 2 1 u0    ] ;
nonHomogeneousInitialCondUdot0 = [ 2 1 udot0 ] ;

numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

plotParamsVector = [2 5 ];
printflag = 2 ;

kres = (Es*A/l0) ;
mres = (rho * A*l0 /2);
omega = sqrt( kres / mres ) ;

if u0 < l0
  omegaReal = omegaN*sqrt(1-ceda^2);
  beta= omegaBar/omegaN ;
  G1=(p0/kres) * (-2*ceda*beta/((1-beta^2)^2+(2*ceda*beta)^2));
  G2=(p0/kres) * ((1-beta^2)/((1-beta^2)^2+(2*ceda*beta)^2)) ;
  A=u0 -G1;
  B=(ceda*omegaN*A-omegaBar*G2)/(omegaReal);

  analyticSolFlag = 1 ;
  analyticFunc = @(t) (A*cos(omegaReal*t)+B*sin(omegaReal*t)).*exp(-ceda*omegaN*t)+ G1*cos(omegaBar*t)+G2*sin(omegaBar*t) ;
  analyticCheckTolerance = 5e-2 ;
else
  error('this analytical solution is not valid for this u0 and l0');
end
% ------------------------------------------------

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;
