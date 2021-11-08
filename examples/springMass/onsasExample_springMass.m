% ------------------------------------
% springmass example
% Notation and analytical based on chapter 3 from
% Ray W. Clough and Joseph Penzien, Dynamics of Structures, Third Edition, 2003
% ------------------------------------

close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );

otherParams.problemName = 'springMass' ;
otherParams.plotsFormat = 'vtk' ;

% spring mass system
k        = 39.47 ;
p0       = 40    ;
c        = 0.1   ;
m        = 1     ;
omegaBar = 4*sqrt(k/m) ;
p0       = 40    ;
u0       = 0.0   ; % initial displacement

% parameters for truss model
l   = 1   ;
A   = 0.1 ;
rho = m * 2 / ( A * l ) ;
E   = k * l /   A       ;

omegaN = sqrt( k / m );
xi     = c / m  / ( 2 * omegaN ) ;
nodalDamping = c ;

freq   = omegaN / (2*pi)      ;
TN     = 2*pi / omegaN        ;
dtCrit = TN / pi              ;

% numerical method params
stopTolDeltau = 1e-10           ;
stopTolForces = 1e-10           ;
stopTolIts    = 30              ;
alphaHHT      = 0;
% ------------------------------------


materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ E 0 ] ;
materials(1).density  = rho ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemTypeGeometry = [2 sqrt(A) sqrt(A) ] ;
elements(2).elemTypeParams = 0 ;

boundaryConds(1).imposDispDofs =  [ 1 3 5 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;

boundaryConds(2).imposDispDofs =  [ 3 5 ] ;
boundaryConds(2).imposDispVals =  [ 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global'                  ;
boundaryConds(2).loadsTimeFact = @(t) p0*sin( omegaBar*t )                    ;
boundaryConds(2).loadsBaseVals = [ 1 0 0 0 0 0 ] ;

% initial conditions
initialConds.nonHomogeneousInitialCondU0    = [ 2 1 u0    ] ;

mesh.nodesCoords = [  0  0  0 ; ...
                      l  0  0 ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0   2   ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 2   ] ;


analysisSettings.methodName    = 'newmark' ;
%md and the following parameters correspond to the iterative numerical analysis settings
analysisSettings.deltaT        =   0.005  ;
analysisSettings.finalTime      =   1*2*pi/omegaN   ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   10   ;
analysisSettings.alphaNM      =   0.25   ;
analysisSettings.deltaNM      =   0.5   ;

[matUsNewmark, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

if (c == 0) && (p0 == 0) % free undamped
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
     ( A * cos( omegaD   * t ) + B  * sin( omegaD   * t ) ) .* exp( -xi * omegaN * t ) ...
    + G1 * cos( omegaBar * t ) + G2 * sin( omegaBar * t ) ;
    analyticCheckTolerance = 5e-2 ;
end

times = 0:analysisSettings.deltaT:(analysisSettings.finalTime+analysisSettings.deltaT) ;

valsAnaly = analyticFunc(times) ;
valsNewmark = matUsNewmark(6+1,:) ;

analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.alphaHHT      =   0   ;

[matUsHHT, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

valsHHT = matUsHHT(6+1,:) ;

verifBooleanNewmark =  ( ( norm( valsAnaly - valsNewmark    ) / norm( valsAnaly   ) ) <  analyticCheckTolerance )
verifBooleanHHT     =  ( ( norm( valsAnaly - valsNewmark     ) / norm( valsAnaly  ) ) <  analyticCheckTolerance )

figure
hold on, grid on
plot(valsAnaly,'b-x')
plot(valsNewmark,'r-o')
plot(valsHHT,'g-s')
