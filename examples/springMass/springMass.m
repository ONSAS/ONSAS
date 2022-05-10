%md# Spring-mass system example
%md
%md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS.m/blob/master/examples/springMass/springMass.m)
%md
%mdIn this example a simple spring-mass system is considered. The notation and the analytical solution are based on chapter 3 from
%mdthe book Dynamics of Structures by Ray W. Clough and Joseph Penzien, Third Edition, 2003.
%md
%md First, the path to the ONSAS folder is added and scalar parameters are set:
% add path
close all, clear all; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters for spring-mass system
k        = 39.47 ;
c        = 0.01  ;
m        = 1     ;
p0       = 40    ; % amplitude of applied load
u0       = 0.1   ; % initial displacement
%md then the parameters for the equivalent truss model are computed:
l   = 1   ;
A   = 0.1 ;
rho = m * 2 / ( A * l ) ;
E   = k * l /   A       ;
%md
omegaN       = sqrt( k / m );
omegaBar     = 4*omegaN ;
xi           = c / m  / ( 2 * omegaN ) ;
freq         = omegaN / (2*pi)      ;
TN           = 2*pi / omegaN        ;
dtCrit       = TN / pi              ;
%md
%md## Analytic solution
%md The analytical solution of the problem is:
%md```math 
%md  u(t) =
%md     ( A_c \cos( \omega_D  t ) + B \sin( \omega_D t ) ) e^{ -\xi \omega_N t } +
%md    G_1  \cos( \bar{\omega} t ) + G_2 \sin( \bar{\omega} t )
%md``` 
%md where the parameters are computed depending on the damping value $c$ and the load $p$:
if (c == 0) && (p0 == 0) % free undamped
  myAnalyticFunc = @(t)   (   u0 * cos( omegaN * t )  ) ;
  analyticCheckTolerance = 2e-1 ;
else
  beta   = omegaBar / omegaN ;  omegaD = omegaN * sqrt( 1-xi^2 ) ; %forced and damped
  G1 = (p0/k) * ( -2 * xi * beta   ) / ( ( 1 - beta^2 )^2 + ( 2 * xi * beta )^2 ) ;
  G2 = (p0/k) * (  1      - beta^2 ) / ( ( 1 - beta^2 )^2 + ( 2 * xi * beta )^2 ) ;
  if u0 < l
    Ac = u0 - G1 ;
    B  =  (xi*omegaN*Ac - omegaBar*G2 ) / (omegaD);
  else
    error('this analytical solution is not valid for this u0 and l0');
  end
  myAnalyticFunc = @(t) ...
     ( Ac * cos( omegaD   * t ) + B  * sin( omegaD   * t ) ) .* exp( -xi * omegaN * t ) ...
    + G1  * cos( omegaBar * t ) + G2 * sin( omegaBar * t ) ;
    analyticCheckTolerance = 5e-2 ;
end
%md 
%md## Numerical solution
%md---------------------
%md The solution is used to validate two different dynamics methods: Newmark and $\alpha-HHT$ considering truss elements with nodal or lumped mass.   
%md
%md### Numerical case 1: truss element model with Newmark method
%md
%md#### Materials
%md
%md A material for the truss is set in order to reproduce a constant matrix at the end node of the truss equal to $m$:
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ E 0 ] ;
materials(1).density         = rho ;
%md
%md### Elements
%md
%md A material for the truss is set in order to reproduce a constant matrix at the end node of the truss equal to $m$:
elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemCrossSecParams = {'circle', [sqrt(4*A/pi) ] } ;
elements(2).massMatType = 'lumped' ;
%md
%md### Boundary conditions
%md
%md The fixed node at the beginning of the truss is set as:
boundaryConds(1).imposDispDofs =  [ 1 3 5 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;
%
%md The end node constrain allows the truss to move in $x$ so:
boundaryConds(2).imposDispDofs =  [ 3 5 ] ;
boundaryConds(2).imposDispVals =  [ 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global'                  ;
boundaryConds(2).loadsTimeFact = @(t) p0*sin( omegaBar*t )                    ;
boundaryConds(2).loadsBaseVals = [ 1 0 0 0 0 0 ] ;
%md
%md### Initial conditions
%md One initial condition is set in $x$ direction and equal to $u_0$:
initialConds.nonHomogeneousUDofs    = [ 1  ] ;
initialConds.nonHomogeneousUVals    = [ u0 ] ;
%md
%md### mesh
%md Only two nodes are considered so the nodes matrix is:
mesh.nodesCoords = [  0  0  0 ; ...
                      l  0  0 ] ;
mesh.conecCell = { } ;
% The first node has no material, the first element of the _elements_ struct, which is `'node'` also the first boundary condition (fixed) and no initial condition is set.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
% The second node has no material, the first element of the _elements_ struct, which is `'node'` also the second boundary condition (x disp free) and the first initial condition ($u_0$) is set.
mesh.conecCell{ 2, 1 } = [ 0 1 2 1   2   ] ;
% Only one element is considered with the first material and the second element of respective structs
mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 2   ] ;
%md
%md and the following parameters correspond to the iterative numerical analysis settings
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   0.005  ;
analysisSettings.finalTime     =   1.2*TN   ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   10   ;
%md
otherParams.problemName = 'springMass' ;
%md The boolean exportFirstMatrices allows to extract global tangent and mass matrix at the initial configuration using: 
global exportFirstMatrices
exportFirstMatrices = true            ;
%md
[matUsNewmark, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
times = 0:analysisSettings.deltaT:(analysisSettings.finalTime+analysisSettings.deltaT) ;
%md
%md
%md### Numerical case 2: nodal mass model with $\alpha$-HHT method
%md
%md The nodalMass field allows to add lumped matrices to a node, since this field is used then the equivalent $\rho$ of the material 1 aforementioned now is set to 0. Although an equal $m$ mass is considered for $u_x$ $u_y$ and $u_z$ at the end node, so: 
materials(1).density   = 0 ;
materials(2).nodalMass = [m m m] ;
%md now the initial condition is added to the end node with the second material:
mesh.conecCell{ 2, 1 } = [ 2 1 2 1   2   ] ;
%md
%md The $\alpha_{HHT}$ method with $\alpha=0$ is the same as Newmark, this is employed to validate results of both methods, as consequence 
analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.alphaHHT      =   0   ;
%md
[matUsHHT, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md 
%md## Verification
%md---------------------
%md The numerical displacements of the end node in $x$ are extracted for different models tested:  
valsNewmark = matUsNewmark(6+1,:) ;
valsHHT  = matUsHHT(6+1,:) ;
%md The analytical solution is evaluated:
valsAnaly = myAnalyticFunc(times) ;
%md The boolean to validate the implementation is evaluated such as:
verifBooleanNewmark =  ( ( norm( valsAnaly - valsNewmark ) / norm( valsAnaly ) ) <  analyticCheckTolerance ) ;
verifBooleanHHT     =  ( ( norm( valsAnaly - valsHHT     ) / norm( valsAnaly ) ) <  analyticCheckTolerance ) ;
verifBoolean = verifBooleanHHT && verifBooleanNewmark ;
%md
%md## Plot verification
%md---------------------figure
hold on, grid on
plot(times, valsAnaly,'b-x')
plot(times, valsNewmark,'r-o')
plot(times, valsHHT,'g-s')
labx = xlabel('time');   laby = ylabel('\lambda(t)') ;
legend( 'analytic', 'truss-Newmark','nodalMass-HHT', 'location','northoutside')
print('output/springMassCheck.png','-dpng')
if exist('../../docs/src/assets/')==7
  % printing plot also to docs directory
  disp('printing plot also to docs directory')
  print('../../docs/src/assets/springMassCheck.png','-dpng')
end
%md
%md```@raw html
%md<img src="../../assets/springMassCheck.png" alt="plot check" width="500"/>
%md```
%md
