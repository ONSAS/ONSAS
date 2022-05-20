%md# Spring-mass-damper system example
%md
%md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS.m/blob/master/examples/springMass/springMass.m)
%md
%md## ODE System definition
%mdIn this example a simple spring-mass-damper system is considered. This simple problem is used to validate different implementations of truss and frame elements submitted to dynamic loads.
%md
%md```@raw html
%md<img src="../../assets/springMassSystem.svg" alt="spring-mass diagram" width="800"/>
%md```
%md## Analytic solution
%mdThe analytical solution is based on chapter 3 from _Dynamics of Structures_ by Ray W. Clough and Joseph Penzien, Third Edition, 2003.
%mdThe analytical solution of the problem is given by:
%md```math
%md  u(t) =
%md     \left( A_c \cos( \omega_D  t ) + B \sin( \omega_D t )  \right) e^{ -\xi \omega_N t } +
%md    G_1  \cos( \bar{\omega} t ) + G_2 \sin( \bar{\omega} t )
%md```
%md where the notation and the parameters interpretation can be seen in the cited reference.
%md
%mdA set of numerical parameters must be defined to compute the analytic solution.
%md
%md We start as all models, clearing the workspace and adding the ONSAS path to the work path.
% clear workspace and add path
close all, clear all; addpath( genpath( [ pwd '/../../src'] ) );
%md The following numeric parameters are considered.
% scalar parameters for spring-mass system
k    = 39.47 ; % spring constant
c    = 0.5   ; % damping parameter
m    = 1     ; % mass of the system
p0   = 40    ; % amplitude of applied load
u0   = 0.1   ; % initial displacement
du0  = 0.0   ; % initial velocity
%md
%md Then other parameters are computed:
%md
omegaN       = sqrt( k / m )           ; % the natural frequency
xi           = c / m  / ( 2 * omegaN ) ;
freq         = omegaN / (2*pi)         ;
TN           = 2*pi / omegaN           ;
dtCrit       = TN / pi                 ;
%md
%md The frequency of the sinusoidal external force is set as:
omegaBar     = 4*omegaN ;
%md
%md The analytic solution can be computed for specific cases as follows:
if (c == 0) && (p0 == 0) % free undamped solution
  myAnalyticFunc = @(t)   (   u0 * cos( omegaN * t )  ) ;
else                     % other cases solution
  beta   = omegaBar / omegaN ;  omegaD = omegaN * sqrt( 1-xi^2 ) ; %forced and damped
  G1 = (p0/k) * ( -2 * xi * beta   ) / ( ( 1 - beta^2 )^2 + ( 2 * xi * beta )^2 ) ;
  G2 = (p0/k) * (  1      - beta^2 ) / ( ( 1 - beta^2 )^2 + ( 2 * xi * beta )^2 ) ;
  if u0 < l
    Ac = u0 - G1 ;  B = (xi*omegaN*Ac - omegaBar*G2 ) / (omegaD);
  else
    error('this analytical solution is not valid for this u0 and l0');
  end
  myAnalyticFunc = @(t) ...
     ( Ac * cos( omegaD   * t ) + B  * sin( omegaD   * t ) ) .* exp( -xi * omegaN * t ) ...
    + G1  * cos( omegaBar * t ) + G2 * sin( omegaBar * t ) ;
end
%md
%md## Numerical solutions
%md
%md The analytic solution is used to validate two numerical solution approaches using different structural physical models, governed by the same ODE.
%md
%md### Numerical case 1: truss element model with Newmark method
%mdIn this case, a truss element is considered, as shown in the figure, with Young modulus, cross-section, area, mass, nodal damping and length corresponding to the parameters considered for the spring-mass-damper system 
%md```@raw html
%md<img src="../../assets/springMassSystem.svg" alt="spring-mass diagram" width="800"/>
%md```
%md
%md The scalar parameters for the equivalent truss model are:
l   = 1                 ;
A   = 0.1               ;
rho = m * 2 / ( A * l ) ;
E   = k * l /   A       ;
%md where the material of the truss was selected to set a mass $m$ at the node $2$.
%md
%md#### Materials
%md
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ E 0 ]          ;
materials(1).density         = rho              ;
%md
%md#### Elements
%md
%md In this case only `'node'` and  `'truss'` elements are considered and the lumped inertial formulation is set for the truss element:
elements(1).elemType = 'node'                                 ;
elements(2).elemType = 'truss'                                ;
elements(2).elemCrossSecParams = {'circle', [sqrt(4*A/pi) ] } ;
elements(2).massMatType = 'lumped'                            ;
%md
%md#### Boundary conditions
%md
%md The node $1$ is fixed, so the boundary condition set is:
boundaryConds(1).imposDispDofs =  [ 1 3 5 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;
%md The node $2$ allows the truss to move in $x$ so the boundary condition set is:
boundaryConds(2).imposDispDofs =  [ 3 5 ] ;
boundaryConds(2).imposDispVals =  [ 0 0 ] ;
%md ant the external load is added into the same boundary condition using:
boundaryConds(2).loadsCoordSys = 'global'                   ;
boundaryConds(2).loadsTimeFact = @(t) p0*sin( omegaBar*t )  ;
boundaryConds(2).loadsBaseVals = [ 1 0 0 0 0 0 ]            ;
%md
%md#### Initial conditions
%md Initial displacement and velocity are set:
initialConds(1).nonHomogeneousUDofs    = 1   ;
initialConds(1).nonHomogeneousUVals    = u0  ;
initialConds(1).nonHomogeneousUdotDofs = 1   ;
initialConds(1).nonHomogeneousUdotVals = du0 ;
%md
%md#### Analysis settings
%md The following parameters correspond to the iterative trapezoidal Newmark method with the following tolerances, time step, tolerances and final time
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   0.005   ;
analysisSettings.finalTime     =   2.5*TN  ;
analysisSettings.stopTolDeltau =   1e-10    ;
analysisSettings.stopTolForces =   1e-10    ;
analysisSettings.stopTolIts    =   10      ;
%md
%md#### OtherParams
%md The nodalDispDamping is added into the model using:
otherParams.nodalDispDamping =   c    ;
%md The name of the problem is:
%md
otherParams.problemName = 'springMass_case1'     ;
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
% Only one element is considered with the first material and the second element setting
mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 2   ] ;
%md
%md Execute ONSAS and save the results:
[matUsNewmark, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md### Numerical case 2: nodal mass model with $\alpha$-HHT method and user loads function
%md
%md#### Material
%md The nodalMass field allows to add lumped matrices to a node, since this field is used, then the equivalent $\rho$ of the `material(1)` aforementioned now is set to 0. Although an equal mass $m$ is considered for $u_x$ $u_y$ and $u_z$ at the node $2$, so:
materials(1).density   = 0       ;
materials(2).nodalMass = [m m m] ;
%md
%md#### Boundary conditions
%md the boundary conditions struct is entirely re-written.
% repeat the BCs for node 1
boundaryConds = { } ;
boundaryConds(1).imposDispDofs =  [ 1 3 5 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;
% repeat the BCs for node 2
boundaryConds(2).imposDispDofs =  [ 3 5 ] ;
boundaryConds(2).imposDispVals =  [ 0 0 ] ;
%md ant the external load is added into the same boundary condition using:
boundaryConds(3).userLoadsFilename = 'myLoadSpringMass' ;
%md where inside the function 'myLoadSpringMass' the external force vector of the structure with 12 = (2x6) entries is computed.
%md
%md now the initial condition is added to the node $2$ with the second material:
mesh.conecCell{ 2, 1 } = [ 2 1 2 1   2  ] ;
%md
%md The $\alpha_{HHT}$ method with $\alpha=0$ is equivalent to Newmark, this is employed to validate results of both methods, as consequence is using:
analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.alphaHHT      =   0        ;
%md
otherParams.problemName = 'springMass_case2'     ;
%md
%md Execute ONSAS and save the results:
[matUsHHT, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md## Verification
%md---------------------
%md The numerical displacements of the node $2$ is extracted for both study cases:
valsNewmark = matUsNewmark(6+1,:) ;
valsHHT     = matUsHHT(6+1,:)     ;
%md The analytical solution is evaluated:
times       = linspace( 0,analysisSettings.finalTime, size(matUsHHT,2) ) ;
valsAnaly   = myAnalyticFunc(times)                                                          ;
%md The boolean to validate the implementation is evaluated such as:
analyticCheckTolerance = 5e-2 ;
verifBooleanNewmark =  ( ( norm( valsAnaly - valsNewmark ) / norm( valsAnaly ) ) <  analyticCheckTolerance ) ;
verifBooleanHHT     =  ( ( norm( valsAnaly - valsHHT     ) / norm( valsAnaly ) ) <  analyticCheckTolerance ) ;
verifBoolean        = verifBooleanHHT && verifBooleanNewmark                                                 ;
%md

%md---------------Bending beam mass-spring system--------
l = 10 ;
d = l/10;
A   = pi * d^2 /4; %0.1 ;
Izz = pi * d ^ 4 / 64 ; 
E = k*l*l*l/(3*Izz) ;   
rho = 2*m/(A*l) ; 
materials(1).hyperElasParams = [ E 0 ] ;
materials(1).density  = rho ;

materials(1).hyperElasModel  = 'linearElastic' ; % 1DrotEngStrain should work as well
%md
elements(1).elemType = 'node' ;
elements(2).elemType = 'frame'; %and not truss
elements(2).massMatType = 'lumped' ;
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = [d] ;
elements(2).elemTypeParams = 0 ;
%md
boundaryConds(1).imposDispDofs =  [ 1 2 3 4 5 6] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 0 0 0] ;
%md
boundaryConds(2).loadsCoordSys = 'global'                  ;
boundaryConds(2).loadsTieFact = @(t) p0*sin( omegaBar*t )        ;
%boundaryConds(2).loadsTimeFact = @(t) 100*p0;
boundaryConds(2).loadsBaseVals = [0 0 1 0 0 0 ] ; %along Y axis
%md An initial displacements $u_0$ is set in $y$ direction:
initialConds(1).nonHomogeneousUDofs   = 3;
initialConds(1).nonHomogeneousUVals = [u0];
initialConds(1).nonHomogeneousUdotDofs   = 3;
initialConds(1).nonHomogeneousUdotVals = [du0];
%md
mesh.nodesCoords = [  0  0  0 ; ...
                      l  0  0 ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 1   2   ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 2   ] ;
%md
otherParams.nodalDispDamping =   c    ;
otherParams.problemName = 'springMass_case3'     ;
%md
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   0.005   ;
analysisSettings.finalTime     =   2.5*TN  ;
%%
%---------------------------------------------------------
ku        = 39.47 ; % spring constant
cu        = 0   ; % damping parameter
m        = 1     ; % mass of the system
force = 0;
p0       = 1;
omegaBar = 4*sqrt(ku/m) ;
Ut0= 0.1; dUt0= 0; ddUt0 = 0;
% parameters for truss model
l   = 100 ;
d = l/20;
A   = pi * d^2 /4; %0.1 ;
%J = pi * d ^ 4 / 64 ; Iyy = J / 2 ; Izz = J / 2 ;
Izz = pi * d ^ 4 / 64 ; %actually
% Y deformation 
E = ku*l*l*l/(3*Izz) ;   % E   = ku * l /   A       ;
rho = 2*m/(A*l) ; % rho = m * 2 / ( A * l ) ;
omegaN = sqrt( ku / m );
freq   = omegaN / (2*pi)      ;
TN     = 2*pi / omegaN        
dtCrit = TN / pi              ;
dt = 0.028                ;
finalTime = 5*TN          ;

xi     = cu / m  / ( 2 * omegaN ) ;
nodalDamping = cu ;

% numerical method params
stopTolDeltau = 1e-10           ;
stopTolForces = 1e-10           ;
stopTolIts    = 30              ;
alphaHHT      = 0               ;
% ------------------------------------


materials(1).hyperElasModel  = '1DrotEngStrain' ;%'1DrotEngStrain' ;
materials(1).hyperElasParams = [ E 0 ] ;
materials(1).density  = rho ;

elements(1).elemType = 'node' ;
elements(2).elemType = 'frame';

elements(2).massMatType = 'lumped' ;
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = [d] ;%[sqrt(4*A/pi) ] ;
elements(2).elemTypeParams = 0 ;

boundaryConds(1).imposDispDofs =  [ 1 2 3 4 5 6] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 0 0 0] ;

%boundaryConds(2).imposDispDofs =  [ 1 5 ] ; % Y displacement
%boundaryConds(2).imposDispVals =  [ 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global'                  ;
% boundaryConds(2).loadsTimeFact = @(t) p0*sin( omegaBar*t )        ;
boundaryConds(2).loadsTimeFact = @(t) force        ;
boundaryConds(2).loadsBaseVals = [0 0 1 0 0 0 ] ;
%boundaryConds(2).userLoadsFilename = 'myLoadSpringMass' ;
% initial conditions
%initialConds.nonHomogeneousInitialCondU0   = [ 2 1 Ut0 ] ;
initialConds(1).nonHomogeneousUDofs   = 3;
initialConds(1).nonHomogeneousUVals = [Ut0];
%Not working
% initialConds(1).nonHomogeneousUdotDofs   = 3;
% initialConds(1).nonHomogeneousUdotVals = [dUt0];
mesh.nodesCoords = [  0  0  0 ; ...
                      l  0  0 ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 1   2   ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 2   ] ;

otherParams.nodalDispDamping = cu ;
analysisSettings.methodName    = 'newmark' ;
%md and the following parameters correspond to the iterative numerical analysis settings
analysisSettings.deltaT        =   dt  ;
analysisSettings.finalTime      =  finalTime  ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   10   ;
analysisSettings.alphaNM      =   0.25   ;
analysisSettings.deltaNM      =   0.5   ;
[matUsBending, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
valsBending = matUsBending(6+3,:) 
times = 0:dt:(finalTime+dt); 
figure(11)
hold on, grid on
plot(times, valsBending,'b-o')
xlabel('times')
ylabel('U')
title(sprintf('dt = %.3d, m = %d, c = %d, ku = %d, f = %d', dt, m, cu, ku, force))
legend('ONSAS mass-spring model', 'Newmark solution','location', 'South' ) ;
box off;
%%
%md------------------------------------------------------
%md## Plot verification
%md---------------------
%md The control displacement $u(t)$ is plotted:
figure
hold on, grid on, spanPlot = 8 ; lw = 2.0 ; ms = 11 ; plotfontsize = 20 ;
plot(times, valsAnaly   ,'b-', 'linewidth', lw,'markersize',ms )
plot(times(1:spanPlot:end), valsNewmark(1:spanPlot:end) ,'ro', 'linewidth', lw,'markersize', ms )
plot(times(1:spanPlot:end), valsHHT(1:spanPlot:end)     ,'gs', 'linewidth', lw,'markersize', ms )
labx = xlabel('t [s]');   laby = ylabel('u(t) [m]') ;
legend( 'analytic', 'truss-Newmark','nodalMass-HHT','Beam oscillator', 'location','southeast')
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
if exist('../../docs/src/assets/')==7
  % printing plot also to docs directory
  disp('printing plot also to docs directory')
  print('../../docs/src/assets/springMassCheckU.png','-dpng')
  else
  print('output/springMassCheckU.png','-dpng')
end
figure(10)
hold on, grid on,
plot(times(1:spanPlot:end), valsBending(1:spanPlot:end)     ,'y*', 'linewidth', lw,'markersize', ms )
%md
%md```@raw html
%md<img src="../../assets/springMassCheckU.png" alt="plot check" width="500"/>
%md```
%md
