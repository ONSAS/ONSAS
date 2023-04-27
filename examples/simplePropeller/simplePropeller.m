%md# Simple Propeller example
%mdIn this example a simple propeller problem, inspired in an example [this article](https://arxiv.org/abs/2204.10545), is considered.
%md
%md## Problem definition
%md
%mdThe propeller has three blades with circular cross section and a uniform constant fluid flow is assumed as shown in the figure. Only lift is considered.
%md
%md```@raw html
%md<img src="../../assets/diagprop.png" alt="plot check" width="500"/>
%md```
%mdSince only lift is considered, an analytic solution can be computed. , first the wind parameters are loaded
%mdThe wind velocity is assumed constant and uniform, given by the function `windVel.m` in the same folder.
%md
%md## Analytical solution
%mdBefore computing the analytic solution the variables are cleared and the ONSAS/src path is added
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%mdThen the fluid velocity is computed
va = feval('windVel', 0,0) ;
%md
%mdThe density and kinematic viscosity are set: 
rhoA = 1.225 ; nuA = 1.6e-5 ;
%mdThe blades are considered to have only lift, with a lift coefficient given by the function `liftCoef.m`
c_l = feval('liftCoef', 0) ;
%mdThe material parameters of the blades correspond to steel
E = 210e9 ;  nu = 0.3 ; rho = 6000 ; 
%mdand the geometric parameters of the blades for length and diameter are set as
l = 3 ; d = 0.1;
%
%mdThen the lift load per unit of length is obtained as: 
fl = 1 / 2 * c_l * rhoA * norm(va)^2 * d ;
%md the total moment $M_x$ in node 1 is given by the sum of the moments for the three blades: 
moment1x = 3 * fl * l * l / 2 ;
%mdand then the angular moment is:
bladeMass = rho * l * pi * d ^2 /4 ; 
Jrho =  3 * 1/3 * bladeMass  * l ^ 2 ; 
%md Then, integrating, the angle $\theta_x$ can be obtained as a function of time.
angleXnode1 = @(t)  moment1x / Jrho / 2 * t .^ 2 ;
%md
%mdIf the the blades are considered considerably stiff and only lift is considered, this rigid-rotation analytic solution can be used to verify the numerical solution. 
%md
%md
%md## Numerical solution
%md
%mdThe material parameters are
materials                 = struct() ;
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]        ;
materials.density         = rho             ;
%md 
%md### elements
%mdThe elements are given by: nodes
elements             = struct() ;
elements(1).elemType = 'node'  ;
%mdframe elements for modelling the blades
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams = {'circle' ; d };
elements(2).massMatType =  'consistent'        ;
%md with the definition of the aerodynamic forces
elements(2).dragCoefFunction = @(beta,Re) 0;
elements(2).liftCoefFunction ='liftCoef';
elements(2).momentCoefFunction = @(beta,Re) 0;
%md and a auxiliar truss element for the pole
elements(3).elemType = 'truss' ;
elements(3).elemCrossSecParams = {'circle' ; 1.5*d };
elements(3).massMatType =  'lumped'        ;
%md 
%md### boundary Conditions
%md The only boundary condition is the one of the center node, with all the dofs fixed except for the rotation $\theta_x$
boundaryConditions             = struct() ;
boundaryConds(1).imposDispDofs = [ 1 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 ] ;
%md 
%md### mesh
%md The mesh is defined
mesh             = struct() ;
mesh.nodesCoords = [ 0        0              0            ; ...
                     0  l*sin( pi )        l*cos( pi )    ; ...
                     0  l*sin( pi/3  )     l*cos( pi/3 )  ; ... 
                     0  l*sin( 4*pi/3 )   -l*cos( 4*pi/3 ); ...
                     -d*.75 0 d ; ...
                     -d*.75 0 -l*1.5 ] ;
%md And a simple connectivity is required since only one frame element is used for each blade
mesh.conecCell         = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1   1   ] ;
mesh.conecCell{ 2, 1 } = [ 1 2 0   1 2 ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0   1 3 ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0   1 4 ] ;
%md These are auxiliar elements (nodes and truss) used to model the pole
mesh.conecCell{ 5, 1 } = [ 0 1 1   5 ] ;
mesh.conecCell{ 6, 1 } = [ 0 1 1   6 ] ;
mesh.conecCell{ 7, 1 } = [ 1 3 0   5 6 ] ;
%md
%md### initial Conditions
%md
%mdhomogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
%md### analysisSettings
%md The analysis settings are set
analysisSettings                        = struct() ;
analysisSettings.finalTime              =   400     ;
analysisSettings.deltaT                 =   5       ;
analysisSettings.methodName             = 'alphaHHT';
analysisSettings.stopTolIts             =   50      ;
analysisSettings.geometricNonLinearAero = true      ;
analysisSettings.stopTolDeltau          =   0       ;
analysisSettings.stopTolForces          =   1e-5    ;
analysisSettings.fluidProps = { rhoA ; nuA ; 'windVel' } ;
%md
%md### otherParams
%md
otherParams.problemName =  'simplePropeller' ;
otherParams.plots_format = 'vtk' ;
%md
%md### Run ONSAS
[ matUs, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
%md
%md## Verification
%md
%md numerical time vector is given by:
timeVec = linspace(0, analysisSettings.finalTime, size(matUs, 2) ) ;
%md numerical rotation is obtained from the ONSAS matUs:
dofAngleXnode1 = 2 ;
angleXnode1Numeric = -matUs(dofAngleXnode1,:) ;
%mdThe analytic rotation is:
angleXnode1Analytic = angleXnode1(timeVec) ;
%md and the norm of the difference is computed and the test is verified
verifBoolean = norm( angleXnode1Numeric - angleXnode1Analytic )  ...
                    < ( norm( angleXnode1Analytic ) * 5e-2 ) ;
%md
%md## Plots
%md
lw = 2.0 ; ms = 10; plotfontsize = 22 ; spanPlotTime = 2 ;
fig1 = figure ;
plot( timeVec(1:spanPlotTime:end), angleXnode1Analytic(1:spanPlotTime:end) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( timeVec(1:spanPlotTime:end), angleXnode1Numeric(1:spanPlotTime:end), 'ko' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time(s)');   laby = ylabel('\theta_x node 1') ;
legend('analytic','numeric', 'location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('simple propeller test')
if length(getenv('TESTS_RUN')) > 0 && strcmp( getenv('TESTS_RUN'), 'yes')
  fprintf('\ngenerating output png for docs.\n')
  print(fig1, 'output/verifPropeller.png','-dpng')
else
  fprintf('\n === NOT in docs workflow. ===\n')
end
%md
%md```@raw html
%md<img src="../../assets/generated/verifPropeller.png" alt="plot check" width="500"/>
%md```
%md
%mdThe obtained simulation is: 
%md```@raw html
%md<img src="https://github.com/ONSAS/ONSAS.m/blob/master/docs/src/assets/propeller.gif?raw=true" alt="propeller animation">
%md```