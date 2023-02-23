%md# Simple Propeller example
%mdIn this example a simple propeller submitted to a constant uniform flow is considered. The geometry is given by three blades with circular cross section.
if ~(length(getenv('TESTS_RUN')) > 0 && strcmp( getenv('TESTS_RUN'), 'yes')), close all, clear all, end
%
addpath( genpath( [ pwd '/../../src'] ) );
%md## Problem definition
%mdThe blades are considered considerable stiff and only lift is considered, thus a rigid rotation analytic solution can be used to verify the numerical solution. 
%md
%mdThe wind velocity is given by the `windVel` function
va = feval('windVel', 0,0) ;
%md and the blades are considered to have only lift, with a lift coefficient given by the function `liftCoef`
c_l = feval('liftCoef', 0) ;
%md the density and kinematic viscosity are given by 
rhoA = 1.225 ; nuA = 1.6e-5 ;
%md
%mdThe material parameters of the blades is given by steel
E = 210e9 ;  nu = 0.3 ; rho = 6000 ; 
%mdand the geometric parameters length and diameter are
l = 3 ; d = 0.1;
%md
%md## Analytical solution
%mdcompute solution using the second cardinal, first the wind parameters are loaded
%md lift load per unit of length: 
fl = 1 / 2 * c_l * rhoA * norm(va) ^ 2 * d ;
%md the total moment induced in node 1 in x direction for is the sum for three blades: 
moment1x = 3 * fl * l * l / 2 ;
%mdand then the angular moment is:
bladeMass = rho * l * pi * d ^2 /4 ; 
Jrho =  3 * 1/3 * bladeMass  * l ^ 2 ; 
angleXnode1 = @(t)  moment1x / Jrho / 2 * t .^ 2 ;
%md
%md## Numerical solution
%md
%mdThe material parameters are
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]        ;
materials.density         = rho             ;
%md 
%md### elements
%mdThe elements are
% nodes
elements(1).elemType = 'node'  ;
% blades
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams = {'circle' ; d };
elements(2).massMatType =  'consistent'        ;
elements(2).elemTypeAero  = [0 d 0 4 0 ] ;
elements(2).aeroCoefs     = { []; 'liftCoef'; []  } ;
% auxiliar element
elements(3).elemType = 'truss' ;
elements(3).elemCrossSecParams = {'circle' ; 1.5*d };
elements(3).massMatType =  'lumped'        ;
%md 
%md### boundary Conditions
%md
boundaryConds(1).imposDispDofs = [ 1 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 ] ;
%md 
%md### mesh
%md
mesh.nodesCoords = [ 0        0              0            ; ...
                     0  l*sin( pi )        l*cos( pi )    ; ...
                     0  l*sin( pi/3  )     l*cos( pi/3 )  ; ... 
                     0  l*sin( 4*pi/3 )   -l*cos( 4*pi/3 ); ...
                     -d*.75 0 d ; ...
                     -d*.75 0 -l*1.5 ] ;
%md
mesh.conecCell         = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1   1   ] ;
mesh.conecCell{ 2, 1 } = [ 1 2 0   1 2 ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0   1 3 ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0   1 4 ] ;
% auxiliar elements
mesh.conecCell{ 5, 1 } = [ 0 1 1   5 ] ;
mesh.conecCell{ 6, 1 } = [ 0 1 1   6 ] ;
mesh.conecCell{ 7, 1 } = [ 1 3 0   5 6 ] ;
%md
%md### initial Conditions
%md
% homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
%md### analysisSettings
%md
analysisSettings.finalTime              =   400     ;
analysisSettings.deltaT                 =   5       ;
analysisSettings.methodName             = 'alphaHHT';
analysisSettings.stopTolIts             =   50      ;
analysisSettings.geometricNonLinearAero = true      ;
%analysisSettings.booleanSelfWeight      = false     ;
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
%md numerical rotation angle for nodal moment case is:
dofAngleXnode1 = 2 ;
angleXnode1Numeric = -matUs(dofAngleXnode1,:) ;
%md analytical rotation angle is:
angleXnode1Analytic = angleXnode1(timeVec) ;
%md
verifBoolean = norm( angleXnode1Numeric - angleXnode1Analytic )  ...
                    < ( norm( angleXnode1Analytic ) * 5e-2 ) ;
%md
%md## Plots
%md
lw = 2.0 ; ms = 10; plotfontsize = 22 ; spanPlotTime = 2 ;
fig1 = figure(1) ;
plot( timeVec(1:spanPlotTime:end), angleXnode1Analytic(1:spanPlotTime:end) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( timeVec(1:spanPlotTime:end), angleXnode1Numeric(1:spanPlotTime:end), 'ko' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time(s)');   laby = ylabel('$\theta_x node 1$') ;
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
%md```@raw html
%md<img src="https://github.com/ONSAS/ONSAS.m/blob/master/docs/src/assets/propeller.gif?raw=true" alt="propeller animation">
%md```