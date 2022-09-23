%md# Simple wind turbine example
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) );
%md
%md# General  problem parameters
%md
%md material scalar parameters are
%md```math
E = 210e9 ;  nu = 0.3 ; rho = 6000 ; G = E / (2 * (1+nu)) ;
%md```
% geometrical scalar parameters
%md```math
l = 3 ; d = 0.1; 
%md```
%md# analytical solution
%mdcompute solution by the second carindal, first the wind parameters are loaded
%md```math
rhoA = 1.225 ; nuA = 1.6e-5;   c_l = feval('liftCoef', 0) ; vwind = feval('windVel', 0,0) ;
%md```
%md lift load per unit of length: 
%md```math
fl = 1 / 2 * c_l * rhoA * norm(vwind) ^ 2 * d ;
%md```math
%md the total moment induced in node 1 in x direction for is the sum for three blades: 
%md```math
moment1x = 3 * fl * l * l / 2 ;
%md```
%mdand then the angular moment is:
%md```math
bladeMass = rho * l * pi * d ^2 /4 ; 
Jrho =  3 * 1/3 * bladeMass  * l ^ 2 ; 
angleXnode1 = @(t)  moment1x / Jrho / 2 * t .^ 2 ;
%md```
%md ## Nodal moment case
%mdThe first case proposed aims to validate inertial force undergoing large rotations on a simple wind turbine example.The moment applied is equal to aerodynamic moment induced by external aerodynamic force
%md 
%md### materials
%md
% Since the example contains only aeroFoone rod the fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]        ;
materials.density         = rho             ;
%md 
%md### elements
%md
%md nodes
elements(1).elemType = 'node'  ;
%md beams
elements(2).elemType = 'frame' ;
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} =  d       ;
elements(2).massMatType =  'consistent'        ;
%md 
%md### boundary Conditions
%md
%mdThe elements are submitted to two different BC settings. The first BC corresponds to a free angle in x condition 
boundaryConds(1).imposDispDofs = [ 1 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 ] ;
%mdthen the nodal moment load equivalent is:
boundaryConds(1).loadsCoordSys = 'global' ;
boundaryConds(1).loadsTimeFact = @(t) moment1x ;
boundaryConds(1).loadsBaseVals = [ 0 -1 0 0 0 0 ] ;
%md 
%md### initial Conditions
%md
% homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
%md### mesh parameters
%md
mesh.nodesCoords = [ 0        0              0            ; ...
                     0  l*sin( pi )        l*cos( pi )    ; ...
                     0  l*sin( pi/3  )     l*cos( pi/3 )  ; ... 
                     0  l*sin( 4*pi/3 )   -l*cos( 4*pi/3 ); ] ;
%md
mesh.conecCell         = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
mesh.conecCell{ 2, 1 } = [ 1 2 0 0   1 2 ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 3 ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 0   1 4 ] ;
%md
%md## analysisSettings
%md
analysisSettings.finalTime              =   450     ;
analysisSettings.deltaT                 =   5       ;
analysisSettings.methodName             = 'alphaHHT';
analysisSettings.stopTolIts             =   50      ;
analysisSettings.geometricNonLinearAero = true      ;
analysisSettings.booleanSelfWeight      = false     ;
analysisSettings.stopTolDeltau          =   0       ;
analysisSettings.stopTolForces          =   1e-5    ;
%md
%md### otherParams
%md
otherParams.problemName = strcat( 'onsasExample_simpleWindTurbine_nodalMoment' ) ;
otherParams.plotsFormat = '' ;
%
%md# Execute ONSAS
[ matUs_nodalMoment, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
%md
%md ## Aeroydinaic force case
%md
%md### mesh parameters
%md
mesh.conecCell{ 3, 1 } = [ 1 3 0 0   1 3 ] ;
mesh.conecCell{ 4, 1 } = [ 1 4 0 0   1 4 ] ;
%md
%md## elements
%md
%mdAdd aerodynamic properties into elements struct:
numGaussPoints = 4 ;
computeAeroTangentMatrix = false ;
elements(2).elemTypeAero  = [0 0 d numGaussPoints computeAeroTangentMatrix ] ;
elements(2).aeroCoefs     = { []; 'liftCoef'; []  } ;
%md second blade in (z,-y) quarter 
elements(3).elemType                = 'frame'                                          ;
elements(3).elemCrossSecParams{1,1} = 'circle'                                         ;
elements(3).elemCrossSecParams{2,1} =  d                                               ;
elements(3).elemTypeAero            = [0 d 0 numGaussPoints computeAeroTangentMatrix ] ;
elements(3).aeroCoefs               = { []; 'liftCoef'; []  }                          ;
elements(3).massMatType             =  'consistent'                                    ;
%md third blade in (z,y) quarter 
elements(4).elemType                = 'frame'                                           ;
elements(4).elemCrossSecParams{1,1} = 'circle'                                          ;
elements(4).elemCrossSecParams{2,1} =  d                                                ;
elements(4).elemTypeAero            = [0 -d 0 numGaussPoints computeAeroTangentMatrix ] ;
elements(4).aeroCoefs               = { [], 'liftCoef', []  }                           ;
elements(4).massMatType             =  'consistent'                                     ;
%md
%md ## boundary Conditions
%md
%md then delete boundary conditions into the struct:
boundaryConds(1).loadsCoordSys = [] ;
boundaryConds(1).loadsTimeFact = [] ;
boundaryConds(1).loadsBaseVals = [] ;
%md
%md ## analysisSettings
%md
%md add wind velocity function into analysis settings struct:
analysisSettings.fluidProps = { rhoA ; nuA ; 'windVel' } ;
%md### otherParams
%md
otherParams.problemName = strcat( 'onsasExample_simpleWindTurbine_aeroForce' ) ;
otherParams.plotsFormat = 'vtk' ;
%md# Execute ONSAS
[ matUs_aeroForce, ~ ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ; 
%md
%md## Verification
%md
%md numerical time vector is given by:
timeVec = linspace(0, analysisSettings.finalTime, size(matUs_nodalMoment, 2) ) ;
%md numerical rotation angle for nodal moment case is:
dofAngleXnode1 = 2 ;
angleXnode1Numeric_NodalMoment = -matUs_nodalMoment(dofAngleXnode1,:) ;
%md numerical rotation angle for nodal moment case is:
angleXnode1Numeric_aeroForce = -matUs_aeroForce(dofAngleXnode1,:) ;
%md analytical rotation angle is:
angleXnode1Analytic = angleXnode1(timeVec) ;
%md
verifBoolean = norm( angleXnode1Numeric_aeroForce - angleXnode1Analytic )  ...
                    < ( norm( angleXnode1Numeric_aeroForce ) * 1e-1 ) ;
%md
%md## Plots
%md
lw = 2.0 ; ms = 10; plotfontsize = 22 ;
spanPlotTime = 2 ;
fig1 = figure(1) ;
plot( timeVec(1:spanPlotTime:end), angleXnode1Analytic(1:spanPlotTime:end) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( timeVec(1:spanPlotTime:end), angleXnode1Numeric_NodalMoment(1:spanPlotTime:end), 'ko' , 'linewidth', lw,'markersize',ms )
plot( timeVec(1:spanPlotTime:end), angleXnode1Numeric_aeroForce(1:spanPlotTime:end), 'rs' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time(s)');   laby = ylabel('$\theta_x node 1$') ;
legend('analytic','nodal moment', 'aero force', 'location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print(fig1, 'output/verifSimpleWindTurbine.png','-dpng')
close(1)