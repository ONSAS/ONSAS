
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
%md## Numerical solution
%mdExample of nolinear equaitons in dynamic analyisis section of Bathe pag 827
%md
%md##Parameters
EA = 1e8 ;
A = 0.1  ;
E   = EA / A ; nu = 0  ;
l0  = 3.0443 ;  T = 4.13 ;
m   = 10 ;      g = 9.80 ;
rho = 2*m / ( A * l0 ) ;
%md
%md#### materials
materials.hyperElasModel  = 'SVK' ;
lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials.hyperElasParams = [ lambda mu ] ;
materials.density = rho ;
%md
%md#### elements
elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemTypeGeometry = [2 sqrt(A) sqrt(A) ] ;
elements(2).elemTypeParams = 0 ;
%md
%md#### boundaryConds
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs =  3 ;
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) 1.0 ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -m*g 0 ] ;
%md
%md#### initialsConds
initialConds                = struct() ;
%md
%md### mesh set
mesh.nodesCoords = [   0  0   l0 ; ...
                      l0  0  l0  ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  2   ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0  1 2 ] ;
%md
%md### common analysis settings
analysisSettings.deltaT        = 0.025 ;
analysisSettings.finalTime     = 3 *T  ;
analysisSettings.stopTolDeltau = 1e-10 ;
analysisSettings.stopTolForces = 1e-10 ;
analysisSettings.stopTolIts    = 30    ;
otherParams.plotsFormat = 'vtk' ;

%md### Analysis case 1: Solution using Newmark
analysisSettings.methodName    = 'newmark' ;
analysisSettings.alphaNM      =   0.25   ;
analysisSettings.deltaNM      =   0.5   ;
otherParams.problemName = 'nonlinearPendulumNewmark';
% ------------------------------------
[matUsNew, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

%md### Analysis case 2: Solution using HHT
analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.alphaHHT      =   -0.05   ;
otherParams.problemName = 'nonlinearPendulumHHT';
% ------------------------------------
[matUsHHT, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

#md## verification

%mdMass displacement in z:
controlDofDispUz = 6 + 5 ;
controlDispsMassNewmarkZ = matUsNew( controlDofDispUz , : ) ;
controlDispsMassHHTZ     = matUsHHT( controlDofDispUz , : ) ;

%mdMass displacement in x:
controlDofDispUx = 6 + 1 ;
controlDispsMassNewmarkX = matUsNew( controlDofDispUx , : ) ;
controlDispsMassHHTX     = matUsHHT( controlDofDispUx , : ) ;

%mdAngle measured the vertical:
angleThetaNew = rad2deg( atan2( ( l0 + controlDispsMassNewmarkX ), -controlDispsMassNewmarkZ ) ) ;
angleThetaHHT = rad2deg( atan2( ( l0 + controlDispsMassHHTX ), -controlDispsMassHHTZ ) )         ;       

timesVec             = (0:length(controlDispsMassNewmarkZ)-1) * analysisSettings.deltaT ;


tolVerifDisp = 1e-2 ;
verifBooleanNew =  ( ( abs( controlDispsMassNewmarkZ(end)  ) / abs( l0 ) ) <  tolVerifDisp ) ;
verifBooleanHHT =  ( ( abs( controlDispsMassHHTZ    (end)  ) / abs( l0 ) ) <  tolVerifDisp ) ;
verifBoolean    = verifBooleanNew && verifBooleanHHT ;


%md### Plots
%md
%mdPlot parameters
MS = 10; LW = 1.5 ;
legendNew = [' Newmark with \delta =' num2str(analysisSettings.deltaNM) ' and \alpha = ' num2str(analysisSettings.alphaNM)];
legendHHT = [' HHT with \alpha =' num2str(analysisSettings.alphaHHT) ];

figure, hold on, grid on
plot( timesVec, -controlDispsMassNewmarkZ, 'b-o','markersize',MS,'linewidth',LW)
plot( timesVec, -controlDispsMassNewmarkZ, 'r-x','markersize',MS,'linewidth',LW)
xlabel('time (s)'), ylabel('mass displacement u_z (m)')
legend( legendNew, legendHHT, 'location','NorthEast')


figure, hold on, grid on
plot( timesVec, angleThetaNew, 'k-o','markersize',MS,'linewidth',LW)
plot( timesVec, angleThetaHHT, 'c-x','markersize',MS,'linewidth',LW)
xlabel('time (s)'), ylabel('mass angle \theta (º)')
legend( legendNew, legendHHT, 'location','NorthEast')

print('./output/output.png','-dpng')
