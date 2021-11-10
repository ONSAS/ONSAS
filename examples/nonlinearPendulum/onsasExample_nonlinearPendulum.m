
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
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
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
analysisSettings.deltaT        = 0.05  ;
analysisSettings.finalTime     = 1* T  ;
analysisSettings.stopTolDeltau = 1e-12 ;
analysisSettings.stopTolForces = 1e-12 ;
analysisSettings.stopTolIts    = 30    ;
otherParams.plotsFormat        = 'vtk' ;

%md### Analysis case 1: Solution using Newmark with truss element and mass lumped according to Bathe problem
analysisSettings.methodName = 'newmark' ;
analysisSettings.alphaNM    =   0.25    ;
analysisSettings.deltaNM    =   0.5     ;
otherParams.problemName     = 'nonlinearPendulumNewmarkTrussBathe';
% ------------------------------------
[matUsCase1, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

%md### Analysis case 2: Solution using HHT with truss element and mass lumped according to Bathe problem
otherParams.plotsFormat     = ''        ;
analysisSettings.methodName = 'alphaHHT';
analysisSettings.alphaHHT   =  0        ;        
otherParams.problemName     = 'nonlinearPendulumHHTTrussBathe';
% ------------------------------------
[matUsCase2, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

runFrameValidation = 1;
    if runFrameValidation
    %md### Analysis case 3: Solution using HHT with truss element and mass conssitent
    analysisSettings.finalTime  = 3.4124;
    analysisSettings.alphaHHT   = -0.05 ;        
    elements(2).elemTypeParams  = 1     ;
    otherParams.problemName     = 'nonlinearPendulumHHTTrussConssitent';
    % ------------------------------------
    [matUsCase3, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

    %md### Analysis case 4: Solution using HHT with truss element and mass conssitent
    elements(2).elemType        = 'frame';
    otherParams.problemName     = 'nonlinearPendulumHHTFrameConssitent';
    % ------------------------------------
    [matUsCase4, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
end

#md### extract control displacements of SVK solution
%mdMass displacement in z:
controlDofDispUz = 6 + 5 ;
controlDispZCase1 = matUsCase1( controlDofDispUz , : ) ;
controlDispZCase2 = matUsCase2( controlDofDispUz , : ) ;
if runFrameValidation
    controlDispZCase3 = matUsCase3( controlDofDispUz , : ) ;
    controlDispZCase4 = matUsCase4( controlDofDispUz , : ) ;
end
%mdMass displacement in x:
controlDofDispUx = 6 + 1 ;
controlDispXCase1 = matUsCase1( controlDofDispUx , : ) ;
controlDispXCase2 = matUsCase2( controlDofDispUx , : ) ;
if runFrameValidation
    controlDispXCase3 = matUsCase3( controlDofDispUx , : ) ;
    controlDispXCase4 = matUsCase4( controlDofDispUx , : ) ;
end
%mdAngle measured from the vertical:
angleThetaCase1 = rad2deg( atan2( ( l0 + controlDispXCase1 ), -controlDispZCase1 ) ) ;
angleThetaCase2 = rad2deg( atan2( ( l0 + controlDispXCase2 ), -controlDispZCase2 ) ) ;
if runFrameValidation
    angleThetaCase3 = rad2deg( atan2( ( l0 + controlDispXCase3 ), -controlDispZCase3 ) ) ;
    angleThetaCase4 = rad2deg( atan2( ( l0 + controlDispXCase4 ), -controlDispZCase4 ) ) ;
end

%mdCompute time of vectors
timesVec12  = (0:length(controlDispZCase1)-1) * analysisSettings.deltaT ;
timesVec34  = (0:length(controlDispZCase3)-1) * analysisSettings.deltaT ;


%md## verification
%md### verif newmark and HHT method using uz(N*T)= 0 with lumped mass matrix
tolVerifDisp = 1e-2 ;
verifBooleanCase1 =  ( abs( controlDispZCase1(end) / l0 ) <  tolVerifDisp ) ;
verifBooleanCase2 =  ( abs( controlDispZCase2(end) / l0 ) <  tolVerifDisp ) ;
%md### verif frame using HHT method vs truss
if runFrameValidation
    verifBooleanCase4 =  norm ( controlDispZCase3(end) - controlDispZCase4(end) ) / l0;
    verifBoolean    = verifBooleanCase1 && verifBooleanCase2 && verifBooleanCase4 ;
else
    verifBoolean    = verifBooleanCase1 && verifBooleanCase2 ;
end


%md### Plots
%md
%mdPlot parameters
MS = 10; LW = 1.5 ;
legendCase1 = [' Truss lumped Bathe Newmark \delta =' num2str(analysisSettings.deltaNM) ' and \alpha = ' num2str(analysisSettings.alphaNM)];
legendCase2 = [' Truss lumped Bathe HHT with \alpha =' num2str(analysisSettings.alphaHHT)];
if runFrameValidation
    legendCase3 = [' Truss conssitent HHT with \alpha =' num2str(analysisSettings.alphaHHT)];
    legendCase4 = [' Frame HHT with \alpha =' num2str(analysisSettings.alphaHHT)];
end
%mdPlot displacement solution of bathe problem
figure, hold on, grid on
plot( timesVec12, -controlDispZCase1, 'k-s' ,'markersize', MS,'linewidth', LW)
plot( timesVec12, -controlDispZCase2, 'bo','markersize', MS,'linewidth', LW)
xlabel('time (s)'), ylabel('mass displacement u_z (m)')
legend( legendCase1, legendCase2, 'location','NorthEast')
title("U_z solution bathe")
print('./output/dispPlot.png','-dpng')

%mdPlot a ngle solution
figure, hold on, grid on
plot( timesVec12, -angleThetaCase1, 'k-s' ,'markersize', MS,'linewidth', LW)
plot( timesVec12, -angleThetaCase2, 'bo','markersize', MS,'linewidth', LW)
xlabel('time (s)'), ylabel('\theta displacement (º)')
legend( legendCase1, legendCase2, 'location','NorthEast')
title("Θ solution bathe")
print('./output/thetaPlot.png','-dpng')

if runFrameValidation
    % mdPlot displacement solution of validation frame
    figure, hold on, grid on
    plot( timesVec34, -controlDispZCase3, 'c-s' ,'markersize', MS,'linewidth', LW)
    plot( timesVec34, -controlDispZCase4, 'ro','markersize', MS,'linewidth', LW)
    xlabel('time (s)'), ylabel('mass displacement u_z (m)')
    legend( legendCase3, legendCase4, 'location','NorthEast')
    title("U_z solution validation frame")
    %mdPlot a ngle solution
    figure, hold on, grid on
    plot( timesVec34, -angleThetaCase3, 'c-s' ,'markersize', MS,'linewidth', LW)
    plot( timesVec34, -angleThetaCase4, 'ro','markersize', MS,'linewidth', LW)
    xlabel('time (s)'), ylabel('\theta displacement (º)')
    legend( legendCase3, legendCase4, 'location','NorthEast')
    title("Θ solution frame")
end