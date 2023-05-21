%md# Plasticity | classical truss, 3 rods
close all; clear;
addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters (N/cm2)
E = 210e3 ;
Kplas = 1093 ;
sigma_Y_0 = 330.3 ;
Fu = 1 ;

% x and z coordinates of node 4
x2 = 300 ;
z2 = -400 ;

materials = struct();
materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements = struct();
elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(3).elemType = 'truss';
elements(4).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(1*4/pi)} ;
elements(3).elemCrossSecParams = { 'circle' , sqrt(1*4/pi)} ;
elements(4).elemCrossSecParams = { 'circle' , sqrt(1*4/pi)} ;

boundaryConds = struct();
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs = 3 ;
boundaryConds(2).imposDispVals = 0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (t<=1)*(Fu)*t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -1 0 ]   ;

mesh = struct();
mesh.nodesCoords = [   0  0   0     ; ...
                      x2  0   0     ; ...
                    2*x2  0   0     ; ...
                      x2  0   z2 ]  ;

% MEBI [Material Element Boundary_Conditions Initial_Conditions]

mesh.conecCell = cell(7,1) ;
mesh.conecCell{ 1, 1 } = [ 0 1 1  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1  2   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 1  3   ] ;
mesh.conecCell{ 4, 1 } = [ 0 1 2  4   ] ;
mesh.conecCell{ 5, 1 } = [ 1 2 0  1 4 ] ;
mesh.conecCell{ 6, 1 } = [ 1 3 0  2 4 ] ;
mesh.conecCell{ 7, 1 } = [ 1 4 0  3 4 ] ;

initialConds                = struct() ;

analysisSettings = struct();
analysisSettings.deltaT        =   1    ;
analysisSettings.finalTime     =   100  ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

analysisSettings.posVariableLoadBC = 2 ;

otherParams = struct();
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

% an Eternal Golden Braid

otherParams.problemName       = 'NRAL_Jirasek' ;
analysisSettings.methodName   = 'arcLength' ;
analysisSettings.finalTime    = 100 ;
analysisSettings.incremArcLen = [0.2*ones(1,33) 0.99*ones(1,33) -0.65*ones(1,34)] ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.1)/100 ;
analysisSettings.posVariableLoadBC = 2 ;

global arcLengthFlag
arcLengthFlag = 2 ;

global dominantDofs
dominantDofs = 6*3+5 ;

global scalingProjection
scalingProjection = -1 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNRAL_Jirasek_Green =  -matUs(6*3+5,:) ;
loadFactorsNRAL_Jirasek_Green  =  loadFactorsMat(:,2) ;

figure
plot( controlDispsNRAL_Jirasek_Green, loadFactorsNRAL_Jirasek_Green, 'linewidth', 1.5)
labx = xlabel('Displacement w(t)');
laby = ylabel('\lambda(t)') ;
legend('\fontsize{12} NRAL-Jirasek-Green','Location','southeast');