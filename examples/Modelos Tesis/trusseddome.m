%md# Plasticity | shallow trussed dome
close all; clear;
addpath( genpath( [ pwd '/../../src'] ) );

% scalar parameters (N/cm2)
E = 200e3 ;
Kplas = 200 ;
sigma_Y_0 = 250 ;
Fu = 1 ;

% x and z coordinates of node 4
x2 = 300 ;
z2 = -400 ;

materials = struct();
materials.hyperElasModel  = 'isotropicHardening' ;

% /\ isotropic hardening include use of logarithmic strain /\

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
analysisSettings.finalTime     =   120  ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

analysisSettings.posVariableLoadBC = 2 ;

otherParams = struct();
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ;

otherParams.problemName       = 'Newton-Raphson_Arc-Length_Logarithmic_Strain_Jirasek' ;
analysisSettings.methodName   = 'arcLength' ;
analysisSettings.finalTime    = 120 ;
analysisSettings.incremArcLen = [0.8/20*ones(1,20) -3*0.8/40*ones(1,40) 2*0.8/60*ones(1,60)] ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.1)/120 ;
analysisSettings.posVariableLoadBC = 2 ;

global arcLengthFlag
arcLengthFlag = 2 ;

global dominantDofs
dominantDofs = 6*3+5 ;

global scalingProjection
scalingProjection = -1 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNRAL_Jirasek_logarithmic_strain =  -matUs(6*3+5,:) ;
loadFactorsNRAL_Jirasek_logarithmic_strain  =  loadFactorsMat(:,2) ;

figure
plot( controlDispsNRAL_Jirasek_logarithmic_strain, loadFactorsNRAL_Jirasek_logarithmic_strain, 'linewidth', 1.5)
labx = xlabel('Displacement w(t)');
laby = ylabel('\lambda(t)') ;
legend('\fontsize{12} NRAL-Jirasek-Logarithmic Strain','Location','southeast');