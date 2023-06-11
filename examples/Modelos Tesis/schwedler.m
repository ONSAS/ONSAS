%md# Plasticity | Schwedler Dome (97 nodes, 264 finite elements)
close all; clear;
addpath( genpath( [ pwd '/../../src'] ) );

% scalar parameters (N,m)
E = 2.0e11 ;
Kplas = 0 ;
sigma_Y_0 = 25000e3 ;
Fu = 1 ;

materials = struct();
materials.hyperElasModel  = 'isotropicHardening' ;

% /\ isotropic hardening include use of logarithmic strain /\

materials.hyperElasParams = [ E Kplas sigma_Y_0 ] ;

elements = struct();
elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(0.0032*4/pi)} ;

boundaryConds = struct();
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) (Fu)*t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ] ;

% MEBI [Material Element Boundary_Conditions Initial_Conditions]

base_msh='' ;

mesh = struct() ;
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( [ base_msh 'Schwedler.msh'] ) ;

initialConds                = struct() ;

analysisSettings = struct();

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;

analysisSettings.posVariableLoadBC = 2 ;

otherParams = struct();
otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 2 ; 

otherParams.problemName       = 'Schwedler Dome' ;
analysisSettings.methodName   = 'arcLength' ;
analysisSettings.finalTime    = 300 ;
analysisSettings.incremArcLen = [2/100*ones(1,100) 2/100*ones(1,100) 0.58/100*ones(1,100)] ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(1)/300 ;
analysisSettings.posVariableLoadBC = 2 ;

global arcLengthFlag
arcLengthFlag = 2 ;

global dominantDofs
dominantDofs = 6*96+5 ;

global scalingProjection
scalingProjection = -1 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNRAL_Jirasek_logarithmic_strain =  -matUs(6*96+5,:) ;
loadFactorsNRAL_Jirasek_logarithmic_strain  =  -loadFactorsMat(:,2)/1000 ;

figure(1)
plot( controlDispsNRAL_Jirasek_logarithmic_strain, loadFactorsNRAL_Jirasek_logarithmic_strain, 'linewidth', 1.5)
labx = xlabel('Displacement w(t)');
laby = ylabel('\lambda(t)') ;
legend('\fontsize{12} NRAL-Jirasek-Logarithmic Strain','Location','southeast');
hold  on;