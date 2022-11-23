%md# Cantilever Modal Analysis
%md
%mdBefore defining the structs, the workspace is cleaned and the ONSAS directory is added to the path
close all, clear all, addpath( genpath( [ pwd '/../../src'] ) );
%md
%mdThe material scalar parameters are set.
E = 200e9 ; nu = 0.3;  rho = 700;
%mdThe cross-section of the beam is rectangular. The widths and other geometry scalar parameters are computed.
l = 10 ; diam = .01 ; 
numElements = 10 ; % Number of elements

tf     = 8     ; % s
deltat = 0.1   ; % s

global exportFirstMatrices

exportFirstMatrices = true

%md### materials
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ]         ;
materials.density         = rho              ;
%md
%md### elements
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%elements(2).elemCrossSecParams = { 'rectangle' , [ty tz] } ;
elements(2).elemCrossSecParams = { 'circle' , [diam] } ;
%md The consistent mass approach is considered for the dynamic analysis
elements(2).massMatType = 'consistent';
%md
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%md and the second corresponds to a time dependant external force
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsTimeFact = @(t) t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ] ;
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;

mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  numElements+1 ] ;

for i=1:numElements
  mesh.conecCell{ i+2, 1 } = [ 1 2 0 0  i i+1 ] ;
end

%md### analysisSettings
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   deltat  ;
analysisSettings.finalTime     =   tf   ;
analysisSettings.stopTolDeltau =   1e-10 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   20   ;
%md
%md## otherParams
otherParams.problemName = 'cantilever_modal_analysis';
otherParams.plots_format = 'vtk' ;
%md ONSAS execution
[coRotMatUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md


addpath('./output');
filename = './output/matrices.mat';
load(filename);

KTred = KT( neumdofs, neumdofs );

eigs(KTred,10,'sa')

stop

Mred = massMat( neumdofs, neumdofs );

Mred = Mred + speye(size(Mred)(1));

%[PHI, OMEGA] = eigs(Mred^(-1)*Kred,10,'sm');

