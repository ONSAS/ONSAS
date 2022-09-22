%md# Deployable ring example
%md---
%md
%mdIn this tutorial...
%mdBefore defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar geometry and material parameters are defined.
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
% material scalar parameters
E  = 200e3 ;   nu = 0.3   ;
% geometrical scalar parameters
R = 120   ;  ty = .6 ;   tz = 6 ;
% the number of elements of the mesh
Nelem = 2*21 ;

materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;

otherParams.problemName = 'deployableRing' ;

elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;

elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [ty tz]     ;

boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;

boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsTimeFact = @(t) 1e2*t ;
boundaryConds(2).loadsBaseVals = [ 0 1 0 0 0 0 ] ;
boundaryConds(2).imposDispDofs = [ 3 4 5 6 ] ;
boundaryConds(2).imposDispVals = [ 0 0 0 0 ] ;

initialConds                = struct() ;

Np = Nelem +1;
angles = linspace(0, 2*pi, Np )' -pi*.5 ; angles(end) = [] ;

mesh.nodesCoords = [ R*(1+sin( angles )) -R*cos(angles) zeros(Np-1,1) ] ;

mesh.conecCell = {};
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  Nelem/2+1 ] ;
for i=1:(Nelem-1),
  mesh.conecCell{ end+1, 1} = [ 1 2 0 0  i i+1 ] ;
end
mesh.conecCell{ end+1, 1} = [ 1 2 0 0  Nelem 1 ] ; % last element of the circle

analysisSettings.methodName    = 'arcLength' ;
analysisSettings.incremArcLen = 1.5           ;
analysisSettings.posVariableLoadBC = 2 ;
analysisSettings.deltaT        =   1  ;
analysisSettings.finalTime     =   1400    ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   20   ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(analysisSettings.deltaT)/100 ;

otherParams.plots_format = 'vtk' ;
otherParams.plots_deltaTs_separation = 5 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


figure
plot(matUs((Nelem/2+1)*6-4,:), loadFactorsMat(:,2))