% EulerColumn
close all, clear all
% add path
addpath( genpath( [ pwd '/../../src'] ) );
% material scalar parameters
E = 30e6 ; % kN/m2
nu = 0.2 ;
% geometrical scalar parameters
l = 5 ; %m
ty = .1 ; %m
tz = .1 ; %m
% the number of elements of the mesh
numElements = 8 ;

% MEBI parameters
% Materials
% ----------------------------------------------------------------------
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
% Elements
% ----------------------------------------------------------------------
% Types
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
% Sections
elements(2).elemCrossSecParams = { 'rectangle', [ ty tz ] } ;
elements(2).massMatType   = 'consistent'       ;

% Boundary conditions
% ----------------------------------------------------------------------
% Pinned support
boundaryConds(1).imposDispDofs = [ 1 3 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
% Roller support
boundaryConds(2).imposDispDofs = [ 1 3 6 ] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;
% Load
P = -1 ;
imp = P/1e4 ;
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 imp P 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;

% Initial conditions
% ----------------------------------------------------------------------
initialConds                = struct() ;


% Mesh
% Nodes coords
% ----------------------------------------------------------------------
mesh.nodesCoords = [ zeros(numElements+1,2) (0:(numElements))'*l/numElements ] ;
% Conec cell
% ----------------------------------------------------------------------
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0 1 ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0 numElements+1 ] ;
for i=1:numElements
  mesh.conecCell{ i+2,1 } = [ 1 2 0 0  i i+1 ] ;
end
% Analysis settings

% Parameters
% ----------------------------------------------------------------------
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   2.5    ;
analysisSettings.finalTime     =   125    ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   20   ;

otherParams.problemName = 'EulerColumn';
otherParams.controlDofs = [ numElements+1  5 ] ;
%otherParams.plots_format = 'vtk' ;

A = ty*tz ;
I = max(ty,tz)*min(ty,tz)^3/12 ;
Pcrit = pi()^2*E*I/l^2 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


Dof      		 = (numElements/2 + 1)*6 - 5 	;
controlDisps =  matUs(Dof, :) 						;
loadFactors  =  loadFactorsMat(:, 2) 			;

lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;
figure
grid on
plot( controlDisps, loadFactors, 'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;


verifBoolean = ( abs(controlDisps(end) - 1.7406)/1.7406 ) < 1e-3
