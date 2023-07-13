% Thermal
close all, clear all
% add path

% uncomment your path to ONSAS
addpath( genpath( [ pwd '/../src'] ) );
% or use an environment variable
%addpath( genpath( [ getenv('ONSAS_PATH') ] ) );

% material scalar parameters
E = 30e6 ; % kN/m2
nu = 0.2 ;
% geometrical scalar parameters
l = .1 ; %m
ty = .02  ; %m
tz = .002 ; %m
% the number of elements of the mesh
numElements = 8 ;

% MEBI parameters
% Materials
% ----------------------------------------------------------------------
materials  = struct();
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ E nu ] ;

materials(2).hyperElasModel  = '1DrotEngStrain' ;
materials(2).hyperElasParams = [ E nu ] ;
materials(2).thermalExpansion = 1e-5 ;
% Elements
% ----------------------------------------------------------------------
% Types
elements  = struct();
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
% Sections
elements(2).elemCrossSecParams = { 'rectangle', [ ty tz ] } ;
elements(2).massMatType   = 'consistent'       ;
elements(3).elemType = 'truss' ;
elements(3).elemCrossSecParams = { 'rectangle', [ ty tz ] } ;

% Boundary conditions
% ----------------------------------------------------------------------
% Pinned support
boundaryConds  = struct();
boundaryConds(1).imposDispDofs = [ 1 2 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
% Roller support
boundaryConds(2).imposDispDofs = [ 3 5 ] ;
boundaryConds(2).imposDispVals = [ 0 0 ] ;

P = -1 ;
imp = P/1e5 ;
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 imp 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;

global temperature
temperature = @(t) -100*t + 2*100*(t-10)*(t>10) ;
% Initial conditions
% ----------------------------------------------------------------------
initialConds                = struct() ;


% Mesh
% Nodes coords
% ----------------------------------------------------------------------
mesh  = struct();
mesh.nodesCoords = [ (0:(numElements))'*l/numElements zeros(numElements+1,2) ] ;
% Conec cell
% ----------------------------------------------------------------------
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 1 ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2  numElements+1 ] ;
for i=1:numElements
  mesh.conecCell{ i+2,1 } = [ 1 2 0  i i+1 ] ;
end
mesh.conecCell{ end+1,1 } = [ 2 3 0  1 numElements+1 ] ;
% Analysis settings

% Parameters
% ----------------------------------------------------------------------
analysisSettings  = struct();
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   .10    ;
analysisSettings.finalTime     =   20    ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   20   ;
%md
otherParams  = struct();
otherParams.problemName = 'thermal';
otherParams.controlDofs = [ numElements+1  5 ] ;
otherParams.plots_format = 'vtk' ;

A = ty*tz ;
I = max(ty,tz)*min(ty,tz)^3/12 ;
Pcrit = pi()^2*E*I/l^2 ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


Dof      		 = (numElements/2 + 1)*6 - 1 	;
controlDisps =  matUs(Dof, :) 						;

lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;
figure
grid on
plot( controlDisps, 'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;


verifBoolean = ( abs(controlDisps(end) - 1.7406)/1.7406 ) < 1e-3 ;
