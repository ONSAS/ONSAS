% Orthogonal grid

close all, clear all, close all
% add path
addpath( genpath( [ pwd '/../../../../repos/ONSAS.m/src'] ) );

% material scalar parameters
E = 30e6 ; % kN/m2
nu = 0.2 ;
% mesh
l = 7 ; %m
width = 1 ;
thk = 0.5 ;
nrows = 2 ;
ncolumns = 2 ;

ty = width ;
tz = thk ; 



% MEBI parameters

% Materials
% ----------------------------------------------------------------------
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
% Elements
% ----------------------------------------------------------------------
% Types
elements(1).elemType = 'node'  ;

% Columns
elements(2).elemType = 'frame' ;

% rows
elements(3).elemType = 'frame' ;

% Sections
% Columns
elements(2).elemTypeGeometry = [2 ty tz ] ;
elements(2).elemTypeParams   = 1          ;

% Rows
elements(3).elemTypeGeometry = [2 ty/2 tz/2 ] ;
elements(3).elemTypeParams   = 1          ;


% Boundary conditions
% ----------------------------------------------------------------------
% Pinned support
boundaryConds(1).imposDispDofs = [ 1 3 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
% Roller support
boundaryConds(2).imposDispDofs = [ 1 3 5 6 ] ;
boundaryConds(2).imposDispVals = [ 0 0 0 0 ] ;
% Load
P = -1 ;
imp = P/100 ;

boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsBaseVals = [ 0 imp 0 0 P 0 ] ;

% Initial conditions
% ----------------------------------------------------------------------
initialConds                = struct() ;

% Mesh

% Nodes coords
% ----------------------------------------------------------------------
mesh.nodesCoords = [ ] ;


mesh.nodesCoords = [ 0 0 0 ; 0 1 0 ; 0 2 0 ; 0 1 1 ; 0 1 -1 ] ; 

% Conec cell
% ----------------------------------------------------------------------
mesh.conecCell = {  } ;

% First BC
mesh.conecCell{ 1, 1 } = [ 0 1 1 0 4 ] ;

% Second BC
mesh.conecCell{ 2, 1 } = [ 0 1 2 0 5 ] ;

% Column elements
mesh.conecCell{ 3, 1 } = [ 1 2 0 0 4 2 ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 0 2 5 ] ;

% Row elements
mesh.conecCell{ 5, 1 } = [ 1 2 0 0 1 2 ] ;
mesh.conecCell{ 6, 1 } = [ 1 2 0 0 2 3 ] ;



% Analysis settings

% Parameters
% ----------------------------------------------------------------------
analysisSettings.methodName    = 'newtonRaphson' ;
boundaryConds(2).loadsTimeFact = @(t) 10*t     ;

analysisSettings.deltaT        =   .1  ;
analysisSettings.finalTime     =   0.1   ;
analysisSettings.incremArcLen  =   0.005   ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;
analysisSettings.posVariableLoadBC = 2 ;
analysisSettings.iniDeltaLamb = 0.1 ;


otherParams.problemName = 'prueba';
%~ otherParams.controlDofs = [ numElements+1  5 ] ;
otherParams.plotsFormat = 'vtk' ;


[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


centralNode = ceil((ncolumns/2+1))*(nrows/2+1) ;
Dof      		 = centralNode*6 - 3 	;
controlDisps =  matUs(Dof, :) 						;
loadFactors  =  loadFactorsMat(:, 2) 			;

lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;
figure
grid on
plot( controlDisps, loadFactors, 'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
