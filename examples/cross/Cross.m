% Orthogonal grid

close all, clear all, close all
% add path
addpath( genpath( [ pwd '/../../src'] ) );

% material scalar parameters
E = 30e6 ; % kN/m2
nu = 0.2 ;
% mesh
ty = .2 ;
tz = .1 ;



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
elements(2).elemTypeGeometry = [2 ty*.5 tz*.5 ] ;
elements(2).elemTypeParams   = 1          ;

% Rows
elements(3).elemTypeGeometry = [2 ty/3 tz/3 ] ;
elements(3).elemTypeParams   = 1          ;


% Boundary conditions
% ----------------------------------------------------------------------
% Pinned support
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
% Roller support
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -10 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) 1     ;

boundaryConds(3).loadsCoordSys = 'global'        ;
boundaryConds(3).loadsBaseVals = [ 0 0 0 0 -20 0 ] ;
boundaryConds(3).loadsTimeFact = @(t) t     ;

% Initial conditions
% ----------------------------------------------------------------------
initialConds                = struct() ;

% Mesh

% Nodes coords
% ----------------------------------------------------------------------
mesh.nodesCoords = [ -1   0  0  ; ...
                      0   0  0  ; ...
                      1   0  0  ; ...
                      0  -1  0  ; ...
                      0   1  0  ] ;

% Conec cell
% ----------------------------------------------------------------------
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1 ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0   3 ] ;

mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 2 ] ;
mesh.conecCell{ 4, 1 } = [ 1 2 0 0   2 3 ] ;

mesh.conecCell{ 5, 1 } = [ 1 2 0 0   4 2 ] ;
mesh.conecCell{ 6, 1 } = [ 1 2 0 0   2 5 ] ;

mesh.conecCell{ 7, 1 } = [ 0 1 3 0   5 ] ;


% Analysis settings
% ----------------------------------------------------------------------
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   .5  ;
analysisSettings.finalTime      =   1   ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   10   ;

otherParams.problemName = 'cross_vtk_test';
otherParams.plotsFormat = 'vtk' ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

controlDisps =  matUs( 6+6, :)

loadFactors  =  loadFactorsMat(:, 2) 			;
stop
lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;
figure
grid on
plot( controlDisps, loadFactors, 'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
