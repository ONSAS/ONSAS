%md## beam Truss pendulum
%md
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) );
E = 210e9 ; nu = 0.25 ;  A = 0.1^2 ;
%
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
materials.density = 8000 ;
%
elements(1).elemType = 'node' ;
%
elements(2).elemType = 'truss' ;
elements(2).elemTypeGeometry = [2 sqrt(A) sqrt(A) ] ;
elements(2).elemTypeParams = 1 ;
%
elements(3).elemType = 'frame' ;
elements(3).elemTypeGeometry = [2 sqrt(A) sqrt(A) ] ;
%md
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%
boundaryConds(2).loadsCoordSys = 'global'         ;
boundaryConds(2).loadsTimeFact = @(t) 1.0         ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -1 0 ] ;
%
initialConds                = struct() ;
%
mesh.nodesCoords = [      0  0     0  ; ...
                          1  0     0  ; ...
                          2  0     0  ] ;
%
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  2   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  3   ] ;
%
mesh.conecCell{ 4, 1 } = [ 1 2 0 0  1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 0  2 3 ] ;
%
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   1e-4 ;
analysisSettings.finalTime      =   1e-2 ;
analysisSettings.alphaNM       =   0.25 ;
analysisSettings.deltaNM       =   0.5  ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;
%
otherParams.problemName = 'beamTrussPendulum_newmark' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

figure
plot( matUs(6+5,:) )
