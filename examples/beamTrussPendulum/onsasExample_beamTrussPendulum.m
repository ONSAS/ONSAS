%md## beam Truss pendulum
%md
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) );

E = 210e9 ; nu = 0.25 ;  A = 0.1^2 ;

materials.hyperElasModel  = { '1DrotEngStrain' } ;
materials.hyperElasParams = { [ E nu ] } ;
materials.density = { 8000 } ;

elements.elemType = { 'node', 'truss', 'frame' } ;

elements.elemTypeGeometry = { [], [2 sqrt(A) sqrt(A) ], [2 sqrt(A) sqrt(A) ] };
elements.elemTypeParams = { [], 1, [] };
%md
boundaryConds.loadsCoordSys = { []              ; 'global'         } ;
boundaryConds.loadsTimeFact = { []              ; @(t) 1.0         } ;
boundaryConds.loadsBaseVals = { []              ; [ 0 0 0 0 -1 0 ] } ;
boundaryConds.imposDispDofs = { [ 1 2 3 4 5 6 ] ; [ ]              } ;
boundaryConds.imposDispVals = { [ 0 0 0 0 0 0 ] ; [ ]              } ;

initialConds                = struct() ;

mesh.nodesCoords = [      0  0     0  ; ...
                          1  0     0  ; ...
                          2  0     0  ] ;

mesh.conecCell = { } ;

mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  2   ] ;
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  3   ] ;

mesh.conecCell{ 4, 1 } = [ 1 2 0 0  1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 3 0 0  2 3 ] ;

analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        = 1e-4 ;
analysisSettings.finalTime     =   1e-3 ;
analysisSettings.alphaNM     =   0.25 ;
analysisSettings.deltaNM     =   0.5 ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10 ;

otherParams.problemName = 'beamTrussPendulum_newmark' ;

%md
%md## Analysis case 1: NR with Rotated Eng Strain
%md------------------
%md In the first case ONSAS is run and the solution at the dof of interest is stored .
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
