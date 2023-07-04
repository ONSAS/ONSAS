%md# Cantilever slab example
%md
%md We start as all models, clearing the workspace and adding the ONSAS path to the work path.
% clear workspace and add path
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%md
%md The following numeric parameters are considered.
% scalar parameters for spring-mass system
%md
E = 200e9 ;
nu = 0.3 ;
Lx = 2 ;
Ly = 2 ;
tz = .01 ;
%md
%md## Numerical solution
%md
%md### Materials
%md
materials                    = struct() ;
materials(1).hyperElasModel  = 'linearElastic' ;
materials(1).hyperElasParams = [ E nu ]          ;
%md
%md### Elements
%md
%md In this case only `'node'` and  `'truss'` elements are considered and the lumped inertial formulation is set for the truss element:
elements             = struct() ;
elements(1).elemType = 'node'                                 ;
elements(2).elemType = 'triangle-plate'                             ;
elements(2).elemCrossSecParams = {'thickness', tz } ;
%md
%md### Boundary conditions
%md
boundaryConds                  = struct() ;
boundaryConds(1).imposDispDofs =  [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 0 0 0 ] ;

boundaryConds(2).loadsCoordSys = 'global'                   ;
boundaryConds(2).loadsTimeFact = @(t) t  ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -1e3 0 ]            ;
%md
%md### mesh
%md Only two nodes are considered so the nodes matrix is:
mesh             = struct() ;
mesh.nodesCoords = [  0   0   0 ; ...
                      Lx  0   0 ; ...
                      Lx  Ly  0 ; ...
                      0   Ly  0 ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1,     1 } = [ 0 1 1    1   ] ;
mesh.conecCell{ end+1, 1 } = [ 0 1 1    2   ] ;
mesh.conecCell{ end+1, 1 } = [ 0 1 2    3   ] ;
mesh.conecCell{ end+1, 1 } = [ 0 1 2    4   ] ;

mesh.conecCell{ end+1, 1 } = [ 1 2 0    1 2 3   ] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0    1 3 4   ] ;

%md### Initial conditions
initialConds                  = struct() ;
%md
%md#### Analysis settings
analysisSettings               = struct() ;
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   0.5   ;
analysisSettings.finalTime     =   1  ;
analysisSettings.stopTolDeltau =   1e-10    ;
analysisSettings.stopTolForces =   1e-10    ;
analysisSettings.stopTolIts    =   15      ;
%md
%md#### OtherParams
%md The nodalDispDamping is added into the model using:
otherParams                  = struct() ;
%md The name of the problem is:
%md
otherParams.problemName = 'cantileverSlab'     ;
%md
%md Execute ONSAS and save the results:
[matUsNewmark, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
