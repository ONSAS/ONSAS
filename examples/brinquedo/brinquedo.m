%md# Plate-element cantilever model
%md
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%md
%md## Scalars
E = 200e9 ;
nu = 0.0 ;
tz = .05 ;
%md
%md## Numerical solution using plate elements
%md
%md### Materials
%md
materials                    = struct() ;
materials(1).modelName  = 'elastic-linear' ;
materials(1).modelParams = [ E nu ] ;
%md
%md### Elements
%md
elements             = struct() ;
elements(1).elemType = 'node' ;
elements(2).elemType = 'triangle-plate' ;
elements(2).elemCrossSecParams = {'thickness', tz } ;
%md
%md### Boundary conditions
%md
boundaryConds                  = struct() ;
boundaryConds(1).imposDispDofs =  [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 0 0 0 ] ;
%
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) t  ;
boundaryConds(2).loadsBaseVals = [1 0 0 0 -1 0 ] ;
%
%md
%md### mesh
%md
mesh = struct() ;
mesh.nodesCoords = [ 0 0 0 ; 1 0 1; 1 1 1; 0 1 1 ];
mesh.conecCell = {
    [ 0 1 1 1 ];
    [ 0 1 1 4 ];
    [ 0 1 2 2 ];
    [ 0 1 2 3 ];
    [ 1 2 0 1 2 3 ];
    [ 1 2 0 2 3 4 ];
}
;

%md### Initial conditions
initialConds                  = struct() ;
%md
%md#### Analysis settings
analysisSettings               = struct() ;
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   1   ;
analysisSettings.finalTime     =   1  ;
analysisSettings.stopTolDeltau =   1e-10    ;
analysisSettings.stopTolForces =   1e-10    ;
analysisSettings.stopTolIts    =   10      ;
%md
%md#### OtherParams
%md The nodalDispDamping is added into the model using:
otherParams                  = struct() ;
%md The name of the problem is:
%md
otherParams.problemName  = 'cantileverPlate' ;
otherParams.plots_format = 'vtk' ;
%md
%md Execute ONSAS and save the results:
[ modelInitSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
%mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelInitSol, modelProperties, BCsData ) ;
%md
%md## verification
