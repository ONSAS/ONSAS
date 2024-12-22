% md# Plate patch test
% md
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
% md
% md
E = 1e3 ;
nu = 0.3 ;
a = 20;
b = 10;
tz = 1 ;
% md
% md## Numerical solution
% md
% md### Materials
% md
materials                    = struct() ;
materials(1).modelName  = 'elastic-linear' ;
materials(1).modelParams = [ E nu ]          ;
% md
% md### Elements
% md
elements             = struct() ;
elements(1).elemType = 'node'                                 ;
elements(2).elemType = 'triangle-plate'                             ;
elements(2).elemCrossSecParams = {'thickness', tz } ;
% md
% md### Boundary conditions
% md
boundaryConds                  = struct() ;
boundaryConds(1).loadsCoordSys = 'global'                   ;
boundaryConds(1).loadsTimeFact = @(t) t  ;
boundaryConds(1).loadsBaseVals = [ 0 a 0 -b 0 0 ]            ; % forces in node 1
boundaryConds(1).imposDispDofs =  [ 1 3 5 ]                  ; % fixed nodes: 1 2 4
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;

boundaryConds(2).loadsCoordSys = 'global'                   ;
boundaryConds(2).loadsTimeFact = @(t) t  ;
boundaryConds(2).loadsBaseVals = [ 0 a 0 b 0 0 ]            ; % forces in node 2
boundaryConds(2).imposDispDofs =  [ 1 3 5 ] ; % fixed nodes: 1 2 4
boundaryConds(2).imposDispVals =  [ 0 0 0 ] ;

boundaryConds(3).loadsCoordSys = 'global'                   ;
boundaryConds(3).loadsTimeFact = @(t) t  ;
boundaryConds(3).loadsBaseVals = [ 0 -a 0 b -2 0 ] ; % forces in node 3

boundaryConds(4).loadsCoordSys = 'global'                   ;
boundaryConds(4).loadsTimeFact = @(t) t  ;
boundaryConds(4).loadsBaseVals = [ 0 -a 0 -b 0 0 ]            ;
boundaryConds(4).imposDispDofs =  [ 1 3 5 ] ; % fixed nodes: 1 2 4
boundaryConds(4).imposDispVals =  [ 0 0 0 ] ;

% md
% md### mesh
% md Only two nodes are considered so the nodes matrix is:
mesh             = struct() ;
mesh.nodesCoords = [  0     0     0 ; ...
                      2*a   0     0 ; ...
                      2*a   2*b   0 ; ...
                      0     2*b   0 ; ...
                      3/4*a 3/2*b 0 ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1,     1 } = [ 0 1 1    1   ] ;
mesh.conecCell{ end+1, 1 } = [ 0 1 2    2   ] ;
mesh.conecCell{ end+1, 1 } = [ 0 1 3    3   ] ;
mesh.conecCell{ end+1, 1 } = [ 0 1 4    4   ] ;

mesh.conecCell{ end+1, 1 } = [ 1 2 0    1 2 5 ] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0    2 3 5 ] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0    3 4 5 ] ;
mesh.conecCell{ end+1, 1 } = [ 1 2 0    1 5 4 ] ;

% md### Initial conditions
initialConds                  = struct() ;
% md
% md#### Analysis settings
analysisSettings               = struct() ;
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   0.5   ;
analysisSettings.finalTime     =   1  ;
analysisSettings.stopTolDeltau =   1e-10    ;
analysisSettings.stopTolForces =   1e-10    ;
analysisSettings.stopTolIts    =   15      ;
% md
% md#### OtherParams
% md The nodalDispDamping is added into the model using:
otherParams                  = struct() ;
% md The name of the problem is:
% md
otherParams.problemName = 'platePatchTest'     ;
otherParams.plots_format = 'vtk'     ;
% md
% md
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

numeric_w3 = matUs(6*2+5, end);
numeric_w5 = matUs(6*4+5, end);

verifBoolean = abs( numeric_w3 - (-12.48) ) < 1e-4 