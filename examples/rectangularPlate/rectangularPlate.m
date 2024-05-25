%md# Rectangular Plate test
%md
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%md Scalar parameters
E = 1e3 ;
nu = 0.0 ;
tz = .1 ;
q = 1 ;

materials = struct() ;
materials.modelName  = 'elastic-linear' ;
materials.modelParams =  [ E nu ]       ;

elements = struct() ;
elements(1).elemType           = 'edge'    ;
elements(1).elemCrossSecParams = tz         ;
elements(2).elemType           = 'triangle-plate';
elements(2).elemCrossSecParams = {'thickness', tz } ;

% -----------------------------------------------------
boundaryConds                  = struct() ;

boundaryConds(1).imposDispDofs =  [ 1 3 5 ]                  ; % fixed nodes: 1 2 4
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;

boundaryConds(2).loadsCoordSys = 'global'                   ;
boundaryConds(2).loadsTimeFact = @(t) t  ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -q 0 ]            ; % forces in node 1

% -----------------------------------------------------
base_msh='';
if strcmp( getenv('TESTS_RUN'),'yes') && isfolder('examples'),
  base_msh=['.' filesep 'examples' filesep 'rectangularPlate' filesep];
end
mesh = struct();
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( [ base_msh 'rectangularPlate.msh'] ) ;

% -----------------------------------------------------
initialConds = struct();
%md
% -----------------------------------------------------
%md#### Analysis parameters
%md
%md The Newton-Raphson method is employed to solve 2 load steps. The ratio between `finalTime` and `deltaT` sets the number of load steps used to evaluate `boundaryConds(3).loadsTimeFact` function:  
analysisSettings = struct() ;
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-12 ;
analysisSettings.stopTolForces = 1.0e-12 ;
analysisSettings.finalTime     = 1       ;
analysisSettings.deltaT        = 1      ;
%md
%md#### Output parameters
%md
otherParams = struct() ;
otherParams.problemName = 'rectangularPlate' ;
otherParams.plots_format = 'vtk' ;
%md The ONSAS software is executed for the parameters defined above and the displacement solution of each load(time) step is saved in `matUs`matrix:
%md
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
%mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

Lx = max( mesh.nodesCoords(:,1) ) ;
Ly = max( mesh.nodesCoords(:,2) ) ;
I = Lx*tz^3/12 ;

q1D = q * Lx ;

theta_max_analy =   q1D*Ly^3/(24*E*I ) ;
delta_max_analy = 5*q1D*Ly^4/(384*E*I) ;

theta_max_numer = abs( max(matUs(2:6:end,2)) ) ;
delta_max_numer = abs( min(matUs(5:6:end,2)) ) ;

verifBoolean = ( ( theta_max_analy - theta_max_numer ) < 1e-2*theta_max_analy ) && ...
               ( ( delta_max_analy - delta_max_numer ) < 1e-2*delta_max_analy ) ;

