%md# Plate-element cantilever model
%md
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%md
%md## Scalars
E = 200e9 ;
nu = 0.0 ;
tz = .05 ;
q  = 1e3 ; % 1kN/m^2 
%md
%md## Numerical solution using plate elements
%md
%md### Materials
%md
materials                    = struct() ;
materials(1).hyperElasModel  = 'linearElastic' ;
materials(1).hyperElasParams = [ E nu ] ;
%md
%md### Elements
%md
elements             = struct() ;
elements(1).elemType = 'edge' ;
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
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -q 0 ] ;
%md
%md### mesh
%md
mesh = struct() ;
base_dir='';
if strcmp( getenv('TESTS_RUN'),'yes') && isfolder('examples'),
  base_dir=['.' filesep 'examples' filesep  'cantileverPlate' filesep];
end
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( [base_dir 'geometry_cantileverPlate.msh'] ) ;
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
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md

Ly = .5;
Lx = 1 ;
qlin = q*Ly;
I = Ly*tz^3/12;
analy_wmax = -qlin*Lx^4/(8*E*I) ;
numer_wmax = min(matUs(5:6:end)) ;
verifBoolean = (abs( analy_wmax - numer_wmax ) / abs(analy_wmax)) < 1e-3 

