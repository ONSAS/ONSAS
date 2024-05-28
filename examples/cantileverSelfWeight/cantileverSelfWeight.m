%md# Solid cantilever self-weight example
%md---
%md
%mdExample intended to validate the self-weight for tetrahedron elements.
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
% add path
addpath( genpath( [ pwd '/../../src'] ) ) ; 
% scalar parameters
E = 200e9 ; nu = 0.3 ; Lx = 2 ; Ly = .02 ; Lz = .2 ; rho = 8e3 ;
%md
%md### MEB parameters
%md
%md#### materials
%md The material of the solid considered is the Saint-Venant-Kirchhoff with Lamé parameters computed as
lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
%md since only one material is considered, a scalar struct is defined as follows
materials = struct();
materials.modelName = 'SVK' ;
materials.modelParams = [ lambda mu ] ;
materials.density = rho ;
%md
%md#### elements
%md In this model two kinds of elements are used: `tetrahedron` for the solid and `triangle` for introducing the external loads. Since two kinds of elements are used, the struct have length 2:
elements = struct();
elements(1).elemType = 'triangle' ;
elements(2).elemType = 'tetrahedron' ;
%md
%md#### boundaryConds
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [1 3 5] ;
boundaryConds(1).imposDispVals = [0 0 0] ;
%md
%md
%md### Mesh
mesh = struct();
base_dir='';
if strcmp( getenv('TESTS_RUN'),'yes') && isfolder('examples'),
  base_dir=['.' filesep 'examples' filesep  'cantileverSelfWeight' filesep];
end
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( [ base_dir 'geometry_cantileverSelfWeight.msh'] ) ;
%md#### initialConds
%md since no initial non-homogeneous initial conditions are used, an empty struct is used .
initialConds = struct() ;
%md### Analysis parameters
%md
analysisSettings = struct() ;
analysisSettings.booleanSelfWeight = true  ;
analysisSettings.methodName        = 'newtonRaphson' ;
analysisSettings.stopTolIts        = 30     ;
analysisSettings.stopTolDeltau     = 1.0e-12 ;
analysisSettings.stopTolForces     = 1.0e-12 ;
analysisSettings.finalTime          = 1   ;
analysisSettings.deltaT            = 1   ;
%md
%md### Output parameters
%otherParams.plots_format = 'vtk' ;
otherParams = struct() ;
otherParams.problemName = 'cantileverSelfWeight' ;
%md
%md Execute ONSAS and save the results:
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
%mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;
%md
%md the report is generated
outputReport( modelProperties.outputDir, modelProperties.problemName )

numericalDeflection = min( matUs(:,2) )

EI = E*Lz^3*Ly/12 ;
analyticalDeflectionLinearTheory = - ( rho*9.806*Ly*Lz * Lx^4 ) / ( 8 * EI )
verifBoolean = abs( analyticalDeflectionLinearTheory - numericalDeflection ) / abs( analyticalDeflectionLinearTheory ) < 0.03 ;
