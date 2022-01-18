%md# Solid cantilever self-weight example
%md---
%md
%mdExample intended to validate the self-weight for tetrahedron elements.
clear all, close all
% add path
addpath( genpath( [ pwd '/../../src'] ) ) ;
% scalar parameters
E = 200e9 ; nu = 0.3 ; Lx = 2 ; Ly = .02 ; Lz = .2 ; rho = 8e3 ;
%md
%md### MEBI parameters
%md
%md#### materials
%md The material of the solid considered is the Saint-Venant-Kirchhoff with Lamé parameters computed as
lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
%md since only one material is considered, a scalar struct is defined as follows
materials.hyperElasModel = 'SVK' ;
materials.hyperElasParams = [ lambda mu ] ;
materials.density = rho ;
%md
%md#### elements
%md In this model two kinds of elements are used: `tetrahedron` for the solid and `triangle` for introducing the external loads. Since two kinds of elements are used, the struct have length 2:
elements(1).elemType = 'triangle' ;
elements(2).elemType = 'tetrahedron' ;
%md
%md#### boundaryConds
boundaryConds(1).imposDispDofs = [1 3 5] ;
boundaryConds(1).imposDispVals = [0 0 0] ;
%
%md
%md#### initialConds
%md since no initial non-homogeneous initial conditions are used, an empty struct is used .
initialConds = struct();
%md
%md### Mesh
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'geometry_cantileverSelfWeight.msh' ) ;

%md### Analysis parameters
%md
analysisSettings.booleanSelfWeight = true  ;
analysisSettings.methodName        = 'newtonRaphson' ;
analysisSettings.stopTolIts        = 30     ;
analysisSettings.stopTolDeltau     = 1.0e-12 ;
analysisSettings.stopTolForces     = 1.0e-12 ;
analysisSettings.finalTime          = 1   ;
analysisSettings.deltaT            = 1   ;
%md
%md### Output parameters
otherParams.plotsFormat = 'vtk' ;
otherParams.problemName = 'cantileverSelfWeight' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


numericalDeflection = min( matUs(:,2) )

EI = E*Lz^3*Ly/12 ;
analyticalDeflectionLinearTheory = - ( rho*9.806*Ly*Lz * Lx^4 ) / ( 8 * EI )
verifBoolean = abs( analyticalDeflectionLinearTheory - numericalDeflection ) / abs( analyticalDeflectionLinearTheory ) < 0.03
