%md# Cantilever with elasto plastic perfect constitutive model
%md
%md## Previous definitions
close all, clear all;
addpath( genpath( [ pwd '/../../src' ] ) ) ; % add ONSAS directory to path
%md
%md## MEBI parameters: Material-Element-BoundaryConditions-InitialConditions
%md
%md### Materials
%md Constitutive model parameters
E  = 210e6 ; nu = 0.3 ; % 
sigmaY = 250e3 ; 
materials(1).hyperElasModel = 'elastoPlasticPerfect' ;
materials(1).hyperElasParams = [ E, nu, sigmaY ] ;
%md### Elements
elements(1).elemType  = 'node'  ;
elements(2).elemType  = 'frame' ;
%md Section properties
ty = 0.1 ; 
tz = 0.1 ; % cross-section widths
elements(2).elemCrossSecParams = { 'rectangle'; [ ty tz ] };
elements(2).massMatType     =  'consistent' ; 
%md### BoundaryConditions
% Supports
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
% Loads
%md Applied nodal loads
P = 1 ; % applied nodal load
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -P 0 ] ;
%md### InitialConditions
%md empty struct
initialConds = struct() ;
%md
%md## Mesh
%md Mesh nodes
L = 1 ; %
nnodesMesh = 21 ;
xcoords = linspace(0,L,nnodesMesh)' ;
ycoords = zeros(length(xcoords),1) ;
zcoords = zeros(length(xcoords),1) ;
 
mesh.nodesCoords = [ xcoords ycoords zcoords ] ;

%md
%md Conec Cell
mesh.conecCell = { } ;
% Auxiliar
vec1 = (1:1:(length(xcoords)-1))' ;
vec2 = (2:1:(length(xcoords)))' ;
vec3 = ones(length(vec1),1) ;
loadedNode = length(xcoords) ;
%md nodes
mesh.conecCell{1, 1 } = [ 0 1 1 0  1 ] ;
mesh.conecCell{2, 1 } = [ 0 1 2 0  loadedNode ] ;
%md and frame elements
for i=1:(nnodesMesh-1) 
	mesh.conecCell{i+2,1} = [ 1 2 0 0 i i+1 ] ;
end
%md
%md Analysis settings
analysisSettings = struct() ;

otherParams.problemName = 'cantileverEP' ;
otherParams.plotsFormat = 'vtk' ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;








