%md# Aerodynamic linear static cantilever beam example
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
addpath( genpath( [ pwd ] ) );
% material scalar parameters
E = 70e9 ;  nu = 0.3 ; rho = 700 ;
% geometrical scalar parameters
l = 20 ; dext = .5 ;  b = 1e-3  ; dint  = dext - 2*b ;
A = pi * (dext^2 - dint^2) / 4  ;
J = pi * (dext^4 - dint^4) / 32 ; Iyy = J/2 ; Izz = Iyy ;
Irho = diag([J I I],3,3);
% the number of elements of the mesh
numElements = 1 ;
%md##Numerical solution
%md### MEBI parameters
%md
%md### materials
%md Since the example contains only one rod the fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
materials.density = rho ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a rectangular-cross section with $t_y$ and $t_z$ cross-section dimensions in $y$ and $z$ directions, then the elemTypeGeometry field is:
elements(2).elemTypeGeometry = [1 A J I I Irho(1,1) Irho(2,2) Irho(3,3)] ;
elements(2).elemTypeAero = [0 dext 0] ;
%md
%md### boundaryConds
%md
%md The elements are submitted to two different BC settings. The first BC corresponds to a welded condition (all 6 dofs set to zero)
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsTimeFact = @(t) 0 *t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 -1 0 0 ] ;
%md the name of the wind velocity function is: 
boundaryConds(3).userWindVel = 'windVel'
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds                = struct() ;
%md
%md### mesh parameters
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%md the following case only differs in the boundary condition and the node number
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        =   0.1  ;
analysisSettings.finalTime     =   0.1  ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;
%md
%md## otherParams
otherParams.problemName = 'aeroLinStaticCantilever';
otherParams.controlDofs = [ numElements+1  4 ] ;
otherParams.plotsFormat = 'vtk' ;
%md## Analysis case 1: NR with Rotated Eng Strain
%md In the first case ONSAS is run and the solution at the dof (angle of node B) of interest is stored:
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
