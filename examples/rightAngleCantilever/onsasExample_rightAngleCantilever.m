%md# Right angle cantilever beam problem
%md---
%md
%mdIn this tutorial, the right angle cantilever beam problem example and its resolution using ONSAS are described. The aim of this example is to validate the dtnamic co-rotational 3D beam implementation by comparing the results provided by ONSAS with solution provided by (T.L Lee & J.M Battini, 2014) and originally solved in (Simo & Vu-Quoc, 1998).  The Octave script of this example is available at [this url](https://github.com/ONSAS/ONSAS.m/blob/master/examples/rightAngleCantilever/onsasExample_rightAngleCantilever.m).
%md
%mdIt is conformed by two identical right-angled bars, where each memberhas a length of $L = 10$ m. The structure is embedded at the base and a force in z direction is applied at the elbow. This force bends nd troses the system into the xy plane, producing free vibrations of wide amplitude. These force acts during two initial seconds, increases linearly until the first second of simulation and then decreases to zero. 
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/cantileverBeam_HTML.svg" alt="structure diagram" width="500"/>
%md```
%md Before defining the input structs all variables are clenad, the figure closed and the src folder added to the workspace path:
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
% material scalar parameters
%mdThe material properties must comply certain equals $GA = EA = 10^6$ and $GJ = EI = 10^3$, therefore the choice of these is obtained synthetically by solving an indeterminate compatible system. For this work the second moments of inertia along the axis z and y and the values of the linear and transverse modulus of elasticity were chosen are:
%md```math
%md E = G = 10^6 A = 1 I = J = 10^−3,
%md```
%mdthese values satisfy the aformentioned equations: 
% material parameters:
E = 1e6 ;  nu = -0.5 ; rho = 1 ;
% geometrical scalar parameters
I = 1e-3 ; J = I ; A = 1 ; L = 10;
%md In addition the dyadic tensor if inertia is 
%md```math
%md I_rho = diag(20,10,10),
%md```
Irho = diag([20 10 10],3,3)
% the number of elements of each meber:
nElemsPerBeam = 10 ;
%md
%md##Numerical solution
%md### MEBI parameters
%md
%mdThe modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md### materials
%md Since the example contains only one type of material the fields of the `materials` struct will have only one entry. Morover, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
materials.density = rho ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the sutrcut geometries, the node has not geometry thus an empty array is assigned, and the truss elements will be set as a synthetical-cross section with properties stated above, subsequently the elemTypeGeometry field is:
elements(2).elemTypeGeometry = [1 A J I I Irho(1,1) Irho(2,2) Irho(3,3)] ;
elements(2).elemTypeParams   = 1;
%md
%md### boundaryConds
%md
%md The elements are submitted to two different BC settings. The first BC corresponds to a welded condition (all 6 dofs set to zero):
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%mdand the second corresponds to an impact nodal force at the L joint, where the target load produces a circular form of the deformed beam.
boundaryConds(2).loadsCoordSys = 'global' ;
constF = 50 ;
boundaryConds(2).loadsTimeFact = @(t) constF*t*(t<1) + (constF*t)*(t>=1)*(t<2) + 0;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ] ;
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds                = struct() ;
%md
%md### mesh parameters
%mdThe coordinates of the nodes are computed by using axuilar coordinates vector that represent the local nodes cooridante if the origin axis where fixed at the origin of each member. Then the mesh are given by the matrix:
auxCoords     = linspace( 0, L, nElemsPerBeam+1 )' ;
mesh.nodesCoords = [ zeros(nElemsPerBeam+1,1)       auxCoords               zeros(nElemsPerBeam+1,1) ; ...
                        -auxCoords(2:end)         ones(nElemsPerBeam,1)*L   zeros(nElemsPerBeam  ,1) ] ;

%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%md the following case only differs in the boundary condition (not displacemnt imposed) at the node 2.
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  nElemsPerBeam+1 ] ;
%md the beam elements are formed by the first material, the second type of element, and no boundary conditions are applied to any element.
for i = 1 : 2*nElemsPerBeam,
  mesh.conecCell{ i+2,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
analysisSettings.deltaT        =    0.01  ;
analysisSettings.finalTime      =   3 ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   30   ;

analysisSettings.methodName    = 'newmark' ;
analysisSettings.alphaNM      =   0.25   ;
analysisSettings.deltaNM      =   0.5   ;
%md
%md## otherParams
otherParams.problemName = 'rightAngleCantilever';
otherParams.controlDofs = [ nElemsPerBeam+1  5 ] ;
otherParams.plotsFormat = 'vtk' ;
%md In the first case ONSAS is run and the solution at the dof (angle of node B) of interest is stored:
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md the control dof to verificate the solution is the displacments where the load is applyed, this corresponds to the following dof number:
angleControlDof      = (numElements+1)*6 - 1;
controlDispsNREngRot =  -matUs(angleControlDof,:) ;
loadFactorsNREngRot  =  loadFactorsMat(:,2) ;