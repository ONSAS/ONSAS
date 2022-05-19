%md# Chain example
%md---
%md
%mdIn this tutorial, the Chain example and its resolution using ONSAS are described. The aim of this example is to demonstrate how to add non initial homogenous conditions and contrast the output by modelling a chain with beam and truss elements loaded by self-weight only . The chain span length $L$ and is fixed at the edges (node $1$ and $2$) constraining linear and angular displacements. A fixed framework it is defined at node $1$ with a axis <span style="color:blue"> $x$ </span> and <span style="color:blue"> $z$, where the gravity is applied along $-z$. 
%md Three different configurations are defined: reference (<span style="color:red"> *Ref.Config* </span>.), initial(<span style="color:grey"> *Init.Config* </span>.), and deformed (<span style="color:black"> *Def.Config* </span>.). The reference configuration is set considering an horizontal chain, on the other hand the  V-shape initial configuration is set rotating $45$ degrees from the reference configuration and finally during the time motion the deformed configuration can be found using the displacements field measured from the reference configuration.    

%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS.docs/master/docs/src/diagramChain.svg" alt="structure diagram" width="500"/>
%md```
%md## Analytic solution
%md--------------------
%md
%md
%md## Numerical solution
%md---------------------
%md
%md
%mdThe Octave script of this example is available at [this url](https://github.com/ONSAS/ONSAS.m/blob/master/examples/chain/onsasExample_chain.m).
%md
%mdBefore defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar auxiliary parameters are defined.
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
% material scalar parameters
E  = 210e3 ; nu  = 0 ; rho = 8050 ;
% gemotric scalar parameters
rotAng = deg2rad(20) ; L = 2 ; b = 0.05 ; % m - width of square section
% the number of elements of the mesh
numElements = 10  ; %must be greater than and a pair value 
%md
%md## Numerical solution: truss case
%md---
%md
%md### MEBI parameters
%md
%mdThe modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md#### materials
%md Since both bars are formed by the same material only one `materials` struct is defined. The constitutive behavior considered in the first analysis case is the Rotated Engineering strain, then the field `hyperElasModel` is set to:
materials.hyperElasModel  = '1DrotEngStrain' ;
%md and in the field `hyperElasParams` a vector with the parameters of the Engineering Strain model is set
materials.hyperElasParams = [ E nu ] ;
%md which in the case of this model are the Young modulus and the Poisson ratio. Then, the material density is set using `denisity` field of materials struct using
materials.density = rho ;
%md
%md#### elements
%md
%mdTwo different types of elements are required to create the truss model: `node` and `truss`, thus, the _elements_ struct will have two entries. The first entry is a node so
elements(1).elemType = 'node' ;
%md and the second a truss
elements(2).elemType = 'truss';
%md for the geometries, the node has no geometry to assign, and the truss elements will be set as a square-cross section, then the elemCrossSecParams field is:
elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [b b] ;
elements(2).massMatType = 'consistent' ;
%md
%md#### initial Conditions
%md homogeneous initial conditions are considered, by the V-shape depicted in the diagram above, consequently the coordinates mathematical expression of the initial configuration shape is built for 0 and l/2 and flipping that vector from l/2 to 0.
yInitConfigCoordsMiddle = linspace(0, L, numElements)(1:end/2 +1 ) / cos(rotAng) ; 
yInitConfigCoords = [ yInitConfigCoordsMiddle flip( yInitConfigCoordsMiddle(1:end-1) ) ] ;
%md since only non homogenous initial conditios are added into $y$ axis the associated degrees of freedom are
dofsYInitCond = ( 3:6:6*(numElements +1) );
% the number of different initial conditions imposed are:
numDifInitialConds = floor(sum(dofsYInitCond>0) /2) ; 
%md first create an empty `initialConds` struct
initialConds = struct() ;
%md the format code to add non homogeneous initial condition is a vector that $i-th$ position contain the degree of freedom and in the following column the value of the initial condition : [dof1 valueInitCond1; dof2 valueInitCond2]. In this case
initialConds = {} ;
for indCond = 1:numDifInitialConds
  initialConds(indCond).nonHomogeneousUDofs = [ 3 ] ;
  initialConds(indCond).nonHomogeneousUVals = [yInitConfigCoordsMiddle(indCond+1) ] ;
end
%md
%md
%md#### boundaryConds
%md
%md The elements are submitted to only one kinematic boundary condition (BC), then the struct _boundaryConds_ will have length one.
%md The nodes $1$ and $3$ are fixed, without loads applied. Since the gravity is going to be added by a boolean later. 
%md So the fixed displacements boundary condition is typed
%md
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
%md
%md### mesh parameters
%mdThe coordinates of the nodes in the reference configuration are given by the matrix:
mesh.nodesCoords = [ ( 0:(numElements) )' * L / numElements  zeros(numElements+1,2) ] ;
%md where the columns 1,2 and 3 correspond to $x$ (a equally spaced vector form 0 to L), $y$ and $z$ null coordinates, respectively, and the row $i$-th corresponds to the coordinates of node $i$.
%md
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is created
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the `elements` struct), and the first entry of the cells of the boundary conditions struct. Fixed nodes MEBI data is set as
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1  ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  numElements + 1 ] ;
%md since all nodes have the first not homogenous initial condition and no material, and boundary condition then
%md the vector of initial conditions indexes is:
initCondVec = [ [1:numDifInitialConds] flip( [ 1:numDifInitialConds-1] ) ] ;
for numIC =  1:length(initCondVec)
  mesh.conecCell{ numIC+2, 1 } = [ 0 1 0  initCondVec(numIC)  numIC+1  ] ;
end
%md the truss elements are formed by the first material, the second type of element, and no boundary conditions are applied to any element, consecutive to the nods data into conecCell we add
for i = 1:numElements
  mesh.conecCell{ numElements + 1 + i,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
%mdAn alpha HHT method algorithm is used to solve the problem so the `methodName` field of _analysisSettings_ struct is set
analysisSettings.methodName    = 'alphaHHT' ;
%md and the following parameters correspond to the iterative numerical analysis settings are
analysisSettings.deltaT        =   0.01 ;
analysisSettings.finalTime     =   5    ;
analysisSettings.stopTolDeltau =   0    ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   30   ;
%md
%Finally the self weight boolean is activated
analysisSettings.booleanSelfWeight = true ;
### otherParams
%mdA name problem and vtk output is introduced as:
otherParams.problemName = 'trussChain'; 
otherParams.plotsFormat = 'vtk'       ;
%mdONSAS code is executed for the input structs aforementioned explained
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;