%md# Static Von Mises Truss example
%md---
%md
%mdIn this tutorial, the Static Von Mises Truss example and its resolutions using ONSAS are described. The aim of this example is to validate the Newton-Raphson method implementation by comparing the results provided by ONSAS with the analytical solution. The Octave script of this example is available at [this url](https://github.com/ONSAS/ONSAS.m/blob/master/examples/staticVonMisesTruss/onsasExample_staticVonMisesTruss.m).
%md
%mdThe structural model is formed by two truss elements with length $L$ as it is shown in the figure, with the node $2$ submitted to a nodal load $P$ and restrained to move in the $x-z$ plane and nodes $1$ and $3$ fixed.
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/tutorials/StaticVonMisesTruss/vonMisesTruss.svg" alt="structure diagram" width="500"/>
%md```
%md
%mdBefore defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar auxiliar parameters are defined.
close all, clear all ;
addpath( genpath( [ pwd '/../../src'] ) );
E = 210e9 ;  A = 2.5e-3 ; ang1 = 65 ; L = 2 ; nu = 0 ;
% x and z coordinates
auxx = cos( ang1*pi/180 ) * L ;  auxz = sin( ang1*pi/180 ) * L ;
%md
%md## MEBI parameters
%md------------------
%md
%mdThe modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md### materials
%md Since both bars are formed by the same material all the fields of the `materials` struct will have only one entry. contains only one vector. The constitutive behavior is the SaintVenantKirchhoff:
materials.hyperElasModel  = { '1DrotEngStrain'} ;
materials.hyperElasParams = { [ E nu ] } ;
%md and the parameters of this model are the Young modulus and Poisson ratio.
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and truss. The nodes will be assigned in the first entry (index $1$) and the truss at the index $2$. The elemType field is then:
elements.elemType = { 'node','truss' } ;
%md for the geometries, the node has no geometry to assign (empty array), and the truss elements will be set as a square-cross section, then the elemTypeGeometry field is:
elements.elemTypeGeometry = { [], [2 sqrt(A) sqrt(A) ] };
elements.elemTypeParams = { [], 1 };
%md
%md### boundaryConds
%md
%md The elements are submitted to two different BC settings. The nodes $1$ and $3$ are fixed without applied loads (first BC), and node $2$ has a constraint in displacement and an applied load (second BC). The load factor function of the second BC is set so that the target load 1.5e8 is reached at 1 second. The density is set to zero, then no inertial effects are considered.
%md
boundaryConds.loadsCoordSys = { []        ; 'global'   } ;
boundaryConds.loadsTimeFact = { []        ; @(t) 1.5e8*t     } ;
boundaryConds.loadsBaseVals = { []        ; [ 0 0 0 0 -1 0 ] } ;
boundaryConds.imposDispDofs = { [ 1 3 5 ] ; 3          } ;
boundaryConds.imposDispVals = { [ 0 0 0 ] ; 0          } ;
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds                = struct() ;
%md
%md## mesh parameters
%md------------------
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [      0  0     0  ; ...
                       auxx  0  auxz  ; ...
                     2*auxx  0     0  ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  3   ] ;
%md the following case only differs in the boundary condition
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  2   ] ;
%md the truss elements are formed by the first material, the second type of element, and no boundary condition is applied.
mesh.conecCell{ 4, 1 } = [ 1 2 0 0  1 2 ] ;
mesh.conecCell{ 5, 1 } = [ 1 2 0 0  2 3 ] ;
%md
%md## analysisSettings
%md------------------
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        = 0.1 ;
analysisSettings.finalTime     =   1 ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10 ;
analysisSettings.finalTime     =   1 ;
%md
%md## otherParams
%md------------------
otherParams.problemName = 'staticVonMisesTruss_NR';
otherParams.plotParamsVector = [3];
otherParams.controlDofs = [2 5 ];
%md
%md## Analysis case 1: NR with Rotated Eng Strain
%md------------------
%md In the first case ONSAS is run and the solution at the dof of interest is stored .
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsNREngRot =  -matUs(11,:) ;
loadFactorsNREngRot  =  loadFactorsMat(:,2) ;
%md and the analytical value of the load factors is computed
analyticLoadFactorsNREngRot = @(w) -2 * E*A* ...
     ( (  (auxz+(-w)).^2 + auxx^2 - L^2 ) ./ (L * ( L + sqrt((auxz+(-w)).^2 + auxx^2) )) ) ...
  .* (auxz+(-w)) ./ ( sqrt((auxz+(-w)).^2 + auxx^2) )  ;
%md
%md
%md## Results verification
%md------------------
%md
%md### numerical verification
[ analyticLoadFactorsNREngRot( controlDispsNREngRot)' loadFactorsNREngRot ]
difLoad = analyticLoadFactorsNREngRot( controlDispsNREngRot)' - loadFactorsNREngRot ;
verifBoolean = ( norm( difLoad ) / norm( loadFactorsNREngRot ) ) <  1e-4
%md### Plots
lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;
figure
plot( controlDispsNREngRot, analyticLoadFactorsNREngRot( controlDispsNREngRot) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( controlDispsNREngRot, loadFactorsNREngRot, 'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
legend( 'analytic', 'NR-RotEng', 'location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('vonMisesTrussCheck.png','-dpng')
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/tutorials/StaticVonMisesTruss/vonMisesTrussCheck.png" alt="plot check" width="500"/>
%md```
%md
