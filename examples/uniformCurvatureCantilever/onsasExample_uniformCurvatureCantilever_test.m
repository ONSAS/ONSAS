%## Uniform curvature cantielever
%#---
%#
%#In this tutorial, the Uniform curvature cantielver example and its resolutions using ONSAS are described. The aim of this example is to validate the static corrotational 3D beam implementation by comparing the results provided by ONSAS with the analytical solution.
%#
%#The structural model is formed by one beam element member as it is shown in the figure, with the node $2$ submitted to a nodal moment $M$ and the node $1$ is restrained to linear and angular displacements in the $x-y-z$.
%#
%#![structure diagram](uniformCurvatureCantilever.svg)
%#
%#The Octave script is available at: [https://github.com/ONSAS/ONSAS.m/blob/master/examples/uniformCurvatureCantilever](turaaaaaaaaa)
%#
%#Before defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar geometry and material parameters are defined.
close all, clear all ;
addpath( [ pwd '/../../src'] );
E = 200e9 ;  nu = 0.3 ;  rho = 0 ;  l = 10 ;  ty = .1 ;  tz = .1 ; Iy = ty*tz^3/12 ;    Iz = tz*ty^3/12 ; finalTime = 1;
%#
%### MEBI parameters
%#------------------
%#
%#The modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%#### materials
%# Since the example contains only one rod the  fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
materials.hyperElasModel  = { '1DrotEngStrain'} ;
materials.hyperElasParams = { [ E nu ] } ;
%# and the parameters of this model are the Lamé parameters
%#
%### elements
%#
%#Two different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements.elemType = { 'node','beam' } ;
%# for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a rectangular-cross section by ty and tz the dimensions on y and z direction, then the elemTypeGeometry field is:
elements.elemTypeGeometry = { [], [2 ty tz ] };
elements.elemTypeParams = { [], 1 };
%#
%### boundaryConds
%#
%# The elements are submitted to two different BC settings. The nodes $A$ is completely fixed and the node $B$ is submitted by a nodal moment load. According to this, the $A$ has a constraint in displacement besides the node $B$ has an applied load. The load factor function of the moment is set so that the target load make the beam into a circle. The density is set to zero, then no inertial effects are considered.
%#
boundaryConds.loadCoordSys = { []        ; 'global'   } ;
boundaryConds.loadTimeFact = { []        ; @(t) E * Iy * 2 * pi / l /finalTime *t   } ;
boundaryConds.loadBaseVals = { []        ; [ 0 0 0 -1 0 0 ] } ;
boundaryConds.impoDispDofs = { [ 1 2 3 4 5 6 ] ; 1          } ;
boundaryConds.impoDispVals = { [ 0 0 0 ] ; 0          } ;
%#
%#
%### initial Conditions
%# homogeneous initial conditions are considered, then an empty struct is set:
initialConds                = struct() ;
%#
%### mesh parameters
%# The mesh is parmetric according to the number of elements selected with ths variable:
numElements = 10 ;
%#The coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [ (0:(Nelem))' * l/numElements'* zeros(Nelem+1,2) ] ;
%#The connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%# then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%# the following case only differs in the boundary condition
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  2   ] ;
%# the beam elements are formed by the first material, the second type of element, and no boundary conditions are applied to any element. Becouse of the parametric implmetetion depending on numElements variable an auxiliar matrix is defined:
auxiliarMatrixConec = [ (ones(numElements,1)*[ 1 2 0 0]) (1:(numElements))' (2:(numElements+1))' ] ;
%# then the conection vectors are  added to the conCell array as follows:
for i:numElements
  mesh.conecCell{ end+1,1 } = auxiliarMatrixConec (i,:) ;
end
%#
%### analysisSettings
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.deltaT        = 0.1 ;
analysisSettings.finalTime     =   finalTime ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10 ;
analysisSettings.finalTime     =   1 ;
%#
%### otherParams
otherParams.problemName = 'uniformCurvatureCantilever';
otherParams.plotParamsVector = [3] ;
otherParams.controlDofs = [Nelem+1  4 ] ;
%### Analysis case 1: NR with Rotated Eng Strain
%# In the first case ONSAS is run and the solution at the dof (angle of node B) of interest is stored:
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%# the control dof to verificate the solution is the node angle B, this corresponds to the following dof number:
angleControlDof = (Nelem+1)*6 - 2;
controlDispsNREngRot =  -matUs(angleControlDof,:) ;
loadFactorsNREngRot  =  loadFactorsMat(:,2) ;
%# and the analytical value of the load factors is computed
analyticLoadFactorsNREngRot = @(w) w * l / ( E * Iy ) ;

%## Results verification
%#
lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;
figure
plot( controlDispsNREngRot, analyticLoadFactorsNREngRot( controlDispsNREngRot) ,'b-x' , 'linewidth', lw,'markersize',ms )
hold on, grid on
plot( controlDispsNREngRot, loadFactorsNREngRot, 'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
legend('analytic','NR-RotEng','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
