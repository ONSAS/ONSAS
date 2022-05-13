%md# Right angle cantilever beam problem
%md---
%md
%mdIn this tutorial, the right angle cantilever beam problem example and its resolution using ONSAS are described. The aim of this example is to validate the dynamic co-rotational 3D beam implementation by comparing the results provided by ONSAS with a solution provided by [(T.L Lee & J.M Battini, 2014)](https://www.sciencedirect.com/science/article/abs/pii/S0045782513003022) and originally solved in (Simo & Vu-Quoc, 1998). The Octave script of this example is available at [this url](https://github.com/ONSAS/ONSAS.m/blob/master/examples/rightAngleCantilever/onsasExample_rightAngleCantilever.m). This is a validation classical example, in corrotational literature and all units discribed in this problem are meaningless in real terms. 
%md
%mdThe example is conformed by two identical right-angled bars, where each member has a length of $L = 10$. The structure is embedded at the base and a force in z direction is applied at the elbow. This force bends nd troses the system into the $x-y$ plane, producing free vibrations of wide amplitude. This force acts during two initial seconds, increases linearly until the first second of simulation and then decreases to zero. 
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/cantileverBeam_HTML.svg" alt="structure diagram" width="500"/>
%md```
%mdBefore defining the input structs all variables are cleaned, the open windows are closed and the source folder added to the workspace path:
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
% material scalar parameters
%mdThe material and geometrical properties must comply certain equals:
%md```math
%md GA = EA = 10^6, $GJ = EI = 10^3,
%md```
%mdand
%md```math
%md I_rho = diag(20,10,10),
%md```
%mdtherefore these magnitudes are synthetically set by solving an indeterminate compatible system. For this work the input properties are:
% material parameters:
E = 1e6 ;  nu = -0.5 ; rho = 1 ;
% geometrical scalar parameters:
I = 1e-3 ; J = I ; A = 1 ; L = 10;
%mdIn addition the dyadic tensor if inertia is 
Irho = diag([20 10 10],3,3);
%and the number of elements of each meber are:
nElemsPerBeam = 10 ;
%md
%md##Numerical solution
%md### MEBI parameters
%md
%mdThe modelling of the structure begins with the definition of the material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md### materials
%mdSince the example contains only one type of material the fields of the `materials` struct will have only one entry. The material density has no effect in this case because no gravity load is included and the dyadic inertia tensor is set manually. Moreover, the constitutive behavior considered is the Rotated Engineering strain, thus the field `hyperElasModel` is:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
materials.density = rho ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The `elemType` field is then:
elements(1).elemType = 'node'     ;
elements(2).elemType = 'frame'    ;
elements(2).massMatType = 'consistent' ;

%mdIn order to add the struct of geometry the assign to the node is an empty input (because it has not geometrical properties), and the truss elements will be set as with synthetical cross section with properties stated above, subsequently the `elemCrossSecParams` field is the:
elements(2).elemCrossSecParams{1,1} = 'generic' ;
elements(2).elemCrossSecParams{2,1} = [A J I I Irho(1,1) Irho(2,2) Irho(3,3)] ;
%md
%md### boundaryConds
%md
%mdThe elements are submitted to two different BC settings. The first BC corresponds to a welded condition (all 6 dofs of the node are equal to zero during all time span) then the `boundaryConds(1)` set is:
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%mdand the second BC corresponds to an impact nodal force at the L-joint, where the target load produces a triangular force in $z$ direction, reaching $50$ N in $t=1$ s and then decrease to zero in $t=2$ s. So `boundaryConds(2)` is:
boundaryConds(2).loadsCoordSys = 'global' ;
constF = 50 ;
boundaryConds(2).loadsTimeFact = @(t) constF*t*(t<1) + (2*constF - constF*(t))*(t>=1)*(t<2);
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ] ;
%md
%md### initial Conditions
%mdNon homogeneous initial conditions are not considered in this example, consequently the `initialConds`  struct is set empty:
initialConds                = struct() ;
%md
%md### mesh parameters
%mdThe coordinates of the nodes are computed by using auxiliary coordinate vector that represents the local nodes coordinate, if the origin axis where fixed at the origin of each member. Thereafter the mesh is given by the following node's matrix:
auxCoords     = linspace( 0, L, nElemsPerBeam+1 )' ;
mesh.nodesCoords = [ zeros(nElemsPerBeam+1,1)       auxCoords               zeros(nElemsPerBeam+1,1) ; ...
                        -auxCoords(2:end)         ones(nElemsPerBeam,1)*L   zeros(nElemsPerBeam  ,1) ] ;

%mdThe connectivity is introduced using the `conecCell`. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of each element (defined in node connectivity). First an empty cell is initialized:
mesh.conecCell = { } ;
%mdthen the nodes with BC conditions are defined, both with material zero (since nodes dont have material). Node _1_ is a node element type so the first dimenssion of the 'elements' struct [0 1< ...] is assigned. This node has null displacemnt and no load imposed, this is the first dimensions of the `boundaryConds` struct [0 1 1< ...], where each array of this dimension contain the loads and displacement BC. Also any non-homogeneous initial condition is considered (then zero is used) [0 1 1 0< ...]. This node attributes are loaded in the first entry {1<,..} of _conecCell_ inside first dimension of `mesh` struct {1,1<}, finally the node is assigned:
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
%mdthe following node only differs in the boundary condition (not displacement and load imposed) this BC is located at dimension 2 of `boundaryConds` struct [0 1 2<...], so the node _2_ is set in the second cell entry {2,...} of the first dimension mesh struct {2,1}:
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  nElemsPerBeam+1 ] ;
%mdthe beam elements are formed by the first material [ 1< ...], the second type of element[ 1 2<...] and no boundary and non homogeneous initial conditions are applied   [ 1 2 0< 0<  ...] are applied to any element. Also each element connect consecutive nodes [ 1 2 0 0 i< i+1< ]
for i = 1 : 2*nElemsPerBeam,
  mesh.conecCell{ i+2,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
%mdA aplha-HHT algorithm is used to solve this problem with the following parameters during $30$ s: 
analysisSettings.deltaT        =   0.25 ;
analysisSettings.finalTime     =   20    ;
analysisSettings.stopTolDeltau =   0    ;
analysisSettings.stopTolForces =   1e-7 ;
analysisSettings.stopTolIts    =   30   ;


analysisSettings.methodName    = 'alphaHHT' ;
%md
### otherParams
%mdA name problem and vtk output is set:
otherParams.problemName = 'rightAngleCantilever'; 
otherParams.plotsFormat = 'vtk' ;
%mdONSAS code is run for the input structs stated above
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md the control dof to validate the solution is $u_z$ and $u_y$ at the loaded node:
dofDispY = (nElemsPerBeam + 1)*6 - 3 ;
dofDispZ = (nElemsPerBeam + 1)*6 - 1 ;
controlDispUy =  matUs(dofDispY,:)  ;
controlDispUz =  matUs(dofDispZ,:)  ;
loadFactorsNREngRot  =  loadFactorsMat(:,2) ;
timeVec = linspace(0, analysisSettings.finalTime, analysisSettings.finalTime / analysisSettings.deltaT + 1 )' ;
%md### Plots
%mdOutput diplacments and load factor function are plotted:
lw = 2.0 ; lw2 = 1.0 ; ms = 11 ; plotfontsize = 18 ;
figure
plot( timeVec, loadFactorsNREngRot ,'b-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time (t)');   laby = ylabel('\lambda (t)') ;
% set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
% set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
% print('output/rightAngleCantileverLoadFactors.png','-dpng')

figure
plot( timeVec, controlDispUy ,'r-x' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time (t)');   laby = ylabel('uA_y ') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
% print('output/rightAngleCantilever_uA_y.png','-dpng')

figure
plot( timeVec, controlDispUz ,'k-x' , 'linewidth', lw,'markersize',ms )
labx = xlabel('time (t)');   laby = ylabel('uA_z ') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
% print('output/rightAngleCantilever_uA_z.png','-dpng')