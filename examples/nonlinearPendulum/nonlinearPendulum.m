%md# Non linear pendulum
%md---
%md
%mdIn this tutorial, the nonlinear pendulum example and its resolution using ONSAS are described. The aim of this example is to validate the implementations of the $\alpha$-HHT and Newmark methods. The problem consists in a truss element with distributed mass submitted to its own weight.
%md
%mdThe examplpe is extracted from [(K.J Bathe 2006)](https://books.google.com.uy/books?hl=es&lr=&id=rWvefGICfO8C&oi=fnd&pg=PR13&dq=Bathe+finite+element+procedures+book&ots=gHFGuStsC0&sig=odG4LfWjFJ3CdkQxhnn1DVZcUzs#v=onepage&q=Bathe%20finite%20element%20procedures%20book&f=false)  
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/nonLinearPendulum_.svg" alt="structure diagram" width="500"/>
%md```
%mdBefore defining the input structs all variables are cleaned, the open windows are closed and the source folder added to the workspace path:
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
% material scalar parameters
%mdThe material and geometrical properties must comply certain equals:
%md```math
%md EA = 10^8, A=0.1, l0 = 3.0443
%md```
%mdand the period solution is
%md```math
%md T = 4.13 s,
%md```
%mdtherefore these magnitudes are included in the code by:
EA = 1e8     ; A = 0.1  ;
E   = EA / A ; nu = 0   ;
l0  = 3.0443 ; T = 4.13 ;
m   = 10     ; g = 9.80 ;
%mdIn the example presented in [(K.J Bathe 2006)] the force $m*g$ is applied at the end node, consequently a syntheticall coherent densitiy to generate the same effect is:
rho = 2*m / ( A * l0 )  ;
%md##Numerical solution
%md### MEBI parameters
%md
%mdThe modelling of the structure begins with the definition of the material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md### materials
%mdSince the first case contains only one type of material the fields of the `materials` struct will have only one entry.
materials.density = rho ;
%md Moreover, the constitutive behavior considered is the Rotated Engineering strain, thus the field `hyperElasModel` is:
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, node and truss. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The `elemType` field is then:
elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
%mdA rectangular $2$ section is considered with sqrt(A)xsqrt(A). However this type of section has no effect in the results, because of the inertial primacy against stiffness terms. Subsequently `elemCrossSecParams` field is:
elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [ sqrt(A) sqrt(A) ] ;
%mdand the according to the literature example the element include conssitent mass matrix
elements(2).massMatType = 'lumped' ;
%md
%md### boundaryConds
%md
%mdThe elements are submitted to two different BC settings. The first BC corresponds to a simple fixed condition (all 3 dipslacments dofs of the node are equal to zero during all time span) then the `boundaryConds(1)` set is:
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
%mdthen, since for the first case no boolean self weight is considered, the weight force load is applied at the end node adding a second boundary condition `boundaryConds(2)`:
boundaryConds(2).imposDispDofs =  3 ;
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) 1.0 ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -m*g 0 ] ;
%md
%md### initial Conditions
%mdNon homogeneous initial conditions are not considered in this example, consequently the `initialConds`  struct is set empty:
initialConds                = struct() ;
%md
%md### mesh parameters
%mdThe coordinates conisdering a mesh of two nodes is:
mesh.nodesCoords = [   0  0   l0 ; ...
                      l0  0  l0  ] ;
%mdThe connectivity is introduced using the _conecCell_ cell. Each entry of the cell (indexed using {}) contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md Then the entry of node $1$ is introduced:
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%md the first MEBI parameter (Material) is set as _zero_ (since nodes dont have material). The second parameter corresponds to the Element, and a _1_ is set since `node` is the first entry of the  `elements.elemType` cell. For the BC index, we consider that node $1$ is simple fixed, then the first index of the `boundaryConds` struct is used. Finally, no specific initial conditions are set for the node (0) and at the end of the vector the number of the node is included (1).
%md A similar approach is used for node $3$,
mesh.conecCell{ 2, 1 } = [ 0 1 2 0  2   ] ;
%md and for node $2$ only the boundary condition is changed, because it is lodaded.
%md Regarding the truss elements, the first material is considered, the second type of element, and no boundary conditions are applied.
mesh.conecCell{ 3, 1 } = [ 1 2 0 0  1 2 ] ;
%md
%md### analysisSettings
%mdA Newmark algorithm is used to solve this problem with the following parameters during one period $T$:
analysisSettings.deltaT        = 0.05  ;
analysisSettings.finalTime     = 1* T  ;
analysisSettings.stopTolDeltau = 1e-12 ;
analysisSettings.stopTolForces = 1e-12 ;
analysisSettings.stopTolIts    = 30    ;
otherParams.plots_format       = 'vtk' ;

%md### Analysis case 1: Solution using Newmark with truss element and mass lumped and the weight force is included by external force according to Bathe problem
analysisSettings.methodName = 'newmark'     ;
otherParams.problemName     = 'nonlinearPendulumNewmarkTrussBathe';
%md In the first case ONSAS is run and the solution at the dof of interest is stored.
[matUsCase1, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
% ------------------------------------
%md### Analysis case 2: Solution using HHT with truss element, mass lumped according to Bathe problem and self weight boolean activated:
analysisSettings.booleanSelfWeight = true ;
%mdand the external weight load is set to zero
boundaryConds(2).loadsTimeFact = @(t) 0.0 ;
%mdIn order to validate HHT numerical method this is executed with $alphaHHT = 0$ which necessary implies that numerical results must be identical to CASE 1, since HHT is equivalent to Newmark if AplhaHHT = 0,
analysisSettings.methodName = 'alphaHHT';
analysisSettings.alphaHHT   =  0        ;        
%mdthe new name of these problem is:
otherParams.problemName     = 'nonlinearPendulumHHTTrussBathe';
[matUsCase2, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
% ------------------------------------
%md### Analysis case 3: Solution using HHT using frame element, mass lumped at the final element and self weight boolean is activated. For this case denisity is null for the rest of the elements and rhoF =  m / ( A * (lNEWelem * l0) ) for the new one in order to produce the same force that is considered by [(K.J Bathe 2006)]
%md
%mdNow the element type number two is a frame:
elements(2).elemType     = 'frame';
%mdand the according to the literature example the element include conssitent mass matrix
elements(2).massMatType  = 'consistent';
%and the fraction of the new element in the pendulum length:
lumpedParam = 0.01 ;
mesh.nodesCoords = [   0                    0    l0 ; ...
                       (1-lumpedParam)*l0   0    l0 ; ... 
                       l0                   0    l0 ] ;
%mdA conec cell row is added for the last element whith the second material, second element, boundary condition and no initial condition:
mesh.conecCell{ 4, 1 } = [ 2 2 0 0  2 3 ] ;
%mdA new torsional constraint must be added beacuse the torsional spin of the frame element: 
boundaryConds(1).imposDispDofs = [ 1 2 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
%mdIn addition no weight must be computed by the body of the pendulum so the material with null density in the body of the frame:
materials(1).density = 0 ;
%mdalso the new material that produce m*g at the end of the frame is:
rhoF = m / ( A * (lumpedParam * l0) ) ;
materials(2).hyperElasModel  = '1DrotEngStrain' ;
materials(2).hyperElasParams = [ E nu ] ;
materials(2).density = rhoF ;
otherParams.problemName     = 'nonlinearPendulumHHTFrame';
% ------------------------------------
[matUsCase3, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

otherParams.problemName     = 'nonlinearPendulumHHTFrameWithSpring';
boundaryConds(1).springDofs = [ 4 ] ;
boundaryConds(1).springVals = [ 1e2 ] ;

[matUsCase4, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

figure
plot(matUsCase3(4,:))
grid on, hold on
plot(matUsCase4(4,:),'r')

%md### extract control displacements
%mdThe mass displacement in z are:
controlDofDispUz = 6 + 5 ;
controlDispZCase1 = matUsCase1( controlDofDispUz , : ) ;
controlDispZCase2 = matUsCase2( controlDofDispUz , : ) ;
controlDofDispZ3 = 12 + 5 ;
controlDispZCase3 = matUsCase3( controlDofDispZ3 , : ) ;
%mdanalogously the mass dipslacement in x are:
controlDofDispUx12 = 6 + 1 ;
controlDispXCase1 = matUsCase1( controlDofDispUx12 , : ) ;
controlDispXCase2 = matUsCase2( controlDofDispUx12 , : ) ;
controlDofDispX3 = 12 + 1 ;
controlDispXCase3 = matUsCase3( controlDofDispX3 , : ) ;

%mdIn order to contrast the solution with the literature refrence the bounce angle measured from the vertical is computed:
angleThetaCase1 = rad2deg( atan2( ( l0 + controlDispXCase1 ), -controlDispZCase1 ) ) ;
angleThetaCase2 = rad2deg( atan2( ( l0 + controlDispXCase2 ), -controlDispZCase2 ) ) ;
angleThetaCase3 = rad2deg( atan2( ( l0 + controlDispXCase3 ), -controlDispZCase3 ) ) ;

%mdTo plot diplsacements against $t$ the time vector is:
timesVec12  = (0:length(controlDispZCase1)-1) * analysisSettings.deltaT ;

%md## verification
%mdFor all cases the displacment at t=T must be close to zero, so then the error is computed as $uN - 0$ / l0:
tolVerifDisp = 1e-2 ;
verifBooleanCase1 =  ( abs( controlDispZCase1(end) / l0 ) <  tolVerifDisp ) ;
verifBooleanCase2 =  ( abs( controlDispZCase2(end) / l0 ) <  tolVerifDisp ) ;
verifBooleanCase3 =  ( abs( controlDispZCase3(end) / l0 ) <  tolVerifDisp ) ;
%md all cases must be verifyed, so then:
verifBoolean    = verifBooleanCase1 && verifBooleanCase2 && verifBooleanCase3;
 
%md### Plots
%md
%mdPlot parameters
MS = 10; LW = 1.5 ;
legendCase1 = [' Truss lumped Bathe Trapezoidal Newmark'];
legendCase2 = [' Truss lumped Bathe HHT with \alpha =' num2str(analysisSettings.alphaHHT)];
legendCase3 = [' Frame HHT with \alpha =' num2str(analysisSettings.alphaHHT)];
%mdPlot displacement solution of bathe problem
figure, hold on, grid on
plot( timesVec12, -controlDispZCase1, 'k-s' ,'markersize', MS,'linewidth', LW)
plot( timesVec12, -controlDispZCase2, 'bo','markersize', MS,'linewidth', LW)
plot( timesVec12, -controlDispZCase3, 'rx','markersize', MS,'linewidth', LW)
xlabel('time (s)'), ylabel('mass displacement u_z (m)')
legend( legendCase1, legendCase2, legendCase3, 'location','NorthEast')
title("U_z solution Bathe")
%print('./output/dispPlot.png','-dpng')

%mdPlot a ngle solution
figure, hold on, grid on
plot( timesVec12, -angleThetaCase1, 'k-s' ,'markersize', MS,'linewidth', LW)
plot( timesVec12, -angleThetaCase2, 'bo','markersize', MS,'linewidth', LW)
plot( timesVec12, -angleThetaCase3, 'rx','markersize', MS,'linewidth', LW)
xlabel('time (s)'), ylabel('\theta displacement (º)')
legend( legendCase1, legendCase2, legendCase3, 'location','NorthEast')
title("Θ solution Bathe")
%print('./output/thetaPlot.png','-dpng')

