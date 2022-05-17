%md# Dynamic Vibration modes of a Beam with fix nodes in both ends.
%md---
%md
%mdIn this tutorial, the dynamic response of a beam with fix nodes in both ends is solved using ONSAS implementing the corotational model and the linear elastic model. The aim of this example is to validate the dynamic response of a beam, to do this the analytic solution is comapred with the ONSAS result and the error is estimated in order to validete the result. The Octave script of this example is available at [this url](https://github.com/ONSAS/ONSAS.m/blob/master/examples/uniformCurvatureCantilever/onsasExample_uniformCurvatureCantilever.m).
%md
%mdThe problem consists in a beam with fixed nodes in both ends. In a selected position of the beam a forced load with time dependency is applied $(F = F_o sin(wt))$, as it is shown in the figure.
%md
%md```@raw html
%md<img src="assets/dynamicBeamHTML.svg" alt="structure diagram" width="500"/>
%md```
%md
%mdBefore defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar geometry, material parameters and load values are defined.
close all, clear all ;
% add path
addpath( genpath( [ pwd '/../../src'] ) );
% External forced load parameters
Fo     = 100; % N
w      = 2;    % rad/s
% Time values
tf     = 8;    % sec
deltat = 0.1; % sec
%md
% Material scalar parameters
E = 200e9 ;  nu = 0.3 ; rho = 700;
%md geometrical scalar parameters
l = 10 ; ty = .3 ;  tz = .1 ;
Iyy = ty*tz^3/12 ;
Izz = tz*ty^3/12 ;
%md Number of elements
numElements = 51 ;
%md
%md## Analytic solution
%md
%md The dynamic displacement of a forced beam describe by the next differential equation
%md```math
%md EI \frac{\partial^4 w}{\partial x^4} + \rho A \frac{\partial^2
%w}{\partial t^2} = f(x,t)
%md```
%md Implementig a solution w(x,t) = W(x)T(t) it is possible to find:
%md```math
%md w(x,t) = \frac{2fo}{\rho A l}\sum_{n=1}^{\infty } \frac{1}{w_{n}^2 - w^2}\sin(\frac{n \pi a}{l})\sin(\frac{n \pi x}{l})\sin(wt)
%md```
%md
%md## External Load application node
appNode     = (numElements+1)/2+1;
appNodePos  = l*(appNode)/numElements;
%md
%md## Analytic solution of a beam with fix nodes in both ends.
%analyticSol  = @(x, t, wn, n) (1/(wn^2 - w^2))*sin(n*pi*a/l)*sin(n*pi*x/l)*sin(w*t);
%md
t  = 0:deltat:tf; % time vector
x  = 0:l/numElements:l; % beam mesh
n  = 1:1:8; % number of nodes
analyticDisV = 0;
analyticDisW = 0;
%md   
wnY = ((n*pi).^2)*sqrt(E*Izz/rho/(ty*tz)/(l^4)); % Natural frecuency direction Y
wnZ = ((n*pi).^2)*sqrt(E*Iyy/rho/(ty*tz)/(l^4)); % Natural frecuency direction Z
%md
for i=1:length(n)
    analyticDisV = analyticDisV + (2*Fo/(rho*ty*tz*l))*sin(i*x*pi/l).*sin(i*pi*appNodePos/l)*(1./(wnY(i)^2 - w^2)).*sin(w.*t)';
    analyticDisW = analyticDisW + (2*Fo/(rho*ty*tz*l))*sin(i*x*pi/l).*sin(i*pi*appNodePos/l)*(1./(wnZ(i)^2 - w^2)).*sin(w.*t)';
end
%md
%md##  Numerical solution
%md### MEBI parameters
%md
%mdThe modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md### Materials
%md Since the example contains only one rod the fields of the `materials` struct will have only one entry. Although, it is considered constitutive behavior according to the SaintVenantKirchhoff law:
%md The first analysis case implements the co-rotational formulation
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu] ;
materials.density = rho;
%md
%md### Elements
%md
%mdTwo different types of elements are considered, node and beam. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the geometries, the node has not geometry to assign (empty array), and the truss elements will be set as a rectangular-cross section with $t_y$ and $t_z$ cross-section dimensions in $y$ and $z$ directions, then the elemTypeGeometry field is:
elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [ty tz]     ;
elements(2).elemTypeParams          = 1           ;
%md the consistent mass type is selected to solve the dynamic response of
%the beam
elements(2).massMatType             = 'consistent';
%md
%md### boundaryConds
%md
%md The elements are submitted to two different BC settings. The first BC
%corresponds to the fixed points
boundaryConds(1).imposDispDofs = [ 1 2 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
%mdand the second corresponds to an external force with time dependency and
%an application frequency $w$
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsTimeFact = @(t) Fo*sin(w*t) ;
boundaryConds(2).loadsBaseVals = [ 0 0 1 0 1 0 ] ;
%md
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds                   = struct() ;
%md
%md### mesh parameters
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%md the following case only differs in the boundary condition and the node number
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  numElements+1 ] ;
%md the following case only differs in the boundary condition and the node number
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  (numElements+1)/2+1 ] ;
%md the beam elements are formed by the first material, the second type of element, and no boundary conditions are applied to any element.
for i=1:numElements,
  mesh.conecCell{ i+3,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   deltat  ;
analysisSettings.finalTime     =   tf   ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;
%md
%md## otherParams
otherParams.problemName = 'coRotationaluniformDynamicBeam';
otherParams.controlDofs = [ numElements/2+1 3 ] ;
otherParams.plotsFormat = 'vtk' ;
%md Beam with simple supported nodes in the ends
[coRotMatUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md The second analysis case implements the linear elastic formulation
materials.hyperElasModel  = 'linearElastic' ;
otherParams.problemName = 'linearElasticuniformDynamicBeam';
otherParams.controlDofs = [ numElements/2+1 3 ] ;
otherParams.plotsFormat = 'vtk' ;
%md
%md## Analysis case 1: Linear Elastic Y
%md Beam with simple supported nodes in the ends
[linElasMatUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md## Analytic solution
%md
%md### Plots
dofYendNode = 6*(appNode) - 3;
dofZendNode = 6*(appNode) - 1;
%mdPlot parameters:
lw = 2.0 ; lw2 = 1.0 ; ms = 11 ; plotfontsize = 18 ;
%md time vector
timeVec = linspace( 0, tf, size(coRotMatUs,2) );
%md
% error estimated for each method in the application point of the
% external time dependency load
diflinearDispUy = linElasMatUs(dofYendNode, :) - analyticDisV(: , appNode);
diflinearDispUz = linElasMatUs(dofZendNode, :) - analyticDisW(: , appNode);
difcoRotDispUy  = coRotMatUs(dofYendNode, :) - analyticDisV(: , appNode);
difcoRotDispUz  = coRotMatUs(dofZendNode, :) - analyticDisW(: , appNode);
%
intAnalyticDisV = 0;
intAnalyticDisW = 0;
errlinearDispUy = 0;
errcoRotDispUy  = 0;
errlinearDispUz = 0;
errcoRotDispUz  = 0;
%
for h = 1:length(diflinearDispUy)-1
    intAnalyticDisV = intAnalyticDisV + (analyticDisV(h+1 , appNode) + analyticDisV(h , appNode))*deltat/2;
    intAnalyticDisW = intAnalyticDisW + (analyticDisW(h+1 , appNode) + analyticDisW(h , appNode))*deltat/2;
    errlinearDispUy = errlinearDispUy + (diflinearDispUy(h+1) + diflinearDispUy(h))*deltat/2;
    errcoRotDispUy  = errcoRotDispUy + (difcoRotDispUy(h+1) + difcoRotDispUy(h))*deltat/2;
    errlinearDispUz = errlinearDispUz + (diflinearDispUz(h+1) + diflinearDispUz(h))*deltat/2;
    errcoRotDispUz  = errcoRotDispUz + (difcoRotDispUz(h+1) + difcoRotDispUz(h))*deltat/2;
end
%md
errlinearDispUy = norm(errlinearDispUy)/norm(intAnalyticDisV);
errcoRotDispUy  = norm(errcoRotDispUy)/norm(intAnalyticDisV);
errlinearDispUz = norm(errlinearDispUz)/norm(intAnalyticDisW);
errcoRotDispUz  = norm(errcoRotDispUz)/norm(intAnalyticDisW);
%the numerical resolution is validated for both method and both directions.
verifBoolean =  ( errlinearDispUy <  1e-1 ) ...
             && ( errcoRotDispUy  <  1e-1 ) ...
             && ( errlinearDispUz <  1e-1 ) ...
             && ( errcoRotDispUz  <  1e-1 );
figure(1), hold on, grid on
%md plot co-rotational solution
plot(timeVec, coRotMatUs(dofYendNode, :),'r-x' , 'linewidth', lw,'markersize',ms )
%md plot linear elastic solution
plot(timeVec, linElasMatUs(dofYendNode, :),'k-o' , 'linewidth', lw,'markersize',ms )
%md plot analytic solution
plot(timeVec, analyticDisV(:,appNode),'b' , 'linewidth', lw,'markersize',ms )
% legends
legend('coRotational_{disp}','linearElastic_{disp}', 'Analytic_{disp}', 'time displacement 0 m',...
       'location', 'eastoutside')
labx = xlabel('time (s)');   laby = ylabel('displacement (m)') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/Uy','-dpng')
%md
figure(2), hold on, grid on
% plot linear cases
plot(timeVec, coRotMatUs(dofZendNode, :),'r-x' , 'linewidth', lw, 'markersize', ms )
% plot co-rotational cases
plot(timeVec, linElasMatUs(dofZendNode, :),'k-o' , 'linewidth', lw, 'markersize', ms )
% plot analytic cases
plot(timeVec, analyticDisW(:, appNode), 'b' , 'linewidth', lw, 'markersize', ms )
% legends
legend('coRotational_{disp}', 'linearElastic_{disp}', 'Analytic_{disp}', 'time displacement 0 m',...
       'location', 'eastoutside')
labx = xlabel('time (s)');   laby = ylabel('displacement (m)') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/Uy','-dpng')
