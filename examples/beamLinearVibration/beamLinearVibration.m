%md# Linear Dynamic Vibration of a Simply Supported Beam
%md
%md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS.m/blob/master/examples/beamLinearVibration/beamLinearVibration.m)
%md
%mdIn this tutorial, the dynamic response of a simply supported beam is computed using ONSAS with the linear elastic and corotational formulations. The aim of this example is to validate the numerical implementations using the analytic solution.
%md
%mdThe problem consists in a beam with fixed nodes in both ends. In a selected position of the beam a forced load with time dependency $(F(t) = f_o \sin(wt))$ is applied in both perpendicular directions, as it is shown in the figure, where $f_o$ is the force magnitude amplitude and $\omega$ is the force angular frequency.
%md
%md```@raw html
%md<img src="../../assets/beamDynamicVibration.svg" alt="structure diagram" width="500"/>
%md```
%md
%mdBefore defining the structs, the workspace is cleaned and the ONSAS directory is added to the path
close all, clear all, addpath( genpath( [ pwd '/../../src'] ) );
%md
%mdExternal load parameters, time values, material and geometric parameters are defined to find both analytic and numerical solutions. The material scalar parameters are
E = 200e9 ; nu = 0.3;  rho = 700;
%mdand the geometrical scalar parameters are
l = 10 ; ty = .3 ;  tz = .1 ;
Iyy = ty*tz^3/12 ;
Izz = tz*ty^3/12 ;
numElements = 10 ; % Number of elements

%md Time and applied force parameters are
Fo     = 100   ; % N
w      = 2     ; % rad/s
tf     = 8     ; % s
deltat = 0.1 ; % s
%md The node of application of the external forces is assumed to be located at the mid-point of the beam  
assert( rem( numElements, 2 ) == 0, 'the number of elements must be even.' )
appNode    = ( numElements ) / 2 + 1       ;
appNodePos = (appNode-1) * l / numElements ;
%md
%md## Analytic solution
%md
%md The governing equations for a beam with uniform cross-section, density and Young modulus with a transversal distributed applied load $q$ is given by
%md```math
%md EI \frac{\partial^4 w}{\partial x^4}(x,t) + \rho A \frac{\partial^2w}{\partial t^2}(x,t) = q(x,t)
%md```
%mdconsidering the time depdency load, $F(x, t)=fo sin(wt)$, and using a Fourier decomposition with defined initial and boundarie conditions following the mathematical process explained on chapter 10 of _Mechanical Vibrations_ (5th Edition; Rao Singiresu) we obtain the analytic solution for the vertical displacement of our problem 
%md```math
%md w(x,t) = \frac{2fo}{\rho A l} \sum_{n=1}^{\infty} \frac{1}{w_{n}^2 - w^2} \sin\left(\frac{n \pi a}{l} \right) \sin\left(\frac{n \pi x}{l} \right)\sin(wt)
%md```
%md
%mdThe solution can be numerically computed setting a mesh of spatial poins and a vector of times
%md
ts = 0:deltat:tf       ; % times vector
xs = 0:l/numElements:l ; % beam mesh
ns = 1:10              ; % modes
%md
%md As well as the vector of natural frequencies
%md
wnY = ( (ns*pi).^2 ) * sqrt(E*Izz/rho/(ty*tz)/(l^4)) ; % Natural frecuency direction Y
wnZ = ( (ns*pi).^2 ) * sqrt(E*Iyy/rho/(ty*tz)/(l^4)) ; % Natural frecuency direction Z
%md
%md The values of the analytic solution are computed as
%md
analyticDisY = 0; analyticDisZ = 0;
analySolPos = appNodePos ;  % the analytic solution is computed at the point of application of the load
for i = 1:length( ns )
  analyticDisY = analyticDisY + (1/(wnY(i)^2 - w^2)) * sin(i*analySolPos/l*pi) * sin(i*pi*appNodePos/l) * sin(w*ts)' ;
  analyticDisZ = analyticDisZ + (1/(wnZ(i)^2 - w^2)) * sin(i*analySolPos/l*pi) * sin(i*pi*appNodePos/l) * sin(w*ts)' ;
end
analyticDisY = analyticDisY * (2*Fo/(rho*ty*tz*l) ) ;   analyticDisZ = analyticDisZ * (2*Fo/(rho*ty*tz*l) ) ;
%md
%md## Numerical solution
%md### MEBI parameters
%md
%mdThe modelling of the structure begins with the definition of the Material-Element-BoundaryConditions-InitialConditions (MEBI) parameters.
%md
%md### materials
%md Since the example contains only one rod and no nodal masses are used, only one `materials` struct is defined. The first analysis is done using the co-rotational formulation
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu]          ;
materials.density         = rho              ;
%md
%md### elements
%md
%mdTwo different types of elements are considered, `node` and `beam`. The nodes will be assigned in the first entry (index $1$) and the beam at index $2$. The elemType field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md for the crossSection, a frame element of rectangular-cross section with $t_y$ and $t_z$ dimensions in $y$ and $z$ directions is set, then the elemTypeGeometry field is:
elements(2).elemCrossSecParams = { 'rectangle' , [ty tz] } ;
%md The consistent mass approach is considered for the dynamic analysis
elements(2).massMatType = 'consistent';
%md
%md### boundaryConditions
%md
%md The elements are submitted to two different BC settings. The first BC corresponds to the fixed points
boundaryConds(1).imposDispDofs = [ 1 2 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
%md and the second corresponds to a time dependant external force
boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsTimeFact = @(t) Fo*sin(w*t) ;
boundaryConds(2).loadsBaseVals = [ 0 0 1 0 1 0 ] ;
%md
%md### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
%md### mesh parameters
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
%mdThe connectivity is introduced using the _conecCell_ cell. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of the nodes of the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first two nodes are defined, both with material zero (since nodes dont have material), the first element type (the first entry of the cells of the _elements_ struct), and the first entry of the cells of the boundary conditions struct. No non-homogeneous initial condition is considered (then zero is used) and finally the node index is included.
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
%md the following case only differs in the boundary condition and the node number
mesh.conecCell{ 2, 1 } = [ 0 1 1 0  numElements+1 ] ;
%md the following case only differs in the boundary condition and the node number
mesh.conecCell{ 3, 1 } = [ 0 1 2 0  appNode ] ;
%md the beam elements are formed by the first material, the second type of element, and no boundary conditions are applied to any element.
for i=1:numElements
  mesh.conecCell{ i+3,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   deltat  ;
analysisSettings.finalTime     =   tf   ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   10   ;
%md
%md## otherParams
otherParams.problemName = 'coRotationaluniformDynamicBeam';
otherParams.plotsFormat = 'vtk' ;
%md ONSAS execution
[coRotMatUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md The second analysis case implements the linear elastic formulation
materials.hyperElasModel  = 'linearElastic' ;
otherParams.problemName = 'linearElasticuniformDynamicBeam';
%md
%md ONSAS execution
[linElasMatUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md### Error estimation
dofYendNode = 6*(appNode) - 3;     dofZendNode = 6*(appNode) - 1;
%md
%md error computed for each method in the application node of the external force
diflinearDispUy = linElasMatUs(dofYendNode, :)' - analyticDisY ;
diflinearDispUz = linElasMatUs(dofZendNode, :)' - analyticDisZ ;
difcoRotDispUy  = coRotMatUs(dofYendNode, :)'   - analyticDisY ;
difcoRotDispUz  = coRotMatUs(dofZendNode, :)'   - analyticDisZ ;
%md
errlinearDispUy = norm( diflinearDispUy, 1 ) / norm( analyticDisY, 1 ) ;
errcoRotDispUy  = norm( difcoRotDispUy , 1 ) / norm( analyticDisY, 1 ) ;
errlinearDispUz = norm( diflinearDispUz, 1 ) / norm( analyticDisZ, 1 ) ;
errcoRotDispUz  = norm( difcoRotDispUz , 1 ) / norm( analyticDisZ, 1 ) ;
%md
%md the numerical resolution is validated for both method and both directions.
verifBoolean =  ( errlinearDispUy <  5e-2 ) && ( errcoRotDispUy  <  5e-2 ) ...
             && ( errlinearDispUz <  5e-2 ) && ( errcoRotDispUz  <  5e-2 );
%md
%md Plot parameters:
lw = 2.0 ; lw2 = 1.0 ; ms = 11 ; plotfontsize = 18 ;
%md plot y-axis linear, co-rotational and analytic result 
figure, hold on, grid on
plot(ts, coRotMatUs(dofYendNode, :),'r-x' , 'linewidth', lw,'markersize',ms )
plot(ts, linElasMatUs(dofYendNode, :),'k-o' , 'linewidth', lw,'markersize',ms )
plot(ts, analyticDisY,'b' , 'linewidth', lw,'markersize',ms )
legend('coRotational_{disp}','linearElastic_{disp}', 'Analytic_{disp}', 'location', 'eastoutside')
labx = xlabel('time (s)');   laby = ylabel('displacement (m)') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/Uy','-dpng')
%md plot z-axis linear, co-rotational and analytic result 
figure, hold on, grid on
plot(ts, coRotMatUs(dofZendNode, :),'r-x' , 'linewidth', lw, 'markersize', ms )
plot(ts, linElasMatUs(dofZendNode, :),'k-o' , 'linewidth', lw, 'markersize', ms )
plot(ts, analyticDisZ, 'b' , 'linewidth', lw, 'markersize', ms )
legend('coRotational_{disp}', 'linearElastic_{disp}', 'Analytic_{disp}', 'location', 'eastoutside')
labx = xlabel('time (s)');   laby = ylabel('displacement (m)') ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('output/Uy.png','-dpng')
%print('../../docs/src/assets/beamDynamicVibrationVerifUy.png','-dpng')
%md
%md```@raw html
%md<img src="../../assets/beamDynamicVibrationVerifUy.png" alt="structure diagram" width="500"/>
%md```
