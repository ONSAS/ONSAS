%md# Aerodynamic non-linear cantilever beam example
%md
%md [![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS.m/blob/master/examples/nonLinearCantileverAero/nonLinearCantileverAero.m)
%md
%md In this tutorial, the aerodynamic non-linear cantilever beam example is solved using ONSAS. The aim of this problem is to validate aerodynamic steady wind loads applied to a cantilever beam undergoing small strains. The aerodynamic force modification due to the beam deformation is considered (drag reconfiguration). Given the aforementioned characteristics and under the hypothesis of small displacements regime a semi-analytic solution is available.   
%md
%md The beam is submitted to a uniform air wind velocity field $V_a$, at 20 degrees and atmospheric pressure, along axis $y$. Due to revolution symmetry of the problem lift $c_l$ and torsional moment $c_m$ coefficients are null. A drag coefficient $c_d=1.2$ is extracted from [this reference](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/experiments-on-the-flow-past-a-circular-cylinder-at-very-high-reynolds-number/7859A6C46BF4B0F43F11F52AE1C60150). 
%mdThe beam has a length $L$ and a circular solid cross section with diameter $d$ as it is shown in the following figure: 
%md
%md```@raw html
%md<img src="../../assets/nonLinearAerodynamics/ilusNonLinearAerodynamics.svg" alt="general dimensions sketch" width="700"/>
%md```
%md
%md## Small displacements 2D case
%md--------------------
%md
%md### Static semi-analytic solution
%md
%md The wind load forces of a generic cross section can be derived within the quasi-steady-theory. Considering a cross section located at $x$, then the projected wind velocity into the transverse deformed plane is then  $V_p=V_acos(\theta_z)$ (the axial drag is neglected). Subsequently a drag force per unit of length $F_d= \frac{1}{2} \rho d c_d ||V_p||^2$ with $\frac{V_p}{||V_p||}$ direction is applied. In order to link the force $F_d$ with the beam deflection, the uniform distributed force along $y$ is computed as $F_y=F_d.cos(\theta_z)$. This leads to the following third order differential equation:
%md
%md * $EI_{zz} \frac{\partial ^3 \theta_z}{\partial x ^3} = -q_0 c_d\cos^3(\theta_z)$
%md in which $q_0 = \frac{1}{2} \rho_f d ||V_a||^2$ and the air density is $\rho_f = 1.225$ kg/m$^3$.
%md### Numerical solution
%md---------------------
%md
%md Before defining the structs, the workspace is cleaned and ONSAS directory is added:
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );
%md material and geometrical parameters are defined:
E = 1e9 ;  nu = 0.3 ; rho = 1800 ; G = E / (2 * (1+nu)) ;
l = 10 ; d = l/100 ; J = pi * d ^ 4 / 32 ; Iyy = J / 2 ; Izz = Iyy ;  
%md the fluid properties are:
rhoA = 1.225 ; nuA = 1.6e-5;
%md next the number of frame elements for all cases is set:
numElements = 10 ;
%md### MEBI parameters
%md
%md### materials
%md Since the example contains only one linear Euler Bernoulli element the fields of the `materials` struct will have only one entry. Although, the constitutive behavior law selected is Saint-Venant-Kirchhoff:
materials.hyperElasModel  = 'linearElastic' ;
materials.hyperElasParams = [ E nu ]        ;
%md note that the use of  this linear elastic element guarantees the left hand side of the differential equation stated above.
%md
%md### elements
%md
%md Two different types of elements are considered, node and frames. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The _elemType_ field is then:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'frame' ;
%md The node has not cross section geometry to assign (an empty array is automatically set). The solid circular cross section is preset in ONSAS, and to load it just use:
elements(2).elemCrossSecParams{1,1} = 'circle' ;
elements(2).elemCrossSecParams{2,1} = d        ;
%md Now the element aerodynamic properties are defined. First the drag coefficient function located at the folder's example is declared into _userDragCoef_ field as:
elements(2).aeroCoefs   = {'dragCoefCircular'; []; [] }   ;
%md in which the second and third components of the vectors are considered empty since no lift and torsional moment is considered.
%md Next the _elemTypeAero_ field contains the information of the chord vector. This vector is defined first considering the orientation of the cross section set up for drag experiments. According to the revolution symmetry of the problem the chord vector orientation has no impact into drag force vector, since $c_d$ is constant any angle of incidence. However the characteristic dimension of the circular cross section is declared into the norm of the chord vector ( first three entries of _elemTypeAero_ field into `elements` struct ) as: 
numGaussPoints           = 4 ;
elements(2).elemTypeAero = [0 d 0 numGaussPoints true ] ;
%md also 4 number of integration Gauss points are employed to compute each element aerodynamic force vector.
%md
%md### boundaryConds
%md
%md Only one welded (6 degrees of freedom are set to zero) boundary condition (BC) is considered:
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%md
%md### initial Conditions
%md Any non-homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct() ;
%md
%md### mesh parameters
%mdThe coordinates of the mesh nodes are given by the matrix:
mesh.nodesCoords = [ (0:(numElements))'*l/numElements  zeros(numElements+1,2) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
%md Next the frame elements MEBI parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC and no IC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end
%md
%md### analysisSettings
%md
%md First the wind velocity and the fluid properties are set into __ field of `analysisSettings` struct. This will apply a external wind loads for each element with _elemTypeAero_ field into the `elements` struct. The name of the wind velocity function must located on the same example path is: 
analysisSettings.fluidProps = {rhoA; nuA; 'windVelNonLinearStatic'} ;
%md Inside that function a linear velocity $v_a = 30*t$ is declared. The final time will be set to $1$ in order to achieve 30 m/s.
%md The geometrical non-linear effects are considered in this case to compute the aerodynamic force. As consequence the wind load forces are computed on the deformed configuration. The field  _geometricNonLinearAero_ into  `analysisSettings` struct is then set to:
analysisSettings.geometricNonLinearAero = true;
%md  note that if this boolean is not declared ONSAS will automatically assign it as true. 
%md since this problem is static, then a N-R method is employed. The convergence of the method is accomplish with ten equal load steps. The time variable for static cases is a load factor parameter that must be configured into the `windVel.m` function. A linear profile is considered for ten equal velocity load steps as:
analysisSettings.deltaT        =   0.1           ;
analysisSettings.finalTime     =   1             ;
analysisSettings.methodName    = 'newtonRaphson' ;
%md Next the maximum number of iterations per load(time) step, the residual force and the displacements tolerances are set to: 
analysisSettings.stopTolDeltau =   1e-6          ;
analysisSettings.stopTolForces =   1e-6          ;
analysisSettings.stopTolIts    =   40            ;
%md
%md### otherParams
%md The name of the problem and vtk format output are selected: 
otherParams.problemName = 'nonLinearCantileverSD2D';
otherParams.plots_format = 'vtk' ;
%md
%md ONSAS software is executed for the parameters above defined and the displacement solution of each load(time) step is saved into matUsSD matrix:
[matUsSD, ~] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md 
%md### Verification
%md---------------------
%md
%md### DifferentialEquations.jl (reconfiguration) solution.
%md
%md DiffEq.jl solves the third order ordinary differential equation for this case by executing [DiffEq.jl script](https://github.com/ONSAS/ONSAS.m/blob/master/examples/nonLinearCantileverAero/DiffEq.jl). Then  [`assembleJuliaSol.m` script](https://github.com/ONSAS/ONSAS.m/blob/master/examples/nonLinearCantileverAero/assembleJuliaSol.m) function is executed to build the julia solution with `mesh` and `elements` struct as:
% [dSolJulia] = assembleJuliaSol(elements,mesh) ;
%md Then the the relevance linear and angular displacements are extracted using:
% ydefJulia = dSolJulia(3:6:end)              ;
% thetaZdefJulia = dSolJulia(6:6:end)         ;
% xdefJulia = linspace(0,l,length(ydefJulia)) ;
%md
%md### Numeric solution
%md
%md The numerical solution is computed:
xref    = mesh.nodesCoords(:,1)       ;
yref    = mesh.nodesCoords(:,2)       ;
zref    = mesh.nodesCoords(:,3)       ;
xdefNum = xref + matUsSD(1:6:end,end) ;
ydefNum = yref + matUsSD(3:6:end,end) ;
thetaZdefNum = matUsSD(6:6:end,end)   ;
%md
%md### Plot verification
%md
%md The plot parameters are:
lw = 2 ; ms = 12 ;
labelTitle = [' Validating solution with ' num2str(numElements) ' elements' ] ;
axislw = 2 ; axisFontSize = 20 ; legendFontSize = 15 ; curveFontSize = 15 ;    
folderPathFigs = './output/figs/' ;
mkdir(folderPathFigs) ;
%md The linear $u_y$ displacements verification is plotted executing:  
fig1 = figure(1) ;
hold on, grid on
plot(xref      , ydefNum  , 'bo' , 'linewidth', lw, 'markersize', ms+5   );
% plot(xdefJulia , ydefJulia, 'b-' , 'linewidth', lw, 'markersize', ms     );
legend('y numeric SD', 'y semi-analytic SD' )
labx=xlabel(' x[m] ');    laby=ylabel('y[m]');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location','northWest' ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig1 = strcat(folderPathFigs, 'linDispSD.png') ;
print(fig1, namefig1,'-dpng') ;
%md
%md```@raw html
%md<img src="../../assets/nonLinearAerodynamics/linDispSD.png" alt="plot check angular displacements" width="500"/>
%md```
%md
%md The angular $\theta_z$ displacements verification is plotted executing:  
fig2 = figure(2) ;
hold on, grid on
plot(xref,      thetaZdefNum,      'bo' , 'linewidth', lw,'markersize', ms+5   );
% plot(xdefJulia, thetaZdefJulia,    'b-' , 'linewidth', lw, 'markersize', ms    );
legend('\theta_z numeric SD', '\theta_z semi-analyitc SD')
labx=xlabel(' x[m] ');    laby=ylabel('Angle[rad]');
% title (labelTitle)
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location','northWest' ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig2 = strcat(folderPathFigs, 'angDispSD.png') ;
print(fig2, namefig2,'-dpng')
%md
%md```@raw html
%md<img src="../../assets/nonLinearAerodynamics/angDispSD.png" alt="plot check angular displacements" width="500"/>
%md```
%md
%md## Large displacements 2D case
%md--------------------
%md Now a large displacements 2D case is solved. The solution is computed using the co-rotational beam element formulation proposed in [this reference](https://www.sciencedirect.com/science/article/abs/pii/S0045782513003022)
%md### Numerical solution static case
%md---------------------
%md### MEBI parameters
%md---------------------
%md
%md### materials
%md In order to reproduce large displacements results the `materials` struct is then changed to: 
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ 1e6 nu ]       ;
materials.density         = rho              ;
%md
%md### elements
%md The element tangent matrices of the consistent inertial force vector are taking into account by the following boolean:
elements(2).massMatType = 'consistent' ;
%md
%md### otherParams
%md The name of this case problem is:
otherParams.problemName = 'nonLinearCantileverLDStatic' ;
%md
%md ONSAS software is executed for the parameters above defined and the displacement solution of each load(time) step is saved into matUsLDStatic matrix:
[matUsLDStatic, ~] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md--------------------
%md### Numerical solution dynamic case
%md---------------------
%md Next a dynamic example considering large displacements motion is addressed to test the convergence of the dynamic solution disregarding any artificial damping. 
%md
%md### analysisSettings  
%md For such propose the wind velocity function name is now: 
analysisSettings.fluidProps = {rhoA; nuA; 'windVelNonLinearDynamic2D'} ;
%md Inside that function a ramp velocity profile $v_a(t) = 5*t*(t<6.6) + 5*t*(t>=6.6)$ is declared. This is an abrupt wind velocity load from 0 to $7$ m/s in $10$ s .
%md
%md Regarding the integration time method scheme, a classic Newmark trapezoidal is set as:  
analysisSettings.deltaT     =  1        ;
analysisSettings.finalTime  =  200      ;
analysisSettings.methodName = 'newmark' ;
%md
%md### otherParams
%md The name of this case problem is:
otherParams.problemName = 'nonLinearCantileverLDDynamic2D' ;
%md
%md ONSAS software is executed for the parameters above defined and the displacement solution for each time step is saved into matUsLDDynamic matrix:
[matUsLDDynamic, ~] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md 
%md### Verification
%md---------------------
%md The numerical dynamic convergence to the static solution is then verified . The degree of freedom selected for such task is $u_y(t)$ of node A. 
%md
%md### Static solution.
%md
%md Extract static numerical time history displacements $u_y$ of node A. First the selected degree of freedom is:
nodeA = numElements + 1                ;
relativeDofUyA = 3                     ;
dofUyA = (nodeA -1)*6 + relativeDofUyA ; 
%md then node A $u_y$ time history accessed by:
UyAStaticSol = matUsLDStatic(dofUyA,:) ;
%md
%md### Dynamic solution.
%md
%md Extract dynamic numerical solution as follows:
UyADynamicSol = matUsLDDynamic(dofUyA,:) ;
%md next, the time vector is given by:
timVecLD = linspace(0, analysisSettings.finalTime, size(matUsLDDynamic,2) ) ;
%md
%md### Verification Plot
%md
%md Create folder to save figures
folderFigs = strcat('./output/', 'figs/') ;
mkdir(folderFigs) ;
%md The linear $u_y$ displacements verification of node  A is finally plotted executing:  
fig3 = figure(3) ;
hold on,  grid on
% legend first point plot
plot(timVecLD(1), UyADynamicSol(1),...
     'color', 'b', 'linewidth', lw, 'linestyle', '-','markersize', ms, 'marker', 'o')
% static solution plot
plot(timVecLD   , UyAStaticSol(end)*ones(length(timVecLD)), 'k:' , 'linewidth', lw, 'markersize', ms     );
% markers plot
plot(timVecLD(1:8:end), UyADynamicSol(1:8:end),...
     'color', 'b', 'linewidth', lw, 'linestyle', 'none','markersize', ms, 'marker', 'o')
% continium line plot
plot(timVecLD, UyADynamicSol,...
     'color', 'b', 'linewidth', lw, 'linestyle', '-', 'marker', 'none')
legend('dynamic LD', 'static LD' )
labx=xlabel('x[m]');    laby=ylabel('y[m]');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location','northEast' ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig3 = strcat(folderPathFigs, 'uyA.png') ;
print(fig3, namefig3,'-dpng') ;
%md
%md```@raw html
%md<img src="../../assets/nonLinearAerodynamics/uyA.png" alt="plot check angular displacements" width="500"/>
%md```
%md
%md## Large displacements 3D case
%md--------------------
%md A large displacements dynamic 3D case is presented as follows. This example is inspired on Vortex shedding  exposed at [![Youtbue Video](https://img.shields.io/badge/script-url-blue)](https://www.youtube.com/watch?v=Lf9Ffj5rGh8&ab_channel=FrederickGosselin)
%md### MEBI parameters
%md---------------------
%md
%md### materials  
%md In order to reproduce large displacements results the `materials` struct is then changed to: 
materials.hyperElasParams = [ 1e8 nu ]       ;
%md
%md### analysisSettings  
%md Regarding the integration time method scheme, a classic $\alpha-HHT$ method is employed. This method is more stable numerically than Newmark, the keen reader is refereed to [this reference](https://onlinelibrary.wiley.com/doi/abs/10.1002/eqe.4290050306):
analysisSettings.methodName = 'alphaHHT';
%md the simulation time is defined such that:  
analysisSettings.deltaT     =  .2  ;
analysisSettings.finalTime  =  120 ;
%md The emulation of the vortex shedding vibration is generated by a synthetic wind velocity composed by two sinusoidal velocities. A low frequency $Vy_a$ along the mean flow direction $y$ and then a high frequency component  $Vz_a$ along $z$. The high frequency component is selected to produce resonance effects between the flow and the beam, thus the high frequency velocity is selected equal to the first mode bending:
freqBendingFirstMode = (1.875)^2 * sqrt( materials.hyperElasParams(1) * Iyy / (materials.density * (pi * d^4 / 4) * l^4) ) ;
analysisSettings.fluidProps = {rhoA; nuA; 'windVelNonLinearDynamic3D'} ;
%md The velocity function componentes are assembled:
timeVecLD3d = linspace(0,analysisSettings.finalTime, ceil(analysisSettings.finalTime / analysisSettings.deltaT + 1) ) ;
windVelY = [] ; windVelZ = [] ;
for timeIndex = timeVecLD3d
    windVelVecTimeIndex = feval(analysisSettings.fluidProps{3}, 0, timeIndex) ;
    windVelY = [windVelY windVelVecTimeIndex(2) ] ;
    windVelZ = [windVelZ windVelVecTimeIndex(3) ] ;
end
%md
%md### otherParams
%md The name of this case problem is:
otherParams.problemName = 'nonLinearCantileverLD3D' ;
%md
%md ONSAS software is executed for the parameters above defined and the displacement solution for each time step is saved into matUsLD3D matrix:
[matUsLD3D, ~] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md### Plots
%md
%md The wind velocity profile implemented is plotted executing:  
fig4 = figure(4) ;
hold on,  grid on
% legend first point plot
plot(timeVecLD3d(1), windVelY(1),...
     'color', 'b', 'linewidth', lw, 'linestyle', '-','markersize', ms, 'marker', 'o')
plot(timeVecLD3d(1), windVelZ(1),...
     'color', 'r', 'linewidth', lw, 'linestyle', '-','markersize', ms, 'marker', '^')
% markers plot
plot(timeVecLD3d(1:10:end), windVelY(1:10:end),...
     'color', 'b', 'linewidth', lw, 'linestyle', 'none','markersize', ms, 'marker', 'o')
plot(timeVecLD3d(1:17:end), windVelZ(1:17:end),...
     'color', 'r', 'linewidth', lw, 'linestyle', 'none','markersize', ms, 'marker', '^')
% continium line plot
plot(timeVecLD3d, windVelY,...
     'color', 'b', 'linewidth', lw, 'linestyle', '-', 'marker', 'none')
plot(timeVecLD3d, windVelZ,...
     'color', 'r', 'linewidth', lw, 'linestyle', '-', 'marker', 'none')
legend('Va_y', 'Va_z' )
labx=xlabel('t[s]');    laby=ylabel('V_a[m/s]');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location','northEast' ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig4 = strcat(folderPathFigs, 'windVel3D.png') ;
axis([0,50])
print(fig4, namefig4,'-dpng') ;
%md
%md```@raw html
%md<img src="../../assets/nonLinearAerodynamics/windVel3D.png" alt="plot check angular displacements" width="500"/>
%md```
%md
%md Then $u_y$ of node A  is computed using:  
UyADynamicSol3D = matUsLD3D(dofUyA,:) ;
%md analogosuly $u_z$ of A node is:  
UzADynamicSol3D = matUsLD3D(dofUyA + 2,:) ;
%md
%md Open figure and plot
fig5 = figure(5) ;
hold on,  grid on
% legend first point plot uy 
plot(timeVecLD3d(1), UyADynamicSol3D(1),...
     'color', 'b', 'linewidth', lw, 'linestyle', '-','markersize', ms, 'marker', 'o')
% legend first point plot uz
plot(timeVecLD3d(1), UzADynamicSol3D(1),...
     'color', 'r', 'linewidth', lw, 'linestyle', '-','markersize', ms, 'marker', '^')
% markers plot uy
plot(timeVecLD3d(1:13:end), UyADynamicSol3D(1:13:end),...
     'color', 'b', 'linewidth', lw, 'linestyle', 'none','markersize', ms, 'marker', 'o')
% continium line plot uy
plot(timeVecLD3d, UyADynamicSol3D,...
     'color', 'b', 'linewidth', lw, 'linestyle', '-', 'marker', 'none')
% markers plot uz
plot(timeVecLD3d(1:23:end), UzADynamicSol3D(1:23:end),...
     'color', 'r', 'linewidth', lw, 'linestyle', 'none','markersize', ms, 'marker', '^')
% continium line plot uz
plot(timeVecLD3d, UzADynamicSol3D,...
     'color', 'r', 'linewidth', lw, 'linestyle', '-', 'marker', 'none')
legend('U_y node A', 'U_z node A' )
labx=xlabel('x[m]');    laby=ylabel('Dispalcements[m]');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location','northEast' ) ;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize ) ;
set(labx, 'FontSize', axisFontSize); set(laby, 'FontSize', axisFontSize) ;
namefig5 = strcat(folderPathFigs, 'uA3D.png') ;
print(fig5, namefig5,'-dpng') ;
%md
%md```@raw html
%md<img src="../../assets/nonLinearAerodynamics/uA3D.png" alt="plot check angular displacements" width="500"/>
%md```
%md
%md Finally a GIF to illustrate the motion amplitude is subsequently presented:
%md
%md```@raw html
%md<img src="../../assets/nonLinearAerodynamics/cyilindricalCantBeam3D.gif" alt="plot check angular displacements" width="500"/>
%md```
%md