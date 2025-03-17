% Copyright 2024, ONSAS Authors (see documentation)
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
%
% md# Reconfiguration cantilever beam example
% md
% md [![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS/blob/master/examples/dragBeamReconfiguration/dragBeamReconfiguration.m)
% md
% md In this tutorial, a cantilever beam submitted to a flow producing drag forces is considered. The main goal is to validate the drag forces computation considering reconfiguration in a problem with large displacements. The example is based on one of the problems considered in [this reference](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/drag-reduction-of-flexible-plates-by-reconfiguration/A3A066DE3886C40B1463508C92D1C81D), where a reference solution is presented and validated with experimental data. The reference solution was generated with the code publicly available in [this repository](https://github.com/lm2-poly/Reconfiguration-Beam/blob/master/reconfiguration.m)
% md
% md The problem consists in a cantilever beam submitted to a fluid flow with uniform velocity $v_a$($x$,$t$) = $v_a$(t) $c_2$, as shown int he Figure bellow. The beam is clamped on the boundary at $x=0$ m, and the span length is $L$. The cross-section of the beam is circular with diameter $d$. For the material of the beam a linear elastic isotropic model is considered, with Young modulus $E$ density $\rho$.
% md
% md```@raw html
% md<img src="../../assets/dragBeamReconfiguration/ilus.svg" alt="general sketch" width="450"/>
% md```
% md## Dimensionless analysis
% md--------------------
% md
% mdThe problem can be studied through the following dimensionless numbers:
% md
% md```math
% md c_y = \frac{\rho_f L^3 v_a^2}{16 E I_{zz}}, \qquad \mathcal{R} = \frac{F}{\frac{1}{2}\rho_f L d c_d v_a^2}
% md```
% md where  $F$ is the global drag force towards $c_2$, $c_y$ is the Cauchy number that describes the ratio between the stiffness of the beam and the flow load and the reconfiguration number $\mathcal{R}$ reflects the geometric nonlinear effect by dividing the drag of the flexible beam to that of a rigid one of the same geometry
% md## Numerical solution
% md---------------------
% md Before defining the structs, the workspace is cleaned and the ONSAS directory is added:
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes')
  clear all;
end
addpath(genpath([pwd '/../../src']));
% md
% md The problem parameters are loaded:
[L, d, Izz, E, nu, rhoS, rhoF, nuF, ~, NR, cycd_vec, uydot_vec] = loadParametersCirc();
% md where $\rho_f$ and $\rho_s$ are the fluid and solid densities respectively, $\nu$ is the fluid eenematic viscosity, and  $NR$ is the number of load steps (or velocity cases solved):
% md
% md
% md The number of elements employed to discretize the beam is:
numElements = 10;
% md
% md### MEB parameters
% md
% md### materials
% md Since the example contains only one material and co-rotational strain element so then `materials` struct is:
materials             = struct();
materials.modelName   = 'elastic-rotEngStr';
materials.modelParams = [E nu];
materials.density     = rhoS;
% md
% md### elements
% md
% md Two different types of elements are considered, node and frames. The nodes will be assigned in the first entry (index $1$) and the beam at the index $2$. The _elemType_ field is then:
elements(1).elemType = 'node';
elements(2).elemType = 'frame';
% md for the geometries, the node has not geometry to assign (empty array), and frame elements will be set as a circular section with $d$ diameter.
elements(2).elemCrossSecParams{1, 1} = 'circle';
elements(2).elemCrossSecParams{2, 1} = [d];
% md where the aerodynamic coefficents and the chord vector are set by default for `'circle'` cross sections type. For the validation case a constant drag coefficients $c_d =1.2" is used, this is defined in the `'dragCircular'` function:
elements(2).aeroCoefFunctions = {'dragCircular', [], []};
% md
% md### boundaryConds
% md
% md Only one welded (6 degrees of freedom are set to zero) boundary condition (BC) is considered:
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];
% md
% md### initial Conditions
% md Any non-homogeneous initial condition (IC) are set for this case, then an empty struct is used:
initialConds = struct();
% md
% md### mesh parameters
% mdThe coordinates of the mesh nodes are given by the matrix:
mesh = struct();
mesh.nodesCoords = [(0:(numElements))' * L / numElements  zeros(numElements + 1, 2)];
% mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEB parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { };
% md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. Finally the node is assigned:
mesh.conecCell{ 1, 1 } = [0 1 1   1];
% md Next the frame elements MEB parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
for i = 1:numElements
  mesh.conecCell{ i + 1, 1 } = [1 2 0   i i + 1];
end
% md
% md### analysisSettings
% md
% md The fluid properties are set into _fluidProps_ field into `analysisSettings` struct. In this field the fluid velocity, viscosity and density are defined, This will apply a external fluid loads according to the quasi-steady theory for each element with aerodynamic coefficients fields into the `elements` struct. The name of the fluid velocity function located on the same example path is introduced as a string `'windVelCircStatic'`:
analysisSettings = struct();
analysisSettings.fluidProps = {rhoF; nuF; 'windVelCircStatic'};
% md The geometrical non-linear effects and the aerodynamic stiffness matrix are considered in this case to compute the aerodynamic force vector:
analysisSettings.geometricNonLinearAero = true;
analysisSettings.computeAeroStiffnessMatrix = true;
% md since this problem is static, then a N-R method is employed. The time step `deltaT` is 1 since the time here is an index in the fluid velocity vector `uydot_vec` untile the length of this vector `NR` is reached:
analysisSettings.deltaT        =   1;
analysisSettings.finalTime     =   NR;
analysisSettings.methodName    = 'newtonRaphson';
% md Next the maximum number of iterations per load(time) step, the residual force and the displacements tolerances are set to (if null tolerance is set the criterion is not considered):
analysisSettings.stopTolDeltau =   0;
analysisSettings.stopTolForces =   1e-8;
analysisSettings.stopTolIts    =   50;
% md
% md### otherParams
% md The name of the problem and vtk format output are selected:
otherParams  = struct();
otherParams.problemName = 'staticReconfigurationCircleUserDefinedDrag';
otherParams.plots_format = 'vtk';
% md
% md### Case 1: validation case with constant $c_d = 1.2$
% md
% md The ONSAS software is executed for the parameters above defined and the displacement solution of each load(time) step is saved in `matUsCase1`matrix:
% md
% md The reaction forces at $x=0$ are stored  in the following vector and declaring a global variable:
global globalReactionForces
globalReactionForces = zeros(6 * analysisSettings.finalTime, 1);
% md
% md The node index where the reaction forces are computed is:
global glboalNodeReactionForces
glboalNodeReactionForces = 1;
% md
% md### Numeric solution
% md
% md Execute ONSAS and save the results:
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUsCase1, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
% md

% md
% md## Verification
% md---------------------
% md The numerical solution is extracted. First the reference coordinaes
xref = mesh.nodesCoords(:, 1);
yref = mesh.nodesCoords(:, 2);
zref = mesh.nodesCoords(:, 3);
% md Then the values of $R$ and $C_y C_d$ are computed :
numLoadSteps = size(matUsCase1, 2);
timeVec = linspace(0, analysisSettings.finalTime, numLoadSteps);
% initialize vectors
Cy = zeros(numLoadSteps - 1, 1);
R  = zeros(numLoadSteps - 1, 1);
C_d = feval(elements(2).aeroCoefFunctions{1}, 0, 0);
% fill them
for windVelStep = 1:numLoadSteps - 1
  % Compute dimensionless magnitudes
  windVel         = feval(analysisSettings.fluidProps{3, :}, 0, timeVec(windVelStep + 1));
  normWindVel     = norm(windVel);
  dirWindVel      = windVel / normWindVel;
  Cy(windVelStep) =  1 / 2 * rhoF * normWindVel^2 * (L)^3 * d / (E * Izz);

  % numeric drag
  FReaction = globalReactionForces((windVelStep - 1) * 6 + 1:windVelStep * 6);
  FDragi = FReaction(3);
  FDRef  = 1 / 2 * rhoF * normWindVel^2 * C_d * d * L;
  R(windVelStep) =  abs(FDragi) / (FDRef);

end
% md
% md### Gosselin et.Al 2010 solution
% md
% resudrag (cycd, R) and def wich contains de deformed configuration for 10^i cycyd values:
base_dir = '';
if strcmp(getenv('TESTS_RUN'), 'yes') && isfolder('examples')
  base_dir = ['.' filesep 'examples' filesep 'dragBeamReconfiguration' filesep];
end
load([base_dir 'Gosselin2010_data.mat'], 'def', 'resudrag');
% md
% md### Validation plots
% md
% md The plot parameters are:
lw = 4;
ms = 8;
axislw = 1;
axisFontSize = 20;
legendFontSize = 15;
curveFontSize = 15;
Gline = 'k-';
ONSASline = 'bo';
ONSASlineBuiltInDrag = 'rx';
% md
% md The modified Cauchy number vs R is plotted:
% md
fig1 = figure(1);
loglog(C_d * Cy, R, ONSASline, 'linewidth', lw, 'markersize', ms);
hold on;
loglog(resudrag(:, 1), resudrag(:, 2), Gline, 'linewidth', lw, 'markersize', ms);
% add legend
legend('ONSAS', 'Gosselin2010');
% set labels legend
labx = xlabel(' c_y*c_d ');
laby = ylabel('R');
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location', 'eastOutside');
% set fonts
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize);
set(labx, 'FontSize', axisFontSize);
set(laby, 'FontSize', axisFontSize);
grid on;
if length(getenv('TESTS_RUN')) > 0 && strcmp(getenv('TESTS_RUN'), 'yes')
  fprintf('\ngenerating output png for docs.\n');
  figure(1);
  print('output/RvsCyCd.png', '-dpng');
else
  fprintf('\n === NOT in docs workflow. ===\n');
end
% md
% md Deformed configurations for different cauchy numbers
% md
fig3 = figure(3);
% plot reference configuration
plot(xref, yref, 'k--', 'linewidth', lw, 'markersize', ms);
hold on;
grid on;
% add markers for the legend
plot(0, 0,  ONSASline, 'linewidth', lw, 'markersize', ms);
plot(0, 0,  ONSASlineBuiltInDrag, 'linewidth', lw, 'markersize', ms);
% plot deformed configurations for different fluid load steps
for nr = 1:analysisSettings.finalTime
  % Numerical deformed coordinates solution
  xdef = xref + matUsCase1(1:6:end, nr + 1);
  ydef = yref + matUsCase1(3:6:end, nr + 1);
  zdef = zref + matUsCase1(5:6:end, nr + 1);
  thetaXdef = matUsCase1(2:6:end, nr + 1);
  thetaYdef = matUsCase1(4:6:end, nr + 1);
  thetaZdef = matUsCase1(6:6:end, nr + 1);
  % Gosselin deformed coordinates solution
  xdefG = def(1, :, nr);
  ydefG = -def(2, :, nr);
  % Plot
  plot(xdef,  ydef,  ONSASline, 'linewidth', lw, 'markersize', ms);
  plot(xdefG, ydefG,  Gline, 'linewidth', lw, 'markersize', ms);
end
% md
% md
% md```@raw html
% md<img src="../../assets/generated/RvsCyCd.png" alt="plot check deformed configurations" width="500"/>
% md```
% md
% md### Case 2: case with drag coefficient formula proposed in [this reference](https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0000722)
% md
% md  Once `_elemCrossSecParams_` is defined in the `elements` struct then the drag lift and pitch moment are defined with bluit-in functions. As consequence if we want to use the default drag coefficient the `dragCoefFunction` field must be set:
elements(2).aeroCoefFunctions = {[], [], []};
% md moreover the drag formulation in [this reference](https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0000722) is valid with $Re < 2$ x $10^5$. Then the final load step, which is equivalent to the index in the velocity vector `uydot_vec`, is set to 5:
analysisSettings.finalTime = NR - 2;
% md
% md### otherParams
% md The name of the problem and vtk format output are selected:
otherParams  = struct();
otherParams.problemName = 'staticReconfigurationCircleBuiltInDrag';
% md
% md### Numeric solution
% md
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUsCase2, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
% md
% md the report is generated
outputReport(modelProperties.outputDir, modelProperties.problemName);

% md
% md Deformed configurations for different cauchy numbers
% md
for nr = 1:analysisSettings.finalTime
  % Numerical deformed coordinates solution
  xdefCase2 = xref + matUsCase2(1:6:end, nr + 1);
  ydefCase2 = yref + matUsCase2(3:6:end, nr + 1);
  zdefCase2 = zref + matUsCase2(5:6:end, nr + 1);
  thetaXdefCase2 = matUsCase2(2:6:end, nr + 1);
  thetaYdefCase2 = matUsCase2(4:6:end, nr + 1);
  thetaZdefCase2 = matUsCase2(6:6:end, nr + 1);
  % Plot
  plot(xdefCase2,  ydefCase2,  ONSASlineBuiltInDrag, 'linewidth', lw, 'markersize', ms);
end
% add legend
legend('Gosselin2010(cd = 1.2)', 'ONSAS(cd = 1.2)', 'ONSAS(Built-in c_d)');
% add legends and labels
labx = xlabel('x [m]');
laby = ylabel('y [m]');
% set fonts
set(legend, 'linewidth', axislw, 'fontsize', legendFontSize, 'location', 'eastOutside');
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize);
set(labx, 'FontSize', axisFontSize);
set(laby, 'FontSize', axisFontSize);
% axis equal
% save fig
if length(getenv('TESTS_RUN')) > 0 && strcmp(getenv('TESTS_RUN'), 'yes')
  fprintf('\ngenerating output png for docs.\n');
  figure(3);
  print('output/defPlots.png', '-dpng');
else
  fprintf('\n === NOT in docs workflow. ===\n');
end
% md
% md```@raw html
% md<img src="../../assets/generated/defPlots.png" alt="plot check deformed configurations" width="500"/>
% md```
% md
% md### Case 3: Dynamic Analysis Considering Vortex-Induced Vibrations (VIV)
% md```@raw html
% md<img src="../../assets/dragBeamReconfiguration/VIVilus.svg" alt="general sketch VIV case" width="450"/>
% md```
% md
% md In this case, we extend the analysis to include dynamic effects and Vortex-Induced Vibrations (VIV) on the cantilever beam, according to the formulation proposed in [(Villié, et.al., 2024)](https://www.sciencedirect.com/science/article/abs/pii/S0889974624000094).
% md
% md### elements
% md
% md The beam is subjected to a constant fluid flow in $c_2$, which induces drag and lift forces. The aerodynamic coefficients for drag and lift are defined using the `dragCircular` and `liftCoefVIV` functions, respectively.
elements(2).aeroCoefFunctions = {'dragCircular', 'liftCoefVIV', []};
% md
% md### initial Conditions
% md Initial conditions are derived from Case 1 at step `NR - 3` of the previous static analysis, inducing VIV effects in the static equilibrium for the same drag forces. Null velocity and accelerations are set by default:
initialConds    = struct();
initialConds.U  = matUsCase1(:,  end - 3);
% md We set initial values for the in-line and cross-flow wake variables:
dofsWakeVariablesPerElement = 2;
elementQ0 = (2 * ((1:numElements)' / numElements) - 1) * 0.001;
initialConds.Q0 = repelem(elementQ0, dofsWakeVariablesPerElement);
elementP0 = (2 * ((1:numElements)' / numElements) - 1) * 0.002;
initialConds.P0 = repelem(elementP0, dofsWakeVariablesPerElement);
% md
% md### analysisSettings
% md
% md The $\alpha$-HHT algorithm is set with the following tolerances, time step, and final time:
analysisSettings = struct();
analysisSettings.finalTime = 0.4;
% analysisSettings.finalTime =   5.0; % to reproduce animation
analysisSettings.deltaT = 0.01;
analysisSettings.methodName = 'alphaHHT';
analysisSettings.stopTolIts = 15;
analysisSettings.stopTolDeltau = 0;
analysisSettings.stopTolForces = 1e-5;
% md The constant velocity field corresponding to `NR - 3` step is set as well as density and viscosity of the fluid:
analysisSettings.fluidProps = {rhoF; nuF; 'windVelCircDynamic'};
% md The following parameters are set to configure the dynamic analysis considering VIV. The `analysisSettings.crossFlowVIVBool` parameter enables the consideration of cross-flow VIV in the analysis. The `analysisSettings.inLineVIVBool` parameter determines whether in-line VIV is considered, and it is set to `true` in this case. Lastly, the `analysisSettings.addedMassBool` parameter accounts for the added mass effect for the fluid forces computation.
analysisSettings.crossFlowVIVBool = true;
analysisSettings.inLineVIVBool = true;
analysisSettings.addedMassBool = true;
% md### otherParams
% md The name of the problem is set and VTK output is exported to generate an animation:
otherParams = struct();
otherParams.problemName = 'vivBeamReconfiguration';
otherParams.plots_format = '';
% md
% md### Numeric solution
% md
% Initialize and solve the dynamic problem
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% Perform the dynamic analysis
[matUsDynamic, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
% md
% md Evolution of displacements of node A
% md
finalNodeIndex  = numElements + 1;
zDisplacement   = matUsDynamic(5:6:end, :);
timeVector      = 0:analysisSettings.deltaT:analysisSettings.finalTime;
% Generate figure
figure;
plot(timeVector, zDisplacement(finalNodeIndex, :), 'b-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Displacement in z [m]');
title('Displacement in z at A');
grid on;
set(gca, 'linewidth', axislw, 'fontsize', curveFontSize);
set(get(gca, 'xlabel'), 'FontSize', axisFontSize);
set(get(gca, 'ylabel'), 'FontSize', axisFontSize);
% Save figure for automatic deployment
if length(getenv('TESTS_RUN')) > 0 && strcmp(getenv('TESTS_RUN'), 'yes')
  fprintf('\ngenerating output png for docs.\n');
  print('output/zDisplacementVIV.png', '-dpng');
else
  fprintf('\n === NOT in docs workflow. ===\n');
end
% md
% md```@raw html
% md<img src="../../assets/generated/zDisplacementVIV.png" alt="Displacemnts in z direction of node A" width="500"/>
% md```
% md
% md## Verification boolean
% md
% The verification boolean is computed considering the deformed configurations, the cycd curve
% deformed coordinates dif norm and the viv movement as well.
vecDifDeform =  [norm(ydef - ydefG(1:numElements * 10:end)'); ...
                 norm(xdef - xdefG(1:numElements * 10:end)')];
% verification boolean deformed
verifBooleanDef =  vecDifDeform <=  2e-2 * L;
% cycd vs R verification boolean is:
verifBooleanR = abs(R(end) - resudrag(end, 2)) <  5e-3;
% viv boolean verification:
verifBooleanV = abs(zDisplacement(end, 38) - (-7.4e-03)) <  (7.4e-03*1e-3);
% The total verifboolean is:
verifBoolean = verifBooleanR && verifBooleanV && all(verifBooleanDef);
disp(verifBooleanR);
disp(verifBooleanV);
disp(verifBoolean);
