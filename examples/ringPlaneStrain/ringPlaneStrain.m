% Copyright 2025, ONSAS Authors (see documentation)
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
% md# Plane strain ring example
% md
% md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS/blob/master/examples/ringPlaneStrain/ringPlaneStrain.m)
% md
% md In this example a hollow cylinder submitted to an internal pressure $p_i$ as shown in diagram depicted below is considered. The length of the cylinder is $L_z $ m and the internal and external radious are $R_i$ and $R_e$, respectively.
% md
% md```@raw html
% md<img src="../../assets/linearCylinderPlaneStrain/ilusCylinderPlaneStrain.svg" alt="linear cylinder diagram" width="500"/>
% md```
% md
% md Before defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar geometry and material parameters are defined:
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes')
  clear all;
end
% add path
addpath(genpath([pwd '/../../src']));
% scalar parameters
% E = 1e6 ; nu = 0.3 ; p = 30e3 ; L = .75 ; Re = 0.15 ; Ri = 0.1 ;
E = 210;
nu = 0.3;
p = 0.01;
L = .75;
global Re
global Ri
Re = 200;
Ri = 100;
% md## Linear analysis
% md### Analytic solution
% md
% md The solution displacement field is extracted from chapter 4 of  (Timoshenko and Goodier, Theory of Elasticity, 3rd edition). The Navier's equation, imposing no temperature variation, no volumetric forces, and considering a radial dispalcement field leads to:
% md```math
% md  \nabla (\nabla . \textbf{u}(r,\theta,z) )  = 0
% md```
% md Due to the symmetry of the problem $\mathbf{\mathit{u_{\theta}}} = 0 $ and also $\mathbf{ \mathit{ \textbf{u} (r,\theta,z) } } = \mathbf{ \mathit{ \textbf{u}(r,z) } } $. Thus, according to the boundary conditions stated above $\mathit{u_z(r,z)=0}$ and the radial displacements field $\mathit{u_r(r)}$ only varies with $r$. Thereafter by imposing the boundary conditions stated above and substituting ($E$, $\nu$) into Lamé parameters ($\lambda=\frac{ E\nu }{(1 + 2\nu )(1 - 2\nu )}$ and $\mu=\frac{ E\nu }{(1 + 2\nu )}$) we obtain:
% md```math
% md u_r(r) = Ar + \dfrac{B}{r}  \\
% md A = \dfrac{(1+\nu)(1-2\nu)R_i^2p_i}{E(R_e^2-R_i^2)}, \quad
% md B = \dfrac{(1+\nu)R_i^2R_e^2p_i}{E(R_e^2-R_i^2)}
% md```
% md### Numerical solution
% md
% md#### MEB parameters
% md
% md##### materials
% md The constitutive behavior of the material considered is isotropic linear elastic.
% md Since only one material is considered, the structs defined for the materials contain only one entry:
materials = struct();
materials.modelName  = 'elastic-linear';
materials.modelParams =  [E nu];
% md
% md##### elements
% md
% md In this plane model, three kinds of elements are used: `triangle` for the solid, `edges` to add pressure loads and `nodes` to set additional boundary conditions for the numerical resolution. Since three kinds of elements are used, the struct has length 3:
elements = struct();
elements(1).elemType           = 'node';
elements(2).elemType           = 'edge';
elements(2).elemCrossSecParams = L;
elements(3).elemType           = 'triangle';
elements(3).elemTypeParams     = 2;
elements(3).elemCrossSecParams = L;
% md where `elemCrossSecParams` field sets the thickness of the edge and `elemTypeParams` sets the plane strain triangle element.
% md
% md##### boundaryConds
% md Three BCs are considered, one corresponding to a load and two for displacements.
% md The first two BCs constrain displacements in $x$ and $y$ global directions respectively:
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [1];
boundaryConds(1).imposDispVals = [0];
boundaryConds(2).imposDispDofs = [3];
boundaryConds(2).imposDispVals = [0];
% md then the third BC corresponds to the normal pressure.
% md This is introduced in `local` coordinates. The first entry is the force along the local x coordinate of the edge
% md (tangent), the second is the moment along that direction and the third is the force towards the normal vector obtained by rotating the tangent vector 90 degrees in global axis $z$:
boundaryConds(3).loadsCoordSys = 'local';
boundaryConds(3).loadsTimeFact = @(t) t;
boundaryConds(3).loadsBaseVals = [0 0 p 0 0 0];
% md
% md
% md#### Mesh
% md The mesh can be read from the msh file. However, if any changes to the mesh are desired, the .geo file can be edited and the msh file can be re-generated using GMSH.
% md
% md```@raw html
% md<img src="../../assets/linearCylinderPlaneStrain/meshCylinderPlaneStrain.png" alt="mesh plot" width="500"/>
% md```
% md The element properties are set using labels into GMSH follwing the MEB nomenclature. First `triangle` elements have linear elastic material so entry $1$ of the _materialṣ_ struct is assigned. Then for both `node` and `edge` elements any material is set.
% md Next displacement boundary conditions are assigned to the element, since the problem is modeled into $x-y$ plane, a constrain to avoid rotation along $z$ is necessary. This is done fixing $y$ and $x$ displacements (using `boundaryConds(1)` and `boundaryConds(2)` as labels) on points 2 3 4 5.
% md Finally the internal pressure is applied on the `edge` elements linked with curves from one to four (Circles 1-4 in Figure). In accordance with the orientation of the curve set in GMSH, the normal vector obtained in local coordinates is $e_r$ so the internal pressure is assigned using `boundaryConds(3)`. Once the mesh is created is read using:
base_msh = '';
if strcmp(getenv('TESTS_RUN'), 'yes') && isfolder('examples')
  base_msh = ['.' filesep 'examples' filesep 'ringPlaneStrain' filesep];
end
mesh = struct();
[mesh.nodesCoords, mesh.conecCell] = meshFileReader([base_msh 'ring.msh']);
% md
% md#### initialConds
% md Any non-homogeneous initial conditions are considered, thereafter an empty struct is set:
initialConds = struct();
% md#### Analysis parameters
% md
% md The Newton-Raphson method is employed to solve 2 load steps. The ratio between `finalTime` and `deltaT` sets the number of load steps used to evaluate `boundaryConds(3).loadsTimeFact` function:
analysisSettings = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.stopTolIts    = 30;
analysisSettings.stopTolDeltau = 1.0e-12;
analysisSettings.stopTolForces = 1.0e-12;
analysisSettings.finalTime     = 1;
analysisSettings.deltaT        = .5;
% md
% md#### Output parameters
% md
otherParams = struct();
otherParams.problemName = 'linear_PlaneStrain';
otherParams.plots_format = 'vtk';
% md The ONSAS software is executed for the parameters defined above and the displacement solution of each load(time) step is saved in `matUs`matrix:
% md
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

% md
% md### Verification
% mdThe numerical and analytic solutions are compared at the final load step for the internal and external surface (since all the elements on the same surface have the same analytic solution):
% internal surface analytic solution
A = (p * (1 + nu) * (1 - 2 * nu) * Ri^2) / (E * (Re^2 - Ri^2));
B = (p * (1 + nu) * Ri^2 * Re^2)   / (E * (Re^2 - Ri^2));
analyticValRi = A * Ri + B / Ri;
% internal surface numerical solution
dofXRi = 1;
numericalRi = matUs(dofXRi, end);
% external surface analytic solution
analyticValRe = A * Re + B / Re;
% external surface numerical solution
dofXRe = (8 - 1) * 6 + 3;
numericalRe = matUs(dofXRe, end);
% md
% md The numerical and analytical solution for the internal and external surface are plotted:
% plot parameters
lw = 2.0;
ms = 11;
plotfontsize = 10;
figure;
hold on;
grid on;
% internal surface
plot(matUs(dofXRi, :), loadFactorsMat(:, 3), 'ro', 'linewidth', lw, 'markersize', ms);
plot(linspace(0, analyticValRi, length(loadFactorsMat(:, 3))), loadFactorsMat(:, 3), 'k-', 'linewidth', lw, 'markersize', ms);
% internal surface
plot(matUs(dofXRe, :), loadFactorsMat(:, 3), 'ro', 'linewidth', lw, 'markersize', ms);
plot(linspace(0, analyticValRe, length(loadFactorsMat(:, 3))), loadFactorsMat(:, 3), 'k-', 'linewidth', lw, 'markersize', ms);
labx = xlabel('Displacement [m]');
laby = ylabel('\lambda(t)');
legend('Numeric', 'Analytic', 'location', 'East');
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize);
set(labx, 'FontSize', plotfontsize);
set(laby, 'FontSize', plotfontsize);
print('output/verifLinearRingPlaneStrain.png', '-dpng');
% md
% md```@raw html
% md<img src="../../assets/linearCylinderPlaneStrain/verifLinearCylinderPlaneStrain.png" alt="verification plot" width="500"/>
% md```
% md
% md Finally the deformed configuration is illustrated:
% md```@raw html
% md<img src="../../assets/linearCylinderPlaneStrain/defLinearCylinderPlaneStrain.png" alt="def plot" width="500"/>
% md```
% md
% md## Elastoplastic analysis
% md### Semi-analytic solution
% md The solution is extracted from Hill (The mathematical theroy of plasticity, 1950).
% md The yielding pressure $p_0$ is defined as,
% md```math
% md Y = \dfrac{2\sigma_{Y,0}}{\sqrt{3}}  \\
% md p_0 = \dfrac{Y}{2}\left(1+\dfrac{R_i^2}{R_e^2}\right).
% md```
% md The radial displacement of the outer surface of the ring is given by,
% md```math
% md u_r(R_e) = \text{if}~~ p \leq p_0 \\
% md \dfrac{2 p R_e}{E \left( \dfrac{R_e^2}{R_i^2-1}\right) }( 1-\nu^2 ) \\
% md \text{else} \\
% md \dfrac{2pR_e}{E\left(\dfrac{R_e^2}{R_i^2-1}\right)}(1-\nu^2)
% md```
% md where $c$ denotes the plastic front surface in the ring and is given by the implicit function,
% md```math
% md \dfrac{p}{Y} = ln\left(\dfrac{c}{R_i}\right) + \dfrac{1}{2}\left(1-\dfrac{c^ 2}{R_e^2}\right)
% md```
% md### Numerical solution
% scalar parameters
E = 210;
nu = 0.3;
H = 0;
sigmaY0 = 0.24;
L = .75;
p = 0.01;
% md
% md
% md### MEB parameters
% md
% md#### materials
% md The constitutive behavior of the material considered is isotropic hardening.
% md Since only one material is considered, the structs defined for the materials contain only one entry:
materials.modelName  = 'isotropicHardening';
materials.modelParams =  [E nu H sigmaY0];
% md
% md#### elements
% md
% md The elements struct is the same as the previous model.
% md
% md#### boundaryConds
% md The BC struct is the same as in the elastic-linear case. However the loadsTimeFact function can be modified to consider unloading as follows.
boundaryConds(3).loadsTimeFact = @(t) t * (t <= 19) + (t - (t - 19) * 2) * (t > 19);
% md
% md#### initialConds
% md Any non-homogeneous initial conditions are considered, thereafter the struc is the same as in the previous example.
% md
% md### Mesh
% md The mesh can be read from the msh file. The same mesh as in the elastic-linear case is considered for this problem.
% md
Conec = myCell2Mat(mesh.conecCell);
elems = size(Conec, 1);
% md
% md### Analysis parameters
% md
% md The Newton-Raphson method is employed to solve 19 load steps. The ratio between `finalTime` and `deltaT` sets the number of load steps used to evaluate `boundaryConds(3).loadsTimeFact` function:
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.stopTolIts    = 30;
analysisSettings.stopTolDeltau = 1.0e-8;
analysisSettings.stopTolForces = 1.0e-6;
analysisSettings.finalTime     = 19;
analysisSettings.deltaT        = 1;
% md
% md### Output parameters
% md
otherParams.problemName = 'EPP_PlaneStrain';
otherParams.plots_format = 'vtk';
% md The ONSAS software is executed for the parameters defined above and the displacement solution of each load(time) step is saved in `matUs`matrix:
% md

[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

% md
% md### Verification
% mdThe numerical and analytic solutions are compared for the external surface (since all the elements on the same surface have the same analytic solution):
%
global Y
%
Y = 2 * sigmaY0 / sqrt(3);
% p0 = Y/2 * (1-a^2/b^2)
p0 = Y / 2 * (1 - Ri^2 / Re^2); % Yielding pressure
%
pressure_vals = loadFactorsMat(:, 3) * p;
%
cvals = zeros(length(pressure_vals), 1);
ubAna = zeros(length(pressure_vals), 1);
%
% Plastic front value
for i = 1:length(cvals)
  p = pressure_vals(i);
  if i == 1
    % val = fsolve(@(c)cValue(c,p,Y,a,b), a) ;
    val = fsolve(@(c)cValue(c, p, Y, Ri, Re), Ri);
  else
    % val = fsolve(@(c)cValue(c,p,Y,a,b), cvals(i-1)) ;
    val = fsolve(@(c)cValue(c, p, Y, Ri, Re), cvals(i - 1));
  end
  cvals(i) = val;
end
%
% Analytic radial displacement at outer surface
for i = 1:length(cvals)
  p = pressure_vals(i);
  if p < p0
    % ubAna(i) = 2*p*b / ( E*( b^2/a^2-1 ) ) * (1-nu^2) ;
    ubAna(i) = 2 * p * Re / (E * (Re^2 / Ri^2 - 1)) * (1 - nu^2);
  else
    c = cvals(i);
    % ubAna(i) = Y*c^2/(E*b) * (1-nu^2) ;
    ubAna(i) = Y * c^2 / (E * Re) * (1 - nu^2);
  end
end
% md
% md### Plots
% plot parameters
lw = 2.0;
ms = 11;
plotFontSize = 10;
fig = figure;
hold on;
grid on;
% node to plot the solution
node = 5;
dofX = node * 6 - 5;
ubNum = matUs(dofX, :);
%
plot(ubNum, pressure_vals, 'b-o', 'linewidth', lw, 'markersize', ms);
plot(ubAna, pressure_vals, 'g-x', 'linewidth', lw, 'markersize', ms);
%
legend ({'FEM', 'Analytic'}, 'location', 'east');
labx = xlabel('u_b');
laby = ylabel('p');
tit = title('p-u_b');
set(labx, 'fontsize', plotFontSize * .8);
set(laby, 'fontsize', plotFontSize * .8);
set(tit, 'fontsize', plotFontSize);
% md
% md The numerical solution is verified for both cases:
% md
analyticCheckTolerance = 1e-2;
verifBoolean = ((numericalRi - analyticValRi) < analyticCheckTolerance) && ...
               ((numericalRe - analyticValRe) < analyticCheckTolerance) && ...
               ((ubNum(end) - ubAna(end)) < analyticCheckTolerance);
% md
