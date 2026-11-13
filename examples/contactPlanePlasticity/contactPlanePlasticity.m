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
% md# Frictionless contact in plane problems with plasticity
% md
% md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS/blob/master/examples/contactPlanePlasticity/contactPlanePlasticity.m)
% md
% md In this example a plane strain block of dimensions $L_x \times L_y$ and thickness $t$ rests
% md on a rigid frictionless floor placed at $y=0$, and is submitted to a uniform pressure $p$
% md applied on its top face. The displacements in the $x$ direction are fixed on both vertical
% md faces, therefore the block undergoes a state of *confined compression* (also called uniaxial
% md strain), for which a closed-form solution is available.
% md
% md The support provided by the floor is introduced using `contact-node` elements, which
% md regularize the unilateral contact constraint with the penalty method. The purpose of the
% md example is twofold: to verify the contact element against the analytic solution of the
% md confined compression problem, and to verify that it behaves consistently when the material
% md of the block yields.
% md
% md## Analytic solution
% md
% md Since $\epsilon_{xx}=0$ (confined) and $\epsilon_{zz}=0$ (plane strain), the only non-zero
% md strain component is $\epsilon_{yy}$, and the plane strain constitutive matrix gives
% md```math
% md \sigma_{yy} = \dfrac{E(1-\nu)}{(1+\nu)(1-2\nu)} \epsilon_{yy} ,
% md\qquad
% md \sigma_{xx} = \dfrac{\nu}{1-\nu}\sigma_{yy} .
% md```
% md Equilibrium with the applied pressure requires $\sigma_{yy} = -p$, hence the vertical
% md displacement of the top face is given by
% md```math
% md u_y(L_y) = -\dfrac{p(1+\nu)(1-2\nu)}{E(1-\nu)} L_y .
% md```
% md Moreover, the resultant of the contact reactions must balance the applied load, that is
% md $R = p L_x t$, independently of the constitutive model considered. Both relations are used
% md below as verification measures.
% md
% md Before defining the structs, the workspace is cleaned and the ONSAS directory is added to
% md the path.
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes') % hidden
  clear all;
end % hidden
% add path
addpath(genpath([pwd '/../../src']));
% md
% md## Scalar parameters
% md
% md geometry and material of the block
E       = 210;
nu      = 0.3;
sigmaY0 = 0.24;
H       = 0;
Lx      = 2;
Ly      = 1;
thickness = 1;
% md the applied pressures used for the elastic and the elasto-plastic analyses
pElastic = 0.2;
pPlastic = 0.45;
% md For the confined compression state the von Mises equivalent stress is proportional to
% md $|\sigma_{yy}|$, thus the pressure at which the block starts yielding can be computed as
kLateral = nu / (1 - nu);
qFactor  = sqrt(1.5) * sqrt((2 * kLateral - 1)^2 + (2 - kLateral)^2) / 3;
pYield   = sigmaY0 / qFactor;
fprintf('\nyielding pressure of the confined block: %8.4f\n', pYield);
% md so that `pElastic` remains in the elastic range while `pPlastic` produces plastic strains.
% md The elasto-plastic pressure is set only moderately above the yielding pressure. For
% md considerably larger pressures the von Mises model implemented in `elementTriangSolid`
% md stops representing the confined compression state, since its yielding function neglects
% md the out of plane deviatoric stress component and therefore bounds the vertical stress that
% md the material can carry, which equilibrium fixes at the applied pressure.
% md
% md The analytic vertical displacement of the top face is proportional to the applied pressure
% md while the block remains elastic, with the compliance given by the expression obtained above
complianceElastic = (1 + nu) * (1 - 2 * nu) / (E * (1 - nu)) * Ly;
% md
% md The penalty parameter of the contact elements is set as a large multiple of the stiffness
% md of the solid. Since the penetration of a penalized contact is inversely proportional to
% md this parameter, two values are used below to verify that convergence rate.
penaltyParam = 1e3 * E * thickness;
% md
% md## Mesh
% md
% md A structured mesh of `nx` by `ny` quadrilateral cells is generated, each of them split into
% md two `triangle` elements.
nx = 8;
ny = 4;
xCoords = linspace(0, Lx, nx + 1);
yCoords = linspace(0, Ly, ny + 1);
% md the coordinates matrix is built looping in the grid, where the node placed at the column
% md `i` and the row `j` of the grid gets the number `(j-1)*(nx+1)+i`
mesh = struct();
mesh.nodesCoords = zeros((nx + 1) * (ny + 1), 3);
for j = 1:(ny + 1)
  for i = 1:(nx + 1)
    mesh.nodesCoords((j - 1) * (nx + 1) + i, :) = [xCoords(i) yCoords(j) 0];
  end
end
% md the nodes placed on the bottom and on the top faces of the block are collected, as well as
% md the nodes of the two vertical faces
bottomNodes = 1:(nx + 1);
topNodes    = (ny * (nx + 1) + 1):((ny + 1) * (nx + 1));
lateralNodes = [1:(nx + 1):(ny * (nx + 1) + 1)  (nx + 1):(nx + 1):((ny + 1) * (nx + 1))];
% md
% md The connectivity cell is assembled in three groups. First the `edge` elements of the top
% md face, used to introduce the pressure load, and the `node` elements used to fix the
% md horizontal displacements of the vertical faces. Both groups are defined with no material,
% md since they are only used to introduce boundary conditions.
conecCell = {};
for i = 1:nx
  conecCell{end + 1} = [0 1 1  topNodes(i) topNodes(i + 1)];
end
for i = 1:length(lateralNodes)
  conecCell{end + 1} = [0 2 2  lateralNodes(i)];
end
% md then the `triangle` elements of the solid, with the counter-clockwise orientation required
% md by the element
for j = 1:ny
  for i = 1:nx
    nodeSW = (j - 1) * (nx + 1) + i;
    nodeSE = nodeSW + 1;
    nodeNE = j * (nx + 1) + i + 1;
    nodeNW = nodeNE - 1;
    conecCell{end + 1} = [1 3 0  nodeSW nodeSE nodeNE];
    conecCell{end + 1} = [1 3 0  nodeSW nodeNE nodeNW];
  end
end
% md Finally two alternative groups are prepared for the support of the bottom face. The first
% md one fixes the vertical displacements using `node` elements and is used to compute the
% md reference solution, while the second one introduces the `contact-node` elements of the
% md rigid floor. Note that the contact elements do have a material assigned, since the
% md parameters of the rigid plane are read from it.
conecCellDiri    = conecCell;
conecCellContact = conecCell;
for i = 1:length(bottomNodes)
  conecCellDiri{end + 1}    = [0 2 3  bottomNodes(i)];
  conecCellContact{end + 1} = [2 4 0  bottomNodes(i)];
end
% md
% md## MEB parameters
% md
% md### materials
% md The first material is the one of the solid, set below for each analysis case. The second
% md one holds the parameters of the rigid plane of the contact elements, given as
% md `[ penaltyParam  x0  y0  nx  ny ]`, where `(x0,y0)` is a point of the plane and `(nx,ny)`
% md is its normal vector, pointing towards the body.
materials = struct();
materials(1).modelName   = 'elastic-linear';
materials(1).modelParams = [E nu];
materials(2).modelName   = 'rigid-plane-penalty';
materials(2).modelParams = [penaltyParam  0 0  0 1];
% md
% md### elements
% md Four kinds of elements are used: `edge` elements to apply the pressure, `node` elements for
% md the imposed displacements, `triangle` elements for the solid of the block, and
% md `contact-node` elements for the rigid floor.
elements = struct();
elements(1).elemType           = 'edge';
elements(1).elemCrossSecParams = thickness;
elements(2).elemType           = 'node';
elements(3).elemType           = 'triangle';
elements(3).elemTypeParams     = 2; % plane strain
elements(3).elemCrossSecParams = thickness;
elements(4).elemType           = 'contact-node';
% md
% md### boundaryConds
% md The first boundary condition is the pressure applied on the top face, introduced as a
% md global load per unit length on the `edge` elements. The remaining two fix the horizontal
% md displacements of the vertical faces and the vertical displacements of the bottom face.
boundaryConds = struct();
boundaryConds(1).loadsCoordSys = 'global';
boundaryConds(1).loadsTimeFact = @(t) t;
boundaryConds(1).loadsBaseVals = [0 0 -pElastic 0 0 0];
%
boundaryConds(2).imposDispDofs = [1];
boundaryConds(2).imposDispVals = 0;
%
boundaryConds(3).imposDispDofs = [3];
boundaryConds(3).imposDispVals = 0;
% md
% md### initialConds
% md No non-homogeneous initial conditions are considered.
initialConds = struct();
% md
% md## Analysis parameters
analysisSettings = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.stopTolIts    = 30;
analysisSettings.stopTolDeltau = 1.0e-10;
analysisSettings.stopTolForces = 1.0e-8;
analysisSettings.finalTime     = 1;
analysisSettings.deltaT        = 0.2;
%
otherParams = struct();
otherParams.plots_format = '';
% md
% md## Case 1: elastic block supported by imposed displacements
% md
% md This first case provides the reference solution, where the floor is replaced by imposed
% md vertical displacements on the bottom face.
mesh.conecCell = conecCellDiri;
otherParams.problemName = 'contactPlanePlasticity_elasticDirichlet';
%
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
% md The first node of the top face is used as control node, and the whole history of its
% md vertical displacement is stored together with the load factors, to be plotted below. The
% md load factors of the pressure are placed in the first column of `loadFactorsMat`, since the
% md pressure was introduced as the first boundary condition.
dofTop = (topNodes(1) - 1) * 6 + 3;
uyTopDirichletElastic = matUs(dofTop, end);
%
uyCurveDiriElastic     = matUs(dofTop, :);
lambdaCurveDiriElastic = loadFactorsMat(:, 1);
% md and compared with the analytic solution of the confined compression problem
uyTopAnalyticElastic = -pElastic * complianceElastic;
errDirichletElastic = abs(uyTopDirichletElastic - uyTopAnalyticElastic) / abs(uyTopAnalyticElastic);
% md
% md## Case 2: elastic block supported by the rigid floor
% md
% md The imposed displacements of the bottom face are now replaced by the contact elements.
mesh.conecCell = conecCellContact;
otherParams.problemName = 'contactPlanePlasticity_elasticContact';
%
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
%
uyTopContactElastic = matUs(dofTop, end);
errContactElastic = abs(uyTopContactElastic - uyTopAnalyticElastic) / abs(uyTopAnalyticElastic);
%
uyCurveContactElastic     = matUs(dofTop, :);
lambdaCurveContactElastic = loadFactorsMat(:, 1);
% md Since the point of the rigid plane and its normal are $(0,0)$ and $(0,1)$, the gap of each
% md contact node coincides with its vertical displacement, and the resultant of the contact
% md reactions can be computed from the penalty relation and compared with the applied load.
% md The gaps are truncated with a `min` since the contact is unilateral: an open gap produces
% md no reaction at all, instead of a tensile one.
dofsBottom = (bottomNodes - 1) * 6 + 3;
uyBottomContactElastic = matUs(dofsBottom, end);
reactionContactElastic = -penaltyParam * sum(min(uyBottomContactElastic, 0));
errReactionElastic = abs(reactionContactElastic - pElastic * Lx * thickness) / (pElastic * Lx * thickness);
% md The maximum penetration of the rigid plane is also stored, to verify below that it
% md decreases as the penalty parameter is increased.
penetrationElastic = abs(min(uyBottomContactElastic));
% md
% md## Case 3: convergence of the penalty regularization
% md
% md The same problem is solved with a penalty parameter ten times larger. Since the penetration
% md of a penalized contact is inversely proportional to this parameter, both the penetration
% md and the error of the top displacement are expected to decrease by about one order of
% md magnitude.
materials(2).modelParams = [10 * penaltyParam  0 0  0 1];
otherParams.problemName = 'contactPlanePlasticity_elasticContactStiff';
%
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
%
uyTopContactStiff = matUs(dofTop, end);
errContactStiff = abs(uyTopContactStiff - uyTopAnalyticElastic) / abs(uyTopAnalyticElastic);
penetrationStiff = abs(min(matUs(dofsBottom, end)));
%
uyCurveContactStiff     = matUs(dofTop, :);
lambdaCurveContactStiff = loadFactorsMat(:, 1);
% md restore the penalty parameter used in the remaining cases
materials(2).modelParams = [penaltyParam  0 0  0 1];
% md
% md## Case 4: elasto-plastic block supported by imposed displacements
% md
% md The material of the block is replaced by the von Mises isotropic hardening model and the
% md applied pressure is increased beyond the yielding pressure computed above. The load is
% md applied in a larger number of steps, as required by the return mapping algorithm.
materials(1).modelName   = 'isotropicHardening';
materials(1).modelParams = [E nu H sigmaY0];
boundaryConds(1).loadsBaseVals = [0 0 -pPlastic 0 0 0];
analysisSettings.deltaT = 0.1;
%
mesh.conecCell = conecCellDiri;
otherParams.problemName = 'contactPlanePlasticity_plasticDirichlet';
%
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
%
uyTopDirichletPlastic = matUs(dofTop, end);
% md the accumulated plastic strain of the solid elements is extracted to verify that the block
% md actually yielded
acumPlasStrainDirichlet = max(cell2mat(modelSolutions{end}.previousStateCell(:, 3)));
%
uyCurveDiriPlastic     = matUs(dofTop, :);
lambdaCurveDiriPlastic = loadFactorsMat(:, 1);
% md
% md## Case 5: elasto-plastic block supported by the rigid floor
% md
% md Finally the same elasto-plastic problem is solved using the contact elements.
mesh.conecCell = conecCellContact;
otherParams.problemName = 'contactPlanePlasticity_plasticContact';
% md the results of this case are also exported to vtk files, to allow the visualization of the
% md deformed configuration of the block resting on the rigid floor
otherParams.plots_format = 'vtk';
%
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
%
uyTopContactPlastic = matUs(dofTop, end);
errPlastic = abs(uyTopContactPlastic - uyTopDirichletPlastic) / abs(uyTopDirichletPlastic);
% md the resultant of the contact reactions must balance the applied load also in the
% md elasto-plastic regime, since global equilibrium does not depend on the constitutive model
uyBottomContactPlastic = matUs(dofsBottom, end);
reactionContactPlastic = -penaltyParam * sum(min(uyBottomContactPlastic, 0));
errReactionPlastic = abs(reactionContactPlastic - pPlastic * Lx * thickness) / (pPlastic * Lx * thickness);
% md and the accumulated plastic strain is extracted as in the previous case
acumPlasStrainContact = max(cell2mat(modelSolutions{end}.previousStateCell(:, 3)));
%
uyCurveContactPlastic     = matUs(dofTop, :);
lambdaCurveContactPlastic = loadFactorsMat(:, 1);
% md
% md## Unloading: release of the rigid floor and remnant strain
% md
% md The two cases below apply the pressure up to its maximum value and then remove it, using a
% md load factor function that grows up to $t=1$ and decreases back to zero at $t=2$.
% md
% md The rigid plane is not removed from the model, since the block has no other vertical
% md support and deleting it would leave the tangent stiffness matrix singular, with a free
% md rigid body translation. It is the contact itself that is released: being a unilateral
% md constraint, the reaction of the rigid plane decays to zero as the pressure is removed, and
% md the block is left simply resting on the floor. Both the release of the contact and the
% md remnant displacement of the control node are verified below.
boundaryConds(1).loadsTimeFact = @(t) t * (t <= 1) + (2 - t) * (t > 1);
analysisSettings.finalTime     = 2;
mesh.conecCell = conecCellContact;
% md Both unloading analyses are exported to vtk files, so that the deformed configuration of
% md the block can be followed along the loading and the unloading branches. The rigid plane is
% md included in the output as a line cell, which allows to see the block pressing the floor
% md during the loading and resting on it, with the remnant deformation, at the end.
otherParams.plots_format = 'vtk';
% md
% md### Case 6: elastic block, loading and unloading
% md
% md In the elastic case the block must recover its initial configuration once the pressure is
% md removed, therefore no remnant displacement is expected.
materials(1).modelName   = 'elastic-linear';
materials(1).modelParams = [E nu];
boundaryConds(1).loadsBaseVals = [0 0 -pElastic 0 0 0];
analysisSettings.deltaT  = 0.2;
otherParams.problemName  = 'contactPlanePlasticity_elasticUnload';
%
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
%
uyCurveElasticUnload     = matUs(dofTop, :);
lambdaCurveElasticUnload = loadFactorsMat(:, 1);
% md the maximum and the remnant displacements of the control node are extracted, the first one
% md at the step where the load factor reaches its maximum and the second one at the last step
[~, indPeakElastic] = max(lambdaCurveElasticUnload);
uyPeakElasticUnload   = uyCurveElasticUnload(indPeakElastic);
uyRemnantElasticUnload = uyCurveElasticUnload(end);
% md and the resultant of the contact reactions at the end of the analysis is computed, to
% md verify that the rigid plane is no longer pressing the block
reactionEndElasticUnload = -penaltyParam * sum(min(matUs(dofsBottom, end), 0));
% md
% md### Case 7: elasto-plastic block, loading and unloading
% md
% md In the elasto-plastic case part of the deformation is not recovered, thus a remnant
% md displacement of the control node must be obtained once the pressure is removed.
materials(1).modelName   = 'isotropicHardening';
materials(1).modelParams = [E nu H sigmaY0];
boundaryConds(1).loadsBaseVals = [0 0 -pPlastic 0 0 0];
analysisSettings.deltaT  = 0.1;
otherParams.problemName  = 'contactPlanePlasticity_plasticUnload';
%
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
%
uyCurvePlasticUnload     = matUs(dofTop, :);
lambdaCurvePlasticUnload = loadFactorsMat(:, 1);
%
[~, indPeakPlastic] = max(lambdaCurvePlasticUnload);
uyPeakPlasticUnload    = uyCurvePlasticUnload(indPeakPlastic);
uyRemnantPlasticUnload = uyCurvePlasticUnload(end);
reactionEndPlasticUnload = -penaltyParam * sum(min(matUs(dofsBottom, end), 0));
% md The remnant displacement is compared with the one that would be obtained if the whole
% md deformation at the maximum load were recovered elastically. Since the block is unloaded
% md from a stress state on the yielding surface, the recovered part is the elastic compliance
% md times the maximum pressure.
uyRecoveredElastic = pPlastic * complianceElastic;
% md
% md## Verification
% md
% md The reference solution obtained with imposed displacements reproduces the analytic solution
% md of the confined compression problem exactly, since the corresponding displacement field is
% md linear and therefore contained in the space of the constant strain triangles used.
verifDirichletElastic = errDirichletElastic < 1e-8;
% md The solution obtained with the contact elements differs from the analytic one by the
% md penetration introduced by the penalty regularization, which for the parameter used is about
% md four orders of magnitude smaller than the displacement of the top face.
verifContactElastic = errContactElastic < 1e-3;
% md Increasing the penalty parameter by one order of magnitude reduces both the error and the
% md penetration by about the same factor, as expected for the penalty method.
verifPenaltyConvergence = (errContactStiff < errContactElastic / 5) && ...
                          (penetrationStiff < penetrationElastic / 5);
% md The resultant of the contact reactions balances the applied load in both the elastic and
% md the elasto-plastic analyses.
verifReactions = (errReactionElastic < 1e-5) && (errReactionPlastic < 1e-5);
% md In the elasto-plastic analyses the block yields, and the solution obtained supporting it
% md with the contact elements matches the one obtained with imposed displacements.
verifPlasticity = (acumPlasStrainDirichlet > 0) && (acumPlasStrainContact > 0);
verifContactPlastic = errPlastic < 1e-3;
% md After the pressure is removed the reaction of the rigid plane vanishes in both unloading
% md analyses, since the contact is a unilateral constraint and cannot pull the block.
verifContactRelease = (abs(reactionEndElasticUnload) < 1e-8 * pElastic * Lx * thickness) && ...
                      (abs(reactionEndPlasticUnload) < 1e-8 * pPlastic * Lx * thickness);
% md The elastic block recovers its initial configuration completely, whereas the elasto-plastic
% md one keeps a remnant displacement.
verifElasticRecovery = abs(uyRemnantElasticUnload) < 1e-8 * abs(uyPeakElasticUnload);
verifPlasticRemnant  = (uyRemnantPlasticUnload < 0) && ...
                       (abs(uyRemnantPlasticUnload) > 0.05 * abs(uyPeakPlasticUnload));
% md Moreover, since the block is unloaded elastically from a stress state on the yielding
% md surface, the remnant displacement must equal the one at the maximum load minus the
% md elastically recovered part, which provides a quantitative check of the remnant strain.
errRemnant = abs(uyRemnantPlasticUnload - (uyPeakPlasticUnload + uyRecoveredElastic)) / ...
             abs(uyRemnantPlasticUnload);
verifRemnantValue = errRemnant < 1e-2;
% md
verifBoolean = verifDirichletElastic && verifContactElastic && verifPenaltyConvergence && ...
               verifReactions && verifPlasticity && verifContactPlastic && ...
               verifContactRelease && verifElasticRecovery && verifPlasticRemnant && ...
               verifRemnantValue;
% md
fprintf('\n--- contact in plane problems with plasticity: verification ---\n');
fprintf('analytic top displacement (elastic)            : %12.6e\n', uyTopAnalyticElastic);
fprintf('imposed displacements support (elastic)        : %12.6e  (rel err %8.2e)\n', ...
        uyTopDirichletElastic, errDirichletElastic);
fprintf('contact support (elastic)                      : %12.6e  (rel err %8.2e)\n', ...
        uyTopContactElastic, errContactElastic);
fprintf('contact support, penalty x10 (elastic)         : %12.6e  (rel err %8.2e)\n', ...
        uyTopContactStiff, errContactStiff);
fprintf('max penetration, penalty and penalty x10       : %12.6e  %12.6e\n', ...
        penetrationElastic, penetrationStiff);
fprintf('contact reaction / applied load (elastic)      : rel err %8.2e\n', errReactionElastic);
fprintf('imposed displacements support (plastic)        : %12.6e\n', uyTopDirichletPlastic);
fprintf('contact support (plastic)                      : %12.6e  (rel err %8.2e)\n', ...
        uyTopContactPlastic, errPlastic);
fprintf('contact reaction / applied load (plastic)      : rel err %8.2e\n', errReactionPlastic);
fprintf('max accumulated plastic strain, diri / contact : %12.6e  %12.6e\n', ...
        acumPlasStrainDirichlet, acumPlasStrainContact);
fprintf('--- unloading ---\n');
fprintf('elastic: peak / remnant displacement           : %12.6e  %12.6e\n', ...
        uyPeakElasticUnload, uyRemnantElasticUnload);
fprintf('plastic: peak / remnant displacement           : %12.6e  %12.6e\n', ...
        uyPeakPlasticUnload, uyRemnantPlasticUnload);
fprintf('plastic: expected remnant (peak + elastic rec) : %12.6e  (rel err %8.2e)\n', ...
        uyPeakPlasticUnload + uyRecoveredElastic, errRemnant);
fprintf('contact reaction after unloading, el / pl      : %12.6e  %12.6e\n', ...
        reactionEndElasticUnload, reactionEndPlasticUnload);
fprintf('verifBoolean                                   : %i\n', verifBoolean);
% md
% md## Plots
% md
% md The load factor of the applied pressure is plotted against the vertical displacement of the
% md control node. Two figures are produced, one for each load level considered, since the load
% md factor of each analysis case is relative to a different reference pressure.
lw = 2.0;
ms = 11;
plotfontsize = 10;
% md
% md### Elastic cases
% md
% md The three elastic cases are compared with the analytic solution of the confined compression
% md problem, which is linear in the applied pressure. The three numerical curves are
% md indistinguishable at the scale of the figure, since the penetration introduced by the
% md penalty regularization is about four orders of magnitude smaller than the displacement of
% md the control node.
lambdaAnalytic = linspace(0, 1, 20);
uyAnalyticElastic = -lambdaAnalytic * pElastic * complianceElastic;
%
figure;
hold on;
grid on;
plot(-uyAnalyticElastic, lambdaAnalytic, 'k-', 'linewidth', lw);
plot(-uyCurveDiriElastic, lambdaCurveDiriElastic, 'bo', 'linewidth', lw, 'markersize', ms);
plot(-uyCurveContactElastic, lambdaCurveContactElastic, 'rx', 'linewidth', lw, 'markersize', ms);
plot(-uyCurveContactStiff, lambdaCurveContactStiff, 'g+', 'linewidth', lw, 'markersize', ms);
labx = xlabel('-u_y of the control node [m]');
laby = ylabel('\lambda(t)');
legend('Analytic', 'Imposed disps support', 'Contact support', 'Contact support, penalty x10', ...
       'location', 'SouthEast');
title(sprintf('Elastic block, p = %4.2f', pElastic));
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize);
set(labx, 'FontSize', plotfontsize);
set(laby, 'FontSize', plotfontsize);
print('output/verifContactPlanePlasticityElastic.png', '-dpng');
% md
% md### Elasto-plastic cases
% md
% md For the elasto-plastic cases no closed-form solution is available, therefore the analytic
% md elastic response is included only as a reference, drawn up to the load factor at which the
% md block starts yielding. Beyond that point the numerical curves depart from the elastic line,
% md and the two support conditions remain superposed.
lambdaYield = pYield / pPlastic;
lambdaAnalyticPlastic = linspace(0, lambdaYield, 20);
uyAnalyticPlastic = -lambdaAnalyticPlastic * pPlastic * complianceElastic;
%
figure;
hold on;
grid on;
plot(-uyAnalyticPlastic, lambdaAnalyticPlastic, 'k-', 'linewidth', lw);
plot(pYield * complianceElastic, lambdaYield, 'ks', 'linewidth', lw, 'markersize', ms);
plot(-uyCurveDiriPlastic, lambdaCurveDiriPlastic, 'bo', 'linewidth', lw, 'markersize', ms);
plot(-uyCurveContactPlastic, lambdaCurveContactPlastic, 'rx', 'linewidth', lw, 'markersize', ms);
labx = xlabel('-u_y of the control node [m]');
laby = ylabel('\lambda(t)');
legend('Analytic elastic reference', 'Yielding pressure', 'Imposed disps support', ...
       'Contact support', 'location', 'SouthEast');
title(sprintf('Elasto-plastic block, p = %4.2f', pPlastic));
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize);
set(labx, 'FontSize', plotfontsize);
set(laby, 'FontSize', plotfontsize);
print('output/verifContactPlanePlasticityPlastic.png', '-dpng');
% md
% md### Loading and unloading
% md
% md The last figure shows the load-unload paths of the two analyses supported by the rigid
% md floor. The elastic block returns to the origin, while the elasto-plastic one unloads along
% md a line parallel to the initial elastic branch and reaches zero load factor with a remnant
% md displacement, which is the permanent deformation left by the plastic flow.
figure;
hold on;
grid on;
plot(-uyCurveElasticUnload, lambdaCurveElasticUnload, 'bo-', 'linewidth', lw, 'markersize', ms);
plot(-uyCurvePlasticUnload, lambdaCurvePlasticUnload, 'rx-', 'linewidth', lw, 'markersize', ms);
plot(-uyRemnantPlasticUnload, 0, 'ks', 'linewidth', lw, 'markersize', ms);
labx = xlabel('-u_y of the control node [m]');
laby = ylabel('\lambda(t)');
legend(sprintf('Elastic block, p = %4.2f', pElastic), ...
       sprintf('Elasto-plastic block, p = %4.2f', pPlastic), 'Remnant displacement', ...
       'location', 'SouthEast');
title('Loading and unloading on the rigid floor');
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize);
set(labx, 'FontSize', plotfontsize);
set(laby, 'FontSize', plotfontsize);
print('output/verifContactPlanePlasticityUnload.png', '-dpng');
