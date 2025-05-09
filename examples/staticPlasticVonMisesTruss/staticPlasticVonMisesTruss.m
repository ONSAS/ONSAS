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
% md# Static Von-Mises Truss example
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes')
  clear all;
end
addpath(genpath([pwd '/../../src']));
% scalar parameters
E = 210e9;
A = 2.5e-3;
ang1 = 65;
L = 2;
Kplas = E * .1;
sigma_Y_0 = 25e6;

% x and z coordinates of node 2
x2 = cos(ang1 * pi / 180) * L;
z2 = sin(ang1 * pi / 180) * L;

materials = struct();
materials.modelName  = 'plastic-rotEngStr';
materials.modelParams = [E Kplas sigma_Y_0];

elements = struct();
elements(1).elemType = 'node';
elements(2).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle', sqrt(A * 4 / pi) };

boundaryConds = struct();
boundaryConds(1).imposDispDofs = [1 3 5];
boundaryConds(1).imposDispVals = [0 0 0];

boundaryConds(2).imposDispDofs =   3;
boundaryConds(2).imposDispVals =  0;
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) 3.0e8 * t;
boundaryConds(2).loadsBaseVals = [0 0 0 0 -1 0];

mesh = struct();
mesh.nodesCoords = [0  0   0; ...
                    x2  0  z2; ...
                    2 * x2  0   0];

mesh.conecCell = cell(5, 1);
mesh.conecCell{ 1, 1 } = [0 1 1  1];
mesh.conecCell{ 2, 1 } = [0 1 1  3];
mesh.conecCell{ 3, 1 } = [0 1 2  2];
mesh.conecCell{ 4, 1 } = [1 2 0  1 2];
mesh.conecCell{ 5, 1 } = [1 2 0  2 3];

initialConds                = struct();

analysisSettings = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.deltaT        =   2e-5;
analysisSettings.finalTime     =   1e-3;

analysisSettings.stopTolDeltau =   1e-8;
analysisSettings.stopTolForces =   1e-8;
analysisSettings.stopTolIts    =   15;

analysisSettings.posVariableLoadBC = 2;

otherParams = struct();
otherParams.problemName = 'static_plastic_von_mises_truss';
otherParams.plots_format = 'vtk';
otherParams.plots_deltaTs_separation = 2;

[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

deltas = -matUs(6 + 5, :)';
eles = sqrt(x2^2 + (z2 - deltas).^2);

valsLin = -2 * (E * (eles - L) ./ L) * A .* (z2 - deltas) ./ eles;

Etan = E * Kplas / (E + Kplas);
sigmas_hard = sigma_Y_0 + Etan * (abs(eles - L) ./ L - sigma_Y_0 / E);

valsFlu = 2 * (sigmas_hard * A) .* ((z2 - deltas) ./ eles);
valsP = min(valsLin, valsFlu);

% softening
Kplas = -E * .05;
materials.modelParams = [E Kplas sigma_Y_0];

analysisSettings.methodName    = 'arcLength';
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.2) / 100;
analysisSettings.incremArcLen = [2e-4 4e-5 * ones(1, 100)];

[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUsB, loadFactorsMatB, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

deltasB = -matUsB(6 + 5, :)';

eles = sqrt(x2^2 + (z2 - deltasB).^2);

valsLin = -2 * E * A * (z2 - deltasB) ./ eles .* (eles - L) ./ L;

Etan = E * Kplas / (E + Kplas);
sigmas_hard = sigma_Y_0 + Etan * (abs(eles - L) ./ L - sigma_Y_0 / E);

valsFlu = 2 * (sigmas_hard * A) .* ((z2 - deltasB) ./ eles);
valsPB = min(valsLin, valsFlu);

verifBoolean = (norm(valsP  - loadFactorsMat(:, 2)) < 1e-4 * norm(valsP)) && ...
               (norm(valsPB - loadFactorsMatB(:, 2)) < 1e-4 * norm(valsPB));

figure;
plot(deltas, valsP, 'b-x');
grid on;
hold on;
plot(deltas, loadFactorsMat(:, 2), 'r-o');
plot(deltasB, loadFactorsMatB(:, 2), 'g-o');
plot(deltasB, valsPB, 'k-*');
legend('analytic-hard', 'numeric-hard', 'analytic-soft', 'numeric-soft');
