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

close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes') % hidden
  clear all;
end % hidden
% add path
addpath(genpath([pwd '/../../src']));
% material scalar parameters
E = 200e9;
nu = 0.3;
% geometrical scalar parameters
l = 10;
ty = 1.0;
tz = .1;
Iy = ty * tz^3 / 12;
Mobj = E * Iy * 2 * pi / l ;
% the number of elements of the mesh
numElements = 10;
materialsL                 = struct();
materialsL.modelName  = 'elastic-linear';
materialsL.modelParams = [E nu];

materialsNL                 = struct();
materialsNL.modelName  = 'elastic-rotEngStr';
materialsNL.modelParams = [E nu];

elements             = struct();
elements(1).elemType = 'node';
elements(2).elemType = 'frame';
elements(2).elemCrossSecParams{1, 1} = 'rectangle';
elements(2).elemCrossSecParams{2, 1} = [ty tz];


boundaryConds                  = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];

boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) Mobj * t;
boundaryConds(2).loadsBaseVals = [0 0 0 -1 0 0];

initialConds = {};

mesh             = struct();
mesh.nodesCoords = [(0:(numElements))' * l / numElements  zeros(numElements + 1, 2)];
mesh.conecCell = { };
mesh.conecCell{ 1, 1 } = [0 1 1   1];
mesh.conecCell{ 2, 1 } = [0 1 2   numElements + 1];
for i = 1:numElements
  mesh.conecCell{ i + 2, 1 } = [1 2 0  i i + 1];
end


analysisSettings               = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.deltaT        =   0.0001;  % TEMPORARY
analysisSettings.finalTime      =   .0004;  % TEMPORARY
analysisSettings.stopTolDeltau =   1e-6;
analysisSettings.stopTolForces =   1e-6;
analysisSettings.stopTolIts    =   10;


otherParams             = struct();
otherParams.problemName = 'uniformCurvatureCantilever-frame';
otherParams.controlDofs = [numElements + 1  4];
otherParams.plots_format = 'vtk';

[modelCurrSol, modelProperties, BCsData] = initONSAS(materialsNL, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

angleControlDof      = (numElements + 1) * 6 - 2;
controlDispsNREngRot =  -matUs(angleControlDof, :);
loadFactorsNREngRot  =  loadFactorsMat(:, 2);
analyticLoadFactorsNREngRot = @(w) E * Iy * w / l;

verifBoolean = norm(analyticLoadFactorsNREngRot(controlDispsNREngRot) - loadFactorsNREngRot')                    < (norm(analyticLoadFactorsNREngRot(controlDispsNREngRot)) * 1e-4);
verifBoolean













% lw = 2.0;
% ms = 11;
% plotfontsize = 22;
% figure;
% plot(controlDispsNREngRot, analyticLoadFactorsNREngRot(controlDispsNREngRot), 'b-x', 'linewidth', lw, 'markersize', ms);
% hold on;
% grid on;
% plot(controlDispsNREngRot, loadFactorsNREngRot, 'k-o', 'linewidth', lw, 'markersize', ms);
% labx = xlabel('Displacement');
% laby = ylabel('\lambda');
% legend('analytic', 'NR-RotEng', 'location', 'North');
% set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize);
% set(labx, 'FontSize', plotfontsize);
% set(laby, 'FontSize', plotfontsize);
% print('output/verifCantileverBeam.png', '-dpng');


