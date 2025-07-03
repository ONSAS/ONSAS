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
% EulerColumn
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes')
  clear all;
end
% add path
addpath(genpath([pwd '/../../src']));
% material scalar parameters
E = 30e6; % kN/m2
nu = 0.2;
nu = 0 ;
% geometrical scalar parameters
l = 5; % m
ty = .1; % m
tz = .1; % m
% the number of elements of the mesh
numElements = 8;

% MEB parameters
% Materials
% ----------------------------------------------------------------------
materials  = struct();
materials.modelName  = 'elastic-rotEngStr';
materials.modelParams = [E nu];
% Elements
% ----------------------------------------------------------------------
% Types
elements  = struct();
elements(1).elemType = 'node';
elements(2).elemType = 'frame';
% Sections
elements(2).elemCrossSecParams = { 'rectangle', [ty tz] };
elements(2).massMatType   = 'consistent';

% Boundary conditions
% ----------------------------------------------------------------------
% Pinned support
boundaryConds  = struct();
boundaryConds(1).imposDispDofs = [1 3 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0];
% Roller support
boundaryConds(2).imposDispDofs = [1 3 6];
boundaryConds(2).imposDispVals = [0 0 0];
% Load
P = -1;
imp = P / 1e4;
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsBaseVals = [0 0 0 imp P 0];
boundaryConds(2).loadsTimeFact = @(t) t;

% Initial conditions
% ----------------------------------------------------------------------
initialConds                = struct();

% Mesh
% Nodes coords
% ----------------------------------------------------------------------
mesh  = struct();
mesh.nodesCoords = [zeros(numElements + 1, 2) (0:(numElements))' * l / numElements];
% Conec cell
% ----------------------------------------------------------------------
mesh.conecCell = { };
mesh.conecCell{ 1, 1 } = [0 1 1 1];
mesh.conecCell{ 2, 1 } = [0 1 2  numElements + 1];
for i = 1:numElements
  mesh.conecCell{ i + 2, 1 } = [1 2 0  i i + 1];
end
% Analysis settings

% Parameters
% ----------------------------------------------------------------------
analysisSettings  = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.deltaT        =   2.5;
analysisSettings.finalTime     =   140;
analysisSettings.stopTolDeltau =   1e-8;
analysisSettings.stopTolForces =   1e-8;
analysisSettings.stopTolIts    =   20;
% md
otherParams  = struct();
otherParams.problemName = 'EulerColumn';
otherParams.controlDofs = [numElements + 1  5];
otherParams.plots_format = 'vtk';
% otherParams.plots_format = 'vtk' ;

A = ty * tz;
I = max(ty, tz) * min(ty, tz)^3 / 12;
Pcrit = pi()^2 * E * I / l^2;

[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
% md

Dof          = (numElements / 2 + 1) * 6 - 5;
controlDisps =  matUs(Dof, :);
loadFactors  =  loadFactorsMat(:, 2);

lw = 2.0;
ms = 11;
plotfontsize = 22;
figure;
grid on;
plot(controlDisps, loadFactors, 'k-o', 'linewidth', lw, 'markersize', ms);
labx = xlabel('Displacement');
laby = ylabel('$\lambda$');

verifBoolean = (abs(controlDisps(end) - 1.7406) / 1.7406) < 1e-3;
% stop


% ====================================================
% shell non linear
% ====================================================
%
% elements             = struct();
% elements(1).elemType = 'edge';
% elements(1).elemCrossSecParams = tz;
% elements(2).elemType = 'triangle-shell';
% elements(2).elemCrossSecParams = {'thickness', tz };
% %
% boundaryConds                  = struct();
% boundaryConds(1).imposDispDofs =  [1 3 5];
% boundaryConds(1).imposDispVals =  [0 0 0];
% %
% boundaryConds(2).imposDispDofs =  [1 3];
% boundaryConds(2).imposDispVals =  [0 0];
% %
% boundaryConds(2).loadsCoordSys = 'global';
% boundaryConds(2).loadsTimeFact = @(t)  t;
% imp = P/1e6
% boundaryConds(2).loadsBaseVals = [0 0 0 imp/(ty*tz) P/(ty*tz) 0];
% %
% mesh = struct();
% base_dir = '';
% if strcmp(getenv('TESTS_RUN'), 'yes') && isfolder('examples')
%   base_dir = ['.' filesep 'examples' filesep  'eulerColumn' filesep];
% end
% [mesh.nodesCoords, mesh.conecCell] = meshFileReader([base_dir 'geometry_column.msh']);
% assert(max(mesh.nodesCoords(:, 3)) == l && max(mesh.nodesCoords(:, 2)) == ty);
% %
% otherParams.problemName = 'EulerColumn-NL-shell';
% otherParams.plots_format = 'vtk';
% %
% analysisSettings               = struct() ;
% % analysisSettings.methodName    = 'newtonRaphson';
% % analysisSettings.deltaT        =   1.25 ;
% % analysisSettings.finalTime     =   140 ;  


% % analysisSettings.methodName   = 'arcLength';
% finalTime = 20;
% analysisSettings.finalTime    = finalTime;
% analysisSettings.deltaT        =   2.5 ;
% analysisSettings.stopTolDeltau =   1e-6 ;
% analysisSettings.stopTolForces =   1e-6 ;
% analysisSettings.stopTolIts    =   15 ;
% % analysisSettings.incremArcLen = [ 0.0001*ones(1, 40) 0.0001/2*ones(1, 40) ];
% % analysisSettings.incremArcLen = 0.00015;
% analysisSettings.iniDeltaLamb = 1 / 100;
% analysisSettings.posVariableLoadBC = 2;


% analysisSettings.incremArcLen = [ ones(1, 19) ];
% % md Sets arcLengthFlag = 2 to secifiy Jirasek constraint method.
% analysisSettings.ALdominantDOF = [73*6-5 1];

% %
% [modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

% node = 73 ;
% Dof_shell         = node * 6 - 5;
% controlDispsShell      = matUs(Dof_shell, :);
% loadFactorsShell  = loadFactorsMat(:, 2)

% close all;
% lw = 2.0;
% ms = 11;
% plotfontsize = 22;
% figure;
% grid on;
% plot(controlDispsShell, loadFactorsShell, 'k-o', 'linewidth', lw, 'markersize', ms);
% labx = xlabel('Displacement');
% laby = ylabel('$\lambda$');
