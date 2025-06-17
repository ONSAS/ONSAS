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
% md# Cantilever problem using plate and shell elements
% md
clc
clear all
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes')
  clear all;
end
addpath(genpath([pwd '/../../src']));


E = 1.2e7 ; % kN/m2
nu = 0.3;
tz = 0.1;
Ly = 1 ; % m
Lx = 10 ; % m
qz  = 1/(Ly*tz) ; % kN/m^2
qx = 0 ;

% md## Numerical solution using plate elements
% md
% md### Materials
% md
materials                    = struct();

materials(1).modelParams = [E nu];
% md
% md### Elements
% md
elements             = struct();
elements(1).elemType = 'edge';
elements(1).elemCrossSecParams = tz;

elements(2).elemTypeParams     = 2;
elements(2).elemType           = 'triangle-shell';
elements(2).elemCrossSecParams = {'thickness', tz };
% md
% md### Boundary conditions
% md
boundaryConds                  = struct();
boundaryConds(1).imposDispDofs =  [1 2 3 4 5 6];
boundaryConds(1).imposDispVals =  [0 0 0 0 0 0];
%
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) t;
boundaryConds(2).loadsBaseVals = [0 0 0 0 -qz 0];
%
% md
% md### mesh
% md
mesh = struct();
base_dir = '';
[mesh.nodesCoords, mesh.conecCell] = meshFileReader([base_dir 'cantilever.msh']);
assert(max(mesh.nodesCoords(:, 1)) == Lx && max(mesh.nodesCoords(:, 2)) == Ly);

% md### Initial conditions
initialConds                  = struct();
% md
% md#### Analysis settings
analysisSettings               = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.deltaT        =   1;
analysisSettings.finalTime     =   20;
analysisSettings.stopTolDeltau =   1e-6;
analysisSettings.stopTolForces =   1e-6;
analysisSettings.stopTolIts    =   10;

otherParams                  = struct();

materials(1).modelName  = 'elastic-rotEngStr';
otherParams.problemName  = 'cantileverPlate_nonLinear';
otherParams.plots_format = 'vtk';

[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);

% nod = 18;
nod = 3;
controldofs_z = nod*6-1 ; % uz
controldofs_x = nod*6-5 ; % ux
controlDisps_z_Shell = -matUs(controldofs_z,:)
controlDisps_x_Shell = matUs(controldofs_x,:)
loadFactors_Shell  =  loadFactorsMat(:, 2)

Iy = Ly*tz^3/12 ;
f_ana = qz*loadFactors_Shell*Ly*tz*Lx^3 / (3*E*Iy) 

% =============================================================================================================================================
stop

% geometrical scalar parameters
l = Lx;
ty = Ly;

Qz = qz*Ly*tz ;

% the number of elements of the mesh
numElements = 10;

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
boundaryConds(2).loadsTimeFact = @(t) t;
boundaryConds(2).loadsBaseVals = [0 0 0 0 -Qz 0];

initialConds = {};

mesh             = struct();
mesh.nodesCoords = [(0:(numElements))' * Lx / numElements  zeros(numElements + 1, 2)];
mesh.conecCell = { };
mesh.conecCell{ 1, 1 } = [0 1 1   1];
mesh.conecCell{ 2, 1 } = [0 1 2   numElements + 1];
for i = 1:numElements
  mesh.conecCell{ i + 2, 1 } = [1 2 0  i i + 1];
end

otherParams             = struct();
otherParams.problemName = 'cantilever_frame';
otherParams.plots_format = 'vtk';

[modelCurrSol, modelProperties, BCsData] = initONSAS(materialsNL, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

nod = numElements+1;
controldofs_z = nod*6-1 ; % uz
controldofs_x = nod*6-5 ; % ux
controlDisps_z_Frame = -matUs(controldofs_z,:)
controlDisps_x_Frame = matUs(controldofs_x,:)
loadFactors_Frame  =  loadFactorsMat(:, 2)

lw = 2.0;
ms = 10;
plotfontsize = 10;
figure;

plot(controlDisps_z_Shell, loadFactors_Shell, 'b-x', 'linewidth', lw, 'markersize', ms);
hold on;
grid on;
plot(controlDisps_x_Shell, loadFactors_Shell, 'g-o', 'linewidth', lw, 'markersize', ms);
plot(controlDisps_z_Frame, loadFactors_Shell, 'k--x', 'linewidth', lw, 'markersize', ms);
plot(controlDisps_x_Frame, loadFactors_Shell, 'r--*', 'linewidth', lw, 'markersize', ms);

labx = xlabel('Displacement');
laby = ylabel('\lambda');