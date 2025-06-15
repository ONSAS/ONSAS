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
qz  = 1 ; % kN/m^2
qx = 0 ;
Ly = 1.0 ; % m
Lx = 10.0 ; % m
% md## Numerical solution using plate elements
% md
% md### Materials
% md
materials                    = struct();
materials(1).modelName  = 'elastic-rotEngStr';
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
boundaryConds(3).loadsCoordSys = 'global';
boundaryConds(3).loadsTimeFact = @(t) t;
boundaryConds(3).loadsBaseVals = [qx 0 0 0 0 0];
% md
% md### mesh
% md
mesh = struct();
base_dir = '';
if strcmp(getenv('TESTS_RUN'), 'yes') && isfolder('examples')
  base_dir = ['.' filesep 'examples' filesep  'cantileverPlate' filesep];
end
[mesh.nodesCoords, mesh.conecCell] = meshFileReader([base_dir 'geometry_cantileverPlate.msh']);
assert(max(mesh.nodesCoords(:, 1)) == Lx && max(mesh.nodesCoords(:, 2)) == Ly);

% md### Initial conditions
initialConds                  = struct();
% md
% md#### Analysis settings
analysisSettings               = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.deltaT        =   0.5;
analysisSettings.finalTime     =   10;
analysisSettings.stopTolDeltau =   1e-6;
analysisSettings.stopTolForces =   1e-6;
analysisSettings.stopTolIts    =   10;

otherParams                  = struct();
otherParams.problemName  = 'cantileverPlate-shell-nonlinear';
otherParams.plots_format = 'vtk';

[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);

% nod = 18;
nod = 3;
controldofs = nod*6-1 ; % uz
controlDisps = -matUs(controldofs,:)
loadFactors  =  loadFactorsMat(:, 2)

lw = 2.0;
ms = 10;
plotfontsize = 10;
figure;
grid on;
plot(controlDisps, loadFactors, 'b-x', 'linewidth', lw, 'markersize', ms);
% hold on;

labx = xlabel('Displacement');
laby = ylabel('\lambda');
