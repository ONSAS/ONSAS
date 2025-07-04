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
% md# Cantilever problem using plate and shell elements
% md
clc
clear all
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes')
  clear all;
end
addpath(genpath([pwd '/../../src']));
% md
% md## Scalars
E = 200e6;
nu = 0.0;
Ly = .5;
Lx = 1;
tz = .05;
I = Ly * tz^3 / 12;
%
qx  = 100 ; % kN/m2
qz  = 200 ; % kN/m2
%
Fz = qz * Ly * tz ;
%
% ====================================================
% plate element
% ====================================================
%
% md
% md## Numerical solution using plate elements
% md
% md### Materials
% md
materials                    = struct();
materials(1).modelName  = 'elastic-linear';
materials(1).modelParams = [E nu];
% md
% md### Elements
% md
elements             = struct();
elements(1).elemType = 'edge';
elements(1).elemCrossSecParams = tz;
elements(2).elemType = 'triangle-plate';
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
boundaryConds(2).loadsBaseVals = [qx 0 0 0 -qz 0];
%
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
analysisSettings.deltaT        =   1;
analysisSettings.finalTime     =   1;
analysisSettings.stopTolDeltau =   1e-6;
analysisSettings.stopTolForces =   1e-6;
analysisSettings.stopTolIts    =   10;
% md
% md#### OtherParams
% md The nodalDispDamping is added into the model using:
otherParams                  = struct();
% md The name of the problem is:
% md
otherParams.problemName  = 'cantileverPlate-plateElem';
otherParams.plots_format = 'vtk';
% md
% md Execute ONSAS and save the results:
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% md
% md## verification
nelem = size(modelProperties.Conec, 1);
matSolic = getInternalForces(modelSolutions{end}.localInternalForces, 1:nelem, {'Mx', 'My', 'Mxy'});
numer_maxMx_plate = max(max(matSolic)) ;
numer_wmax_plate = min(matUs(5:6:end)) ;
% md
analy_maxMx = Fz * Lx / Ly ;
analy_wmax = -Fz * Lx^3 / (3 * E * I);
analy_dxmax = qx * Lx / E;
%
% ====================================================
% triangle CST element
% ====================================================
%
elements(2).elemType           = 'triangle';
elements(2).elemTypeParams     = 2;
elements(2).elemCrossSecParams = tz;
otherParams.problemName  = 'cantileverPlate-CSTElem';

[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
%
numer_dxmax_CST = max(matUs(1:6:end));
%
% ====================================================
% linear shell element
% ====================================================
%
elements(2).elemType           = 'triangle-shell';
elements(2).elemCrossSecParams = {'thickness', tz };
otherParams.problemName  = 'cantileverPlate-shell-linear';
otherParams.plots_format = 'vtk';
%
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
%
numer_dxmax_linear_shell = max(matUs(1:6:end));
numer_wmax_linear_shell  = min(matUs(5:6:end));
%
% ====================================================
% non linear shell element
% ====================================================
%
materials(1).modelName  = 'elastic-rotEngStr';
otherParams.problemName  = 'cantileverPlate-shell-nonlinear';
otherParams.plots_format = 'vtk';
%
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);

node = 10 ;
numer_dxmax_nonlin_shell = matUs(node*6-5,end);
numer_wmax_nonlin_shell  = matUs(node*6-1,end);


% ====================================================
% non linear frame element
% ====================================================
%
numElements = 5 ;
%
elements             = struct();
elements(1).elemType = 'node';
elements(2).elemType = 'frame';
%
elements(2).elemCrossSecParams{1, 1} = 'rectangle';
elements(2).elemCrossSecParams{2, 1} = [Ly tz];
%
boundaryConds                  = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];
%
Px = qx*Ly*tz ;
Pz = qz*Ly*tz ;
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) t;
boundaryConds(2).loadsBaseVals = [Px 0 0 0 -Pz 0];
%
initialConds = {};
%
mesh             = struct();
mesh.nodesCoords = [(0:(numElements))' * Lx / numElements  zeros(numElements + 1, 2)];
mesh.conecCell = { };
mesh.conecCell{ 1, 1 } = [0 1 1   1];
mesh.conecCell{ 2, 1 } = [0 1 2   numElements + 1];

for i = 1:numElements
  mesh.conecCell{ i + 2, 1 } = [1 2 0  i i + 1];
end
%
otherParams             = struct();
otherParams.problemName = 'cantileverPlate-NL-frame';
%
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

node = numElements+1;
dofs=node*6-5:node*6;

numer_dxmax_nonlin_frame = matUs(node*6-5,end);
numer_wmax_nonlin_frame = matUs(node*6-1,end);

[ numer_dxmax_linear_shell numer_dxmax_nonlin_shell numer_dxmax_nonlin_frame ];
[ numer_wmax_linear_shell numer_wmax_nonlin_shell numer_wmax_nonlin_frame ];

err_w = (numer_wmax_nonlin_shell-numer_wmax_nonlin_frame)/numer_wmax_nonlin_frame *100;
err_x = (numer_dxmax_nonlin_shell-numer_dxmax_nonlin_frame)/numer_dxmax_nonlin_frame *100;

analy_wmax;
analy_dxmax;

% md
verifBoolean = (abs(analy_wmax - numer_wmax_plate) / abs(analy_wmax))  < 5e-3 && ...
             (abs(analy_maxMx - numer_maxMx_plate) / abs(analy_maxMx)) < 1e-2 && ...
            (abs(analy_dxmax - numer_dxmax_CST) / abs(analy_dxmax)) < 1e-3 && ...
            (abs(analy_wmax  - numer_wmax_linear_shell) / abs(analy_wmax)) < 5e-3 && ...
            (abs(analy_dxmax - numer_dxmax_linear_shell) / abs(analy_dxmax)) < 1e-3 && ...
            (abs(analy_wmax  - numer_wmax_nonlin_shell) / abs(analy_wmax)) < 5e-3 && ...
            (abs(analy_wmax - numer_wmax_nonlin_frame) / abs(analy_wmax)) < 5e-3;


