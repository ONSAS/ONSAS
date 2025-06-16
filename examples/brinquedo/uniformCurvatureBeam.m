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
tz = .05;
qx  = 100 ; % kN/m^2
qz  = 10 ; % kN/m^2
% md
Ly = .5;
Lx = 1;

analy_maxMx = qz * Lx / 2;
qlin = qz * Ly;
I = Ly * tz^3 / 12;
analy_wmax = -qlin * Lx^4 / (8 * E * I);
analy_dxmax = qx * Lx / E;


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
analysisSettings.deltaT        =   1;
analysisSettings.finalTime     =   1;
analysisSettings.stopTolDeltau =   1e-6;
analysisSettings.stopTolForces =   1e-6;
analysisSettings.stopTolIts    =   10;
% md
% md#### OtherParams
% md The nodalDispDamping is added into the model using:
otherParams                  = struct();

% ===========================================================================================================================================
elements(2).elemType           = 'triangle-shell';
elements(2).elemTypeParams     = 2;
elements(2).elemCrossSecParams = {'thickness', tz };
otherParams.problemName  = 'cantileverPlate-shell-linear';
otherParams.plots_format = 'vtk';

[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% stop
numer_dxmax_linear_shell = max(matUs(1:6:end))
numer_wmax_linear_shell  = min(matUs(5:6:end))

% ===========================================================================================================================================
materials(1).modelName  = 'elastic-rotEngStr';
otherParams.problemName  = 'cantileverPlate-shell-nonlinear';
otherParams.plots_format = 'vtk';

[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
numer_dxmax_nonlin_shell = max(matUs(1:6:end))
numer_wmax_nonlin_shell  = min(matUs(5:6:end))

nod1 = 16;
nod2 = 38;
nod3 = 15; 
nodes_coords = mesh.nodesCoords([nod1;nod2;nod3],:) ;
Us_1 = matUs((nod1-1)*6+1:nod1*6,end);
Us_2 = matUs((nod2-1)*6+1:nod2*6,end);
Us_3 = matUs((nod3-1)*6+1:nod3*6,end);

[Us_1 Us_2 Us_3 ];
Us = [Us_1 ; Us_2 ; Us_3] ;

[ numer_dxmax_linear_shell numer_dxmax_nonlin_shell numer_wmax_linear_shell numer_wmax_nonlin_shell ]
analy_wmax
analy_dxmax



% geometrical scalar parameters
l = Lx;
ty = Ly;

Qz = 

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
boundaryConds(2).loadsTimeFact = @(t) t;
boundaryConds(2).loadsBaseVals = [0 0 0 -1 0 0];

initialConds = {};

mesh             = struct();
mesh.nodesCoords = [(0:(numElements))' * Lx / numElements  zeros(numElements + 1, 2)];
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


