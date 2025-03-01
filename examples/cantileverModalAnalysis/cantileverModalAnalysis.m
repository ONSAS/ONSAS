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
% md# Cantilever Modal Analysis
% md
% mdBefore defining the structs, the workspace is cleaned and the ONSAS directory is added to the path
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes')
  clear all;
end
addpath(genpath([pwd '/../../src']));
% md
% mdThe material scalar parameters are set.
E = 200e9;
nu = 0.3;
rho = 700;
% mdThe cross-section of the beam is rectangular. The widths and other geometry scalar parameters are computed.
l = 10;
diam = .01;
numElements = 16; % Number of elements

tf     = 0.3; % s
deltat = 0.1; % s

% md### materials
materials             = struct();
materials.modelName   = 'elastic-rotEngStr'; % 'elastic-linear';%
materials.modelParams = [E nu];
materials.density     = rho;
% md
% md### elements
elements = struct();
elements(1).elemType = 'node';
elements(2).elemType = 'frame';
% elements(2).elemCrossSecParams = { 'rectangle' , [ty tz] } ;
elements(2).elemCrossSecParams = { 'circle', diam };
% md The consistent mass approach is considered for the dynamic analysis
elements(2).massMatType = 'consistent';
% md
boundaryConds                  = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];
% md and the second corresponds to a time dependant external force
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) t;
boundaryConds(2).loadsBaseVals = [0 0 1 0 0 0];
% md
% md### initial Conditions
% md homogeneous initial conditions are considered, then an empty struct is set:
initialConds = struct();
% md
mesh = struct();
mesh.nodesCoords = [(0:(numElements))' * l / numElements  zeros(numElements + 1, 2)];

mesh.conecCell = { };
mesh.conecCell{ 1, 1 } = [0 1 1  1];
mesh.conecCell{ 2, 1 } = [0 1 2  numElements + 1];

for i = 1:numElements
  mesh.conecCell{ i + 2, 1 } = [1 2 0  i i + 1];
end

% md### analysisSettings
analysisSettings = struct();
analysisSettings.methodName    = 'newmark';
analysisSettings.deltaT        =   deltat;
analysisSettings.finalTime     =   tf;
analysisSettings.stopTolDeltau =   1e-10;
analysisSettings.stopTolForces =   1e-8;
analysisSettings.stopTolIts    =   20;
analysisSettings.modalAnalysisBoolean = true;
% md
% md## otherParams
otherParams = struct();
otherParams.problemName = 'cantilever_modal_analysis';
otherParams.exportFirstMatrices = true;
% md ONSAS execution
% mdFirst the input structs are converted to structs with the model information
[modelCurrSol, modelProperties, BCsData] = ONSAS_init(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[coRotMatUs, loadFactorsMat, modelSolutions] = ONSAS_solve(modelCurrSol, modelProperties, BCsData);
% md
% md the report is generated
outputReport(modelProperties.outputDir, modelProperties.problemName);

% md
% load matrices file
load([pwd filesep 'output' filesep 'matrices.mat']);

KTred = KT(neumdofs, neumdofs);
Mred  = massMat(neumdofs, neumdofs);

%%
s = linspace(0, l, numElements + 1)';
z = linspace(0, l, numElements * 10 + 1)';
[a, b] = eig(full(KTred), full(Mred));
eigenvalues = flipud(diag(b));

numer_modes = [zeros(6, (numElements) * 6); fliplr(a)];

% Analytical solution
betavect = [1.8751 4.6941 7.8547 10.9955 14.13717];
sigvect  = [0.73409 1.01846 0.9992 1.00003 1];
modes = [1 2 3 4 5];
modes_errors = zeros(size(modes));
close all;
for mode = modes
  figure(mode);
  num_disp_mode = numer_modes(5:6:end, mode);
  num_disp_mode = num_disp_mode ./ max(abs(num_disp_mode)) * sign(num_disp_mode(end));
  plot(s, num_disp_mode, 'rx', 'linewidth', 1.2, 'markersize', 12);
  hold on;
  grid on;
  title(['mode ' num2str(mode)]);
  % freqStruct = (betavect.^2)*sqrt(E*I/(ms+ma)/l^4)/(2*pi);
  S1 = cosh(betavect(mode) * z / l) - cos(betavect(mode) * z / l) - sigvect(mode) * sinh(betavect(mode) * z / l) + sigvect(mode) * sin(betavect(mode) * z / l);
  S1 = S1 ./ max(abs(S1)) * sign(S1(end));
  spanplot = 10;
  S1_points = S1(1:spanplot:end);
  modes_errors(mode) = norm(S1_points - num_disp_mode) / norm(S1_points);
  plot(z(1:spanplot:end), S1_points, 'bo', 'linewidth', 1.0, 'markersize', 11);
  plot(z, S1, 'b', 'linewidth', 1.0, 'markersize', 11);
  legend('ONSAS', 'analytical');
end

verifBoolean = (modes_errors(1) < 1e-3)  &&  (modes_errors(2) < 1e-3);
