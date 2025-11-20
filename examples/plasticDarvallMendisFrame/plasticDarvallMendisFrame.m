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

%  Darvall-Mendis Frame Analysis / Softening Hinges
close all;
clear all;

addpath(genpath([pwd '/../src/'])); % add ONSAS directory to path
addpath(genpath([pwd '/../../ONSAS/src'])); % add ONSAS directory to path

global TZERO

% assumed XY plane

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

% geometry
L  = 3.048;            % m
L1 = 0.55 * L;

Inertia = 0.001;       % m^4

A = 0.103;            % m^2

% material
E = 20.68e6;          % kN/m^2 [kPa]                                                                                                                                                % kN/m^2 [kPa]

nu = 0.3;     % Poisson's ratio

a = -1e-12;
% a = -0.04 ;
% a = -0.06 ;
% a = -0.0718 ;

kh1 = 30000;               % kN.m^2
kh2 = 3000;

Ks = a * E * Inertia / L1 * 10;    % kN.m

Mc_c = 158.18;             % kN.m
My_c = 158.18;
Mu_c = 158.18;

Mc_b = 169.48;             % kN.m
My_b = 169.48;
Mu_b = 169.48;

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

materials = struct();
materials(1).modelName = 'plastic-2Dframe';
materials(1).modelParams = [E Mc_c My_c Mu_c kh1 kh2 Ks nu];
materials(2).modelName = 'plastic-2Dframe';
materials(2).modelParams = [E Mc_b My_b Mu_b kh1 kh2 Ks nu];
% elements
elements = struct();
elements(1).elemType  = 'node';
elements(2).elemType  = 'frame';
elements(2).elemCrossSecParams = {'generic'; [A 1 Inertia Inertia] };
% Boundary Conditions
% Supports
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];

boundaryConds(3).imposDispDofs = [2 4 5];
boundaryConds(3).imposDispVals = [0 0 0];

% Loads
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsBaseVals = [0 0 -1 0 0 0];
boundaryConds(2).loadsTimeFact = @(t) t;

boundaryConds(2).imposDispDofs = [2 4 5];
boundaryConds(2).imposDispVals = [0 0 0];

% Mesh
% Mesh nodes
mesh = struct();
mesh.nodesCoords = [0      0       0; ...
                    0      L / 2     0; ...
                    0      L       0; ...
                    L1 / 2   L       0; ...
                    L1     L       0; ... % loaded node
                    (L - L1) / 2 + L1     L       0; ...
                    L      L       0; ...
                    L      L / 2     0; ...
                    L      0       0];
% Conec Cell
mesh.conecCell = { };
% nodes
mesh.conecCell{1, 1 } = [0 1 1   1];
mesh.conecCell{9, 1 } = [0 1 1   9];

mesh.conecCell{2, 1 } = [0 1 3   2];
mesh.conecCell{3, 1 } = [0 1 3   3];
mesh.conecCell{4, 1 } = [0 1 3   4];
mesh.conecCell{5, 1 } = [0 1 2   5]; % loaded node
mesh.conecCell{6, 1 } = [0 1 3   6];
mesh.conecCell{7, 1 } = [0 1 3   7];
mesh.conecCell{8, 1 } = [0 1 3   8];

% and frame elements
mesh.conecCell{10, 1 }  = [1 2 0   1 2];
mesh.conecCell{11, 1 }  = [1 2 0   2 3];
mesh.conecCell{12, 1 }  = [2 2 0   3 4];
mesh.conecCell{13, 1 }  = [2 2 0   4 5];
mesh.conecCell{14, 1 }  = [2 2 0   5 6];
mesh.conecCell{15, 1 }  = [2 2 0   6 7];
mesh.conecCell{16, 1 }  = [1 2 0   7 8];
mesh.conecCell{17, 1 }  = [1 2 0   8 9];

% InitialConditions
% empty struct
initialConds = struct();

% Analysis settings
analysisSettings                    = {};
analysisSettings.methodName         = 'arcLength';
analysisSettings.deltaT             = 1;
# analysisSettings.incremArcLen       = [1e-5 * ones(1, 2100)];
analysisSettings.incremArcLen       = [1e-5 * ones(1, 2)];
analysisSettings.finalTime          = length(analysisSettings.incremArcLen);
analysisSettings.iniDeltaLamb       = 1;
analysisSettings.posVariableLoadBC  = 2;
analysisSettings.stopTolDeltau      = 1e-14;
analysisSettings.stopTolForces      = 1e-8;
analysisSettings.stopTolIts         = 30;
analysisSettings.ALdominantDOF      = [4 * 6 + 3 -1];

%
otherParams = struct();
otherParams.problemName = 'plasticFrameDarvallMendis';
% otherParams.plots_format = 'vtk' ;

[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);

[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);

displacements_perfect = matUs((4) * 6 + 3, :) * 100; % node with vertical load applied
loadfactors_perfect = loadFactorsMat(:, 2);

Nelem = 8;

moments_hist = cell(Nelem, 6, length(modelSolutions));
for ii = 1:Nelem
  for i = 1:length(modelSolutions)
    aux = modelSolutions{i}.localInternalForces(ii);
    moments_hist(ii, :, i) = { aux.Mz; aux.Mz2; aux.Mz3; aux.tM; aux.Mz_integrado_izq; aux.Mz_integrado_der };
  end
end
Mn1_numericONSAS = moments_hist(:, 1, :);
Mn2_numericONSAS = moments_hist(:, 2, :);
Mn3_numericONSAS = moments_hist(:, 3, :);
tMn_numericONSAS = moments_hist(:, 4, :);
Mflizq_numericONSAS = moments_hist(:, 5, :);
Mflder_numericONSAS = moments_hist(:, 6, :);

fprintf('\n');
fprintf('Number of softening hinges / (%d elements)', Nelem);
fprintf('\n');
fprintf('[');
fprintf('%g, ', TZERO(1:end - 1));
fprintf('%g]\n', TZERO(end));

matUs


stop

fprintf(' a = 0\n');
fprintf(' First hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4) * 6 + 3, TZERO(4)) * 100, loadFactorsMat(TZERO(4), 2));
fprintf(' Second hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4) * 6 + 3, TZERO(7)) * 100, loadFactorsMat(TZERO(7), 2));
fprintf(' Third hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4) * 6 + 3, TZERO(2)) * 100, loadFactorsMat(TZERO(2), 2));

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

% geometry
L  = 3.048;            % m
L1 = 0.55 * L;

Inertia = 0.001;       % m^4

A = 0.103;            % m^2

% material
E = 20.68e6;          % kN/m^2 [kPa]                                                                                                                                                % kN/m^2 [kPa]

nu = 0.3;     % Poisson's ratio

a = -0.04;
% a = -0.06 ;
% a = -0.0718 ;

kh1 = 30000;               % kN.m^2
kh2 = 3000;

Ks = a * E * Inertia / L1 * 10;    % kN.m

Mc_c = 158.18;             % kN.m
My_c = 158.18;
Mu_c = 158.18;

Mc_b = 169.48;             % kN.m
My_b = 169.48;
Mu_b = 169.48;

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

materials = struct();
materials(1).modelName = 'plastic-2Dframe';
materials(1).modelParams = [E Mc_c My_c Mu_c kh1 kh2 Ks nu];
materials(2).modelName = 'plastic-2Dframe';
materials(2).modelParams = [E Mc_b My_b Mu_b kh1 kh2 Ks nu];
% elements
elements = struct();
elements(1).elemType  = 'node';
elements(2).elemType  = 'frame';
elements(2).elemCrossSecParams = {'generic'; [A 1 Inertia Inertia] };
% Boundary Conditions
% Supports
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];

boundaryConds(3).imposDispDofs = [2 4 5];
boundaryConds(3).imposDispVals = [0 0 0];

% Loads
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsBaseVals = [0 0 -1 0 0 0];
boundaryConds(2).loadsTimeFact = @(t) t;

boundaryConds(2).imposDispDofs = [2 4 5];
boundaryConds(2).imposDispVals = [0 0 0];

% Mesh
% Mesh nodes
mesh = struct();
mesh.nodesCoords = [0      0       0; ...
                    0      L / 2     0; ...
                    0      L       0; ...
                    L1 / 2   L       0; ...
                    L1     L       0; ... % loaded node
                    (L - L1) / 2 + L1     L       0; ...
                    L      L       0; ...
                    L      L / 2     0; ...
                    L      0       0];
% Conec Cell
mesh.conecCell = { };
% nodes
mesh.conecCell{1, 1 } = [0 1 1   1];
mesh.conecCell{9, 1 } = [0 1 1   9];

mesh.conecCell{2, 1 } = [0 1 3   2];
mesh.conecCell{3, 1 } = [0 1 3   3];
mesh.conecCell{4, 1 } = [0 1 3   4];
mesh.conecCell{5, 1 } = [0 1 2   5]; % loaded node
mesh.conecCell{6, 1 } = [0 1 3   6];
mesh.conecCell{7, 1 } = [0 1 3   7];
mesh.conecCell{8, 1 } = [0 1 3   8];

% and frame elements
mesh.conecCell{10, 1 }  = [1 2 0   1 2];
mesh.conecCell{11, 1 }  = [1 2 0   2 3];
mesh.conecCell{12, 1 }  = [2 2 0   3 4];
mesh.conecCell{13, 1 }  = [2 2 0   4 5];
mesh.conecCell{14, 1 }  = [2 2 0   5 6];
mesh.conecCell{15, 1 }  = [2 2 0   6 7];
mesh.conecCell{16, 1 }  = [1 2 0   7 8];
mesh.conecCell{17, 1 }  = [1 2 0   8 9];

% InitialConditions
% empty struct
initialConds = struct();

% Analysis settings
analysisSettings                    = {};
analysisSettings.methodName         = 'arcLength';
analysisSettings.deltaT             = 1;
analysisSettings.incremArcLen       = [1e-5 * ones(1, 2500)];
analysisSettings.finalTime          = length(analysisSettings.incremArcLen);
analysisSettings.iniDeltaLamb       = 1;
analysisSettings.posVariableLoadBC  = 2;
analysisSettings.stopTolDeltau      = 1e-14;
analysisSettings.stopTolForces      = 1e-8;
analysisSettings.stopTolIts         = 30;
analysisSettings.ALdominantDOF      = [4 * 6 + 3 -1];

%
otherParams = struct();
otherParams.problemName = 'plastic_2dframe';
% otherParams.plots_format = 'vtk' ;

[modelCurrSol, modelProperties, BCsData] = ONSAS_init(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);

[matUs, loadFactorsMat, modelSolutions] = ONSAS_solve(modelCurrSol, modelProperties, BCsData);

displacements = matUs((4) * 6 + 3, :) * 100; % node with vertical load applied
loadfactors = loadFactorsMat(:, 2);

Nelem = 8;

moments_hist = cell(Nelem, 6, length(modelSolutions));
for ii = 1:Nelem
  for i = 1:length(modelSolutions)
    aux = modelSolutions{i}.localInternalForces(ii);
    moments_hist(ii, :, i) = { aux.Mz; aux.Mz2; aux.Mz3; aux.tM; aux.Mz_integrado_izq; aux.Mz_integrado_der };
  end
end
Mn1_numericONSAS = moments_hist(:, 1, :);
Mn2_numericONSAS = moments_hist(:, 2, :);
Mn3_numericONSAS = moments_hist(:, 3, :);
tMn_numericONSAS = moments_hist(:, 4, :);
Mflizq_numericONSAS = moments_hist(:, 5, :);
Mflder_numericONSAS = moments_hist(:, 6, :);

fprintf('\n');
fprintf('Number of softening hinges / (%d elements)', Nelem);
fprintf('\n');
fprintf('[');
fprintf('%g, ', TZERO(1:end - 1));
fprintf('%g]\n', TZERO(end));

fprintf(' a = -0.04\n');
fprintf(' First hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4) * 6 + 3, TZERO(4)) * 100, loadFactorsMat(TZERO(4), 2));
fprintf(' Second hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4) * 6 + 3, TZERO(7)) * 100, loadFactorsMat(TZERO(7), 2));

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

% a = -0.04 ;
a = -0.06;
% a = -0.0718 ;

kh1 = 30000;               % kN.m^2
kh2 = 3000;

Ks = a * E * Inertia / L1 * 10;    % kN.m

Mc_c = 158.18;             % kN.m
My_c = 158.18;
Mu_c = 158.18;

Mc_b = 169.48;             % kN.m
My_b = 169.48;
Mu_b = 169.48;

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

materials = struct();
materials(1).modelName = 'plastic-2Dframe';
materials(1).modelParams = [E Mc_c My_c Mu_c kh1 kh2 Ks nu];
materials(2).modelName = 'plastic-2Dframe';
materials(2).modelParams = [E Mc_b My_b Mu_b kh1 kh2 Ks nu];
% elements
elements = struct();
elements(1).elemType  = 'node';
elements(2).elemType  = 'frame';
elements(2).elemCrossSecParams = {'generic'; [A 1 Inertia Inertia] };
% Boundary Conditions
% Supports
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];

boundaryConds(3).imposDispDofs = [2 4 5];
boundaryConds(3).imposDispVals = [0 0 0];

% Loads
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsBaseVals = [0 0 -1 0 0 0];
boundaryConds(2).loadsTimeFact = @(t) t;

boundaryConds(2).imposDispDofs = [2 4 5];
boundaryConds(2).imposDispVals = [0 0 0];

% Mesh
% Mesh nodes
mesh = struct();
mesh.nodesCoords = [0      0       0; ...
                    0      L / 2     0; ...
                    0      L       0; ...
                    L1 / 2   L       0; ...
                    L1     L       0; ... % loaded node
                    (L - L1) / 2 + L1     L       0; ...
                    L      L       0; ...
                    L      L / 2     0; ...
                    L      0       0];
% Conec Cell
mesh.conecCell = { };
% nodes
mesh.conecCell{1, 1 } = [0 1 1   1];
mesh.conecCell{9, 1 } = [0 1 1   9];

mesh.conecCell{2, 1 } = [0 1 3   2];
mesh.conecCell{3, 1 } = [0 1 3   3];
mesh.conecCell{4, 1 } = [0 1 3   4];
mesh.conecCell{5, 1 } = [0 1 2   5]; % loaded node
mesh.conecCell{6, 1 } = [0 1 3   6];
mesh.conecCell{7, 1 } = [0 1 3   7];
mesh.conecCell{8, 1 } = [0 1 3   8];

% and frame elements
mesh.conecCell{10, 1 }  = [1 2 0   1 2];
mesh.conecCell{11, 1 }  = [1 2 0   2 3];
mesh.conecCell{12, 1 }  = [2 2 0   3 4];
mesh.conecCell{13, 1 }  = [2 2 0   4 5];
mesh.conecCell{14, 1 }  = [2 2 0   5 6];
mesh.conecCell{15, 1 }  = [2 2 0   6 7];
mesh.conecCell{16, 1 }  = [1 2 0   7 8];
mesh.conecCell{17, 1 }  = [1 2 0   8 9];

% InitialConditions
% empty struct
initialConds = struct();

% Analysis settings
analysisSettings                    = {};
analysisSettings.methodName         = 'arcLength';
analysisSettings.deltaT             = 1;
analysisSettings.incremArcLen       = [1e-5 * ones(1, 2500)];
analysisSettings.finalTime          = length(analysisSettings.incremArcLen);
analysisSettings.iniDeltaLamb       = 1;
analysisSettings.posVariableLoadBC  = 2;
analysisSettings.stopTolDeltau      = 1e-14;
analysisSettings.stopTolForces      = 1e-8;
analysisSettings.stopTolIts         = 30;
analysisSettings.ALdominantDOF      = [4 * 6 + 3 -1];

%
otherParams = struct();
otherParams.problemName = 'plastic_2dframe';
% otherParams.plots_format = 'vtk' ;

[modelCurrSol, modelProperties, BCsData] = ONSAS_init(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);

[matUs, loadFactorsMat, modelSolutions] = ONSAS_solve(modelCurrSol, modelProperties, BCsData);

displacementsII = matUs((4) * 6 + 3, :) * 100; % node with vertical load applied
loadfactorsII = loadFactorsMat(:, 2);

% print('-f1','../tex/imagenes/darvallmendisframeplot.png','-dpng') ;
fprintf(' a = -0.06\n');
fprintf(' First hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4) * 6 + 3, TZERO(4)) * 100, loadFactorsMat(TZERO(4), 2));
fprintf(' Second hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4) * 6 + 3, TZERO(7)) * 100, loadFactorsMat(TZERO(7), 2));
% fprintf(' Third hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4)*6+3,TZERO(2))*100,loadFactorsMat(TZERO(2),2)) ;

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

% a = -0.04 ;
% a = -0.06 ;
a = -0.0718;

kh1 = 30000;               % kN.m^2
kh2 = 3000;

Ks = a * E * Inertia / L1 * 10;    % kN.m

Mc_c = 158.18;             % kN.m
My_c = 158.18;
Mu_c = 158.18;

Mc_b = 169.48;             % kN.m
My_b = 169.48;
Mu_b = 169.48;

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

materials = struct();
materials(1).modelName = 'plastic-2Dframe';
materials(1).modelParams = [E Mc_c My_c Mu_c kh1 kh2 Ks nu];
materials(2).modelName = 'plastic-2Dframe';
materials(2).modelParams = [E Mc_b My_b Mu_b kh1 kh2 Ks nu];
% elements
elements = struct();
elements(1).elemType  = 'node';
elements(2).elemType  = 'frame';
elements(2).elemCrossSecParams = {'generic'; [A 1 Inertia Inertia] };
% Boundary Conditions
% Supports
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];

boundaryConds(3).imposDispDofs = [2 4 5];
boundaryConds(3).imposDispVals = [0 0 0];

% Loads
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsBaseVals = [0 0 -1 0 0 0];
boundaryConds(2).loadsTimeFact = @(t) t;

boundaryConds(2).imposDispDofs = [2 4 5];
boundaryConds(2).imposDispVals = [0 0 0];

% Mesh
% Mesh nodes
mesh = struct();
mesh.nodesCoords = [0      0       0; ...
                    0      L / 2     0; ...
                    0      L       0; ...
                    L1 / 2   L       0; ...
                    L1     L       0; ... % loaded node
                    (L - L1) / 2 + L1     L       0; ...
                    L      L       0; ...
                    L      L / 2     0; ...
                    L      0       0];
% Conec Cell
mesh.conecCell = { };
% nodes
mesh.conecCell{1, 1 } = [0 1 1   1];
mesh.conecCell{9, 1 } = [0 1 1   9];

mesh.conecCell{2, 1 } = [0 1 3   2];
mesh.conecCell{3, 1 } = [0 1 3   3];
mesh.conecCell{4, 1 } = [0 1 3   4];
mesh.conecCell{5, 1 } = [0 1 2   5]; % loaded node
mesh.conecCell{6, 1 } = [0 1 3   6];
mesh.conecCell{7, 1 } = [0 1 3   7];
mesh.conecCell{8, 1 } = [0 1 3   8];

% and frame elements
mesh.conecCell{10, 1 }  = [1 2 0   1 2];
mesh.conecCell{11, 1 }  = [1 2 0   2 3];
mesh.conecCell{12, 1 }  = [2 2 0   3 4];
mesh.conecCell{13, 1 }  = [2 2 0   4 5];
mesh.conecCell{14, 1 }  = [2 2 0   5 6];
mesh.conecCell{15, 1 }  = [2 2 0   6 7];
mesh.conecCell{16, 1 }  = [1 2 0   7 8];
mesh.conecCell{17, 1 }  = [1 2 0   8 9];

% InitialConditions
% empty struct
initialConds = struct();

% Analysis settings
analysisSettings                    = {};
analysisSettings.methodName         = 'arcLength';
analysisSettings.deltaT             = 1;
analysisSettings.incremArcLen       = [1e-5 * ones(1, 1580)];
analysisSettings.finalTime          = length(analysisSettings.incremArcLen);
analysisSettings.iniDeltaLamb       = 1;
analysisSettings.posVariableLoadBC  = 2;
analysisSettings.stopTolDeltau      = 1e-14;
analysisSettings.stopTolForces      = 1e-8;
analysisSettings.stopTolIts         = 30;
analysisSettings.ALdominantDOF      = [4 * 6 + 3 -1];

%
otherParams = struct();
otherParams.problemName = 'plastic_2dframe';
% otherParams.plots_format = 'vtk' ;

[modelCurrSol, modelProperties, BCsData] = ONSAS_init(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);

[matUs, loadFactorsMat, modelSolutions] = ONSAS_solve(modelCurrSol, modelProperties, BCsData);

displacementsIII = matUs((4) * 6 + 3, :) * 100; % node with vertical load applied
loadfactorsIII = loadFactorsMat(:, 2);

% Plots

lw = 2;
ms = 6;
plotfontsize = 14;
step = 75;

f1 = figure('Name', 'Frame / Plasticity (load factors)', 'NumberTitle', 'off');
hold on;
grid on;

plot(abs(displacements_perfect(1:step:end)), loadfactors_perfect(1:step:end), '-s', 'linewidth', lw, 'markersize', ms, "Color", "#77AC30");
plot(abs(displacements(1:step:end)), loadfactors(1:step:end), '-s', 'linewidth', lw, 'markersize', ms, "Color", "#D95319");
plot(abs(displacementsII(1:step:end)), loadfactorsII(1:step:end), '-o', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120");
plot(abs(displacementsIII(1:step:end)), loadfactorsIII(1:step:end), '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD");

labx = xlabel('Desplazamiento u_v (cm)');
laby = ylabel('\lambdaF (kN)');

legend('ONSAS \lambdaF(u_v), a = 0', 'ONSAS \lambdaF(u_v), a = -0,04', 'ONSAS \lambdaF(u_v), a = -0,06', 'ONSAS \lambdaF(u_v), a = -0,0718', 'location', 'southeast');

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize);
set(labx, 'FontSize', plotfontsize);
set(laby, 'FontSize', plotfontsize);
title('Darvall-Mendis Frame / Rotulas de Ablandamiento');

print('-f1', '../tex/imagenes/darvallmendisframeplot.png', '-dpng');
fprintf(' a = -0.0718\n');
fprintf(' First hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4) * 6 + 3, TZERO(4)) * 100, loadFactorsMat(TZERO(4), 2));
% fprintf(' Second hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4)*6+3,TZERO(7))*100,loadFactorsMat(TZERO(7),2)) ;
% fprintf(' Third hinge\n Displacement u = %g\n Load Factor lambdaF = %g\n', matUs((4)*6+3,TZERO(2))*100,loadFactorsMat(TZERO(2),2)) ;
