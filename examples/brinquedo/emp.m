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

close all; clear all;
addpath(genpath([pwd '/../../src']));
% md
% md## Scalars
E = 100e9;
nu = 0.0;
tz = .05;

Ly = 1;
Lx = 0.5;


materials                   = struct();
materials(1).modelParams    = [E nu];

elements             = struct();

elements(2).elemType = 'triangle-shell';
elements(2).elemCrossSecParams = {'thickness', tz };

boundaryConds                  = struct();
boundaryConds(1).imposDispDofs =  [1 2 3 4 5 6];
boundaryConds(1).imposDispVals =  [0 0 0 0 0 0];

boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) t;

Pz = -10 ;
Py = 0 ;

boundaryConds(2).loadsBaseVals = [0 0 Py 0 Pz 0];

% Mesh struct
mesh = struct();

mesh.nodesCoords = [ -Lx/2  0       0  ; 
                      Lx/2  0       0  ; 
                      0     Ly      0  ];

elements(1).elemType = 'node';

mesh.conecCell = {} ;
mesh.conecCell{1,1} = [ 0 1 1  3  ];
mesh.conecCell{2,1} = [ 0 1 2  1  ];
mesh.conecCell{3,1} = [ 0 1 2  2  ];
mesh.conecCell{4,1} = [ 1 2 0  1 2 3  ];

% Initial conditions
initialConds                  = struct();

% Analysis settings

analysisSettings               = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.deltaT        =   1;
analysisSettings.finalTime     =   1;
analysisSettings.stopTolDeltau =   1e-10;
analysisSettings.stopTolForces =   1e-10;
analysisSettings.stopTolIts    =   2;

otherParams                  = struct();
otherParams.problemName  = 'extensionNonLin';
otherParams.plots_format = 'vtk';

nnodes = 3 ;
fprintf('====================================================================================\n')
fprintf('First Case \n')
fprintf('====================================================================================\n')


materials(1).modelName      = 'elastic-linear';
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
dz_max_L = min(matUs(5:6:end));
tx_L     = min(matUs(2:6:end));
dy_max_L = max(matUs(3:6:end));
% stop

Us = matUs(:,end) ;
nodes_coords = mesh.nodesCoords([1;2;3],:) ;
[fsL,KL,~] = internalForcesLinearShellTriangle(reshape( nodes_coords', 1,9 ), Us , 'elastic-linear', [ E nu], tz);

fL = fsL{1};

rotMat = cell(3,1) ;
rotMat(:) = eye(3) ;
[fsNL,KNL,~] = internalForcesShellTriangle(reshape( nodes_coords', 1,9 ), Us , 'elastic-rotEngStr', [ E nu], tz, rotMat);
fsNL = fsNL{1} ;

[ fL fsNL ]




materials(1).modelName  = 'elastic-rotEngStr';
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
dz_max_NL = min(matUs(5:6:end));
tx_NL     = min(matUs(2:6:end));
dy_max_NL = max(matUs(3:6:end));


% [fs, ks, fintLocCoord, DD] = internalForcesPlateTriangle(reshape( nodes_coords', 1,9 ), Us, 'elastic-linear', [ E nu], tz);
% [fsL,KL,~, ~,Kb] = internalForcesLinearShellTriangle(reshape( nodes_coords', 1,9 ), Us , 'elastic-linear', [ E nu], tz);

% ks_J_L = DD;
% ks_f_L = Kb;

% ks_J_L ./ ks_f_L

% dz_max_NL = min(matUs(5:6:end));
% tx_NL     = min(matUs(2:6:end));
% dy_max_NL = max(matUs(3:6:end));

% [ dz_max_L dz_max_NL  dy_max_L dy_max_NL tx_L tx_NL ]

