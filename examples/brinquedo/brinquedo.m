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
E = 2;
nu = 0.0;
tz = .1;

Lz = .5;
Lx = 1;

materials                    = struct();
materials(1).modelName  = 'elastic-linear';
materials(1).modelParams = [E nu];

elements             = struct();
elements(1).elemType = 'node';
elements(2).elemType = 'triangle-shell';
elements(2).elemCrossSecParams = {'thickness', tz };

boundaryConds                  = struct();
boundaryConds(1).imposDispDofs =  [1 2 3 4 5 6];
boundaryConds(1).imposDispVals =  [0 0 0 0 0 0];

boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) t;
boundaryConds(2).loadsBaseVals = [0 0 1e-8 0 0 0];
# boundaryConds(2).loadsBaseVals = [0 0 0 0 1e-8 0];
# boundaryConds(2).loadsBaseVals = [0 1e-8 0 0 0 0];

mesh = struct();
mesh.nodesCoords = [ 0 0 0 ; Lx 0 0; Lx/2 0 Lz ];

mesh.conecCell = {}
mesh.conecCell{1,1} = [ 0 1 1  1  ];
mesh.conecCell{2,1} = [ 0 1 1  2  ];
mesh.conecCell{3,1} = [ 0 1 2  3  ];

mesh.conecCell{4,1} = [ 1 2 0  1 2 3  ];

% md### Initial conditions
initialConds                  = struct();
% md
% md#### Analysis settings
analysisSettings               = struct();
analysisSettings.methodName    = 'newtonRaphson';
analysisSettings.deltaT        =   1;
analysisSettings.finalTime     =   1;
analysisSettings.stopTolDeltau =   1e-10;
analysisSettings.stopTolForces =   1e-10;
analysisSettings.stopTolIts    =   10;
% md
% md#### OtherParams
% md The nodalDispDamping is added into the model using:
otherParams                  = struct();
% md The name of the problem is:
% md
otherParams.problemName  = 'brinquedoLin';
otherParams.plots_format = 'vtk';

# % md Execute ONSAS and save the results:
# [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
# %
# % mdAfter that the structs are used to perform the numerical time analysis
# [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);

# controldofs = 2*6+[ 2 3 ] ;
# dispsLIN_rotx_dispy = matUs( controldofs,end )

# disps_LIN_disz = matUs( 2*6+ 5,end )

# otherParams.problemName  = 'brinquedo_naolinear';




materials(1).modelName  = 'elastic-rotEngStr';


# [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
# %
# % mdAfter that the structs are used to perform the numerical time analysis
# [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);

# controldofs = 2*6+[ 2 3 ] ;
# disps_NONLIN_rotx_dispy = matUs( controldofs,end )

# disps_NONLIN_disz = matUs( 2*6+ 5,end )


a = '-----------------------------------------------------'
#elemDisps = [zeros(12,1); zeros(4,1); 1e-4; 0 ]
elemDisps_rotx = [zeros(12,1); 0;  0e-4; 0; 0; 0; 0 ]

[fsNL,KNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_rotx , 'elastic-linear', [ E nu], tz);

fsNL = fsNL{1};
ksNL = KNL{1}

[fsL,KL,~] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_rotx , 'elastic-linear', [ E nu], tz);

fsL = fsL{1};
ksL = KL{1}

dif_K_rotx = (ksL - ksNL)  ./ ksNL

fsLvsfsNL_rotx = [ fsL fsNL fsL./fsNL ]

dif = fsL - fsNL;

% b = '-----------------------------------------------------'

% norm(dif)


% elemDisps_desy = [zeros(12,1); 0; 0; 1e-4; 0; 0; 0 ]


% [fsNL,KlNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_desy , 'elastic-linear', [ E nu], tz);

% fsNL = fsNL{1};

% [fsL,KeL,~] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_desy , 'elastic-linear', [ E nu], tz);

% fsL = fsL{1};

% fsLvsfsNL_desy = [ fsL fsNL fsL./fsNL ]

% dif = fsL - fsNL;

% norm(dif)
