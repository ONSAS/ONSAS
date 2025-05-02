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

Ly = 2;
Lx = 2;


materials                   = struct();
materials(1).modelName      = 'elastic-linear';
materials(1).modelParams    = [E nu];

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

% Mesh struct
mesh = struct();

mesh.nodesCoords = [ -Lx/2  0       0  ; 
                      Lx/2  0       0  ; 
                      0     Ly      0 ];

mesh.conecCell = {} ;
mesh.conecCell{1,1} = [ 0 1 1  1  ];
mesh.conecCell{2,1} = [ 0 1 1  2  ];
mesh.conecCell{3,1} = [ 0 1 2  3  ];
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
analysisSettings.stopTolIts    =   10;


otherParams                  = struct();
otherParams.problemName  = 'brinquedoLin';
otherParams.plots_format = 'vtk';





a = '-----------------------------------------------------'

%  
nnodes = 3 ;
elemDisps_case1 = zeros(nnodes*6,1) ;

dof_ux_1 = 1*6-5 ;
dof_ux_2 = 2*6-5 ;

ux_1 = -1e-6 ;
ux_2 =  1e-6 ;

elemDisps_case1(dof_ux_1,1) = ux_1 ;
elemDisps_case1(dof_ux_2,1) = ux_2 ;

elemDisps_case1 ;

% First Case - Node 1 & 2 stretch parallel to line 1-2
% Linear
[fsL,KL,~] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case1 , 'elastic-linear', [ E nu], tz);

fsL = fsL{1};
ksL = KL{1};

% Non-Linear
materials(1).modelName  = 'elastic-rotEngStr';
%[fsNL,KNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_rotx , 'elastic-rotEngStr', [ E nu], tz);
[fsNL,KNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case1 , 'elastic-rotEngStr', [ E nu], tz);

fsNL = fsNL{1};
ksNL = KNL{1};

% Difference
dif_f_case1 = fsL ./ fsNL
dif_K_case1 = (ksL - ksNL)  ./ ksNL

sym_kL = issymmetric(ksL, 1e-15)
sym_kNL = issymmetric(ksNL,1e-15)

aux = ksNL-ksNL'
% fsLvsfsNL_rotx = [ fsL fsNL fsL./fsNL ]

%dif = fsL - fsNL;
