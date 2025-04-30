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
clc; close all; clear all;
addpath(genpath([pwd '/../../src']));
% md
% md## Scalars
E = 10000;
nu = 0.0;
tz = .1;

lx = 10 ;
ly = 1.0 ; 

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
P = -.0001
boundaryConds(2).loadsBaseVals = [0 0 0 0 P 0];

mesh = struct();
base_dir = '';
if strcmp(getenv('TESTS_RUN'), 'yes') && isfolder('examples')
  base_dir = ['.' filesep 'examples' filesep  'brinquedo' filesep];
end
[mesh.nodesCoords, mesh.conecCell] = meshFileReader([base_dir 'geometry_cantilever.msh']);
assert(max(mesh.nodesCoords(:, 1)) == lx && max(mesh.nodesCoords(:, 2)) == ly);


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
 [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
# %
 % mdAfter that the structs are used to perform the numerical time analysis
 [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);

nnodes = size(mesh.nodesCoords,1) ;
us_5L=matUs( ((nnodes-2)*6+1):(nnodes-1)*6, end );
us_6L=matUs( ((nnodes-1)*6+1):(nnodes)*6, end );

% % ========================================================
% otherParams.problemName  = 'brinquedo_naolinear';
% materials(1).modelName  = 'elastic-rotEngStr';


%  [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%  %
%  % mdAfter that the structs are used to perform the numerical time analysis
%  [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);

I = ly*tz^3/12 ;
fz = 3*P*(lx)^3/(3*E*I) 
thetax = 3*P*(lx)^2/(2*E*I) 

Mx = 3*P*lx*2
Fy = 3*P 

uy_dof = [ 5*6-3 ] ;
rotx_dof = 4*6+2 ;

disps_LIN_uy = matUs( uy_dof,end )
disps_LIN_rotx = matUs( rotx_dof,end )

% disps_NLIN_uy = matUs( uy_dof,end )
% disps_NLIN_rotx = matUs( rotx_dof,end )

% us_5NL=matUs( ((nnodes-2)*6+1):(nnodes-1)*6, end );
% us_6NL=matUs( ((nnodes-1)*6+1):(nnodes)*6, end );

% [us_5L us_5NL]
% [us_6L us_6NL]

#elemDisps = [zeros(12,1); zeros(4,1); 1e-4; 0 ]
% elemDisps_rotx = [zeros(12,1); 0;  1e-4; 0; 0; 0; 0 ]
% u = [0;  -thetax; 0; 0; -fz; 0] ;

% elemDisps_m = [zeros(12,1); u ; u ] ;

% [fsNL,~,~,KlNL] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,12 )(1:9), elemDisps_m(1:18) , 'elastic-linear', [ E nu], tz);

% fsNL = fsNL{1};

% [fsL,~,~,KeL] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,12 )(1:9), elemDisps_m(1:18) , 'elastic-linear', [ E nu], tz);

% fsL = fsL{1};

% fsLvsfsNL_rotx = [ fsL fsNL fsL./fsNL ]

% dif = fsL - fsNL;

% norm(dif)


% elemDisps_desy = [zeros(12,1); 0; 0; 1e-4; 0; 0; 0 ]


% fsNL = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_desy , 'elastic-linear', [ E nu], tz);

% fsNL = fsNL{1};

% fsL = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_desy , 'elastic-linear', [ E nu], tz);

% fsL = fsL{1};

% fsLvsfsNL_desy = [ fsL fsNL fsL./fsNL ]

% dif = fsL - fsNL;

% norm(dif)
