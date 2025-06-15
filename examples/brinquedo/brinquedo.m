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
E = 2100000000;
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

nnodes = 3 ;
% fprintf('====================================================================================\n')
% fprintf('First Case - Node 1 & 2 stretch parallel to line 1-2\n')
% fprintf('====================================================================================\n')

% elemDisps_case1 = zeros(nnodes*6,1) ;
% dof_ux_1 = 1*6-5 ;
% dof_ux_2 = 2*6-5 ;
% ux_1 = -1e-8 ;
% ux_2 =  1e-8 ;

% elemDisps_case1(dof_ux_1,1) = ux_1 ;
% elemDisps_case1(dof_ux_2,1) = ux_2 ;

% % Linear
% [fsL,KL,~] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case1 , 'elastic-linear', [ E nu], tz);
% fsL = fsL{1};
% ksL = KL{1};

% % Non-Linear
% materials(1).modelName  = 'elastic-rotEngStr';
% [fsNL,KNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case1 , 'elastic-rotEngStr', [ E nu], tz);
% fsNL = fsNL{1};
% ksNL = KNL{1};

% % Difference
% dif_f_case1 = fsL ./ fsNL
% dif_K_case1 = ksL ./ ksNL

% eig(ksL)
% stop

% fprintf('====================================================================================\n')
% fprintf('Second Case - Node 3 stretch y-dir\n')
% fprintf('====================================================================================\n')

% elemDisps_case2 = zeros(nnodes*6,1) ;
% dof_uy_3 = 3*6-3 ;

% uy_3 = 1e-4 ;

% elemDisps_case2(dof_uy_3,1) = uy_3 ;

% % Linear
% [fsL,KL,~] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case2 , 'elastic-linear', [ E nu], tz);
% fsL = fsL{1};
% ksL = KL{1};

% % Non-Linear
% materials(1).modelName  = 'elastic-rotEngStr';
% [fsNL,KNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case2 , 'elastic-rotEngStr', [ E nu], tz);
% fsNL = fsNL{1};
% ksNL = KNL{1};

% % Difference
% dif_f_case2 = fsL ./ fsNL
% dif_K_case2 = ksL ./ ksNL

%  fprintf('====================================================================================\n')
%  fprintf('Second Case 2.1 - Node 2 stretch x-dir // Node 1 stretch x-dir\n')
%  fprintf('====================================================================================\n')

%  elemDisps_case21 = zeros(nnodes*6,1) ;
%  dof_ux_1 = 2*6-5 ;

%  ux_1 = 1e-4 ;

%  elemDisps_case21(dof_ux_1,1) = ux_1 ;

%  % Linear
%  [fsL,KL,~] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case21 , 'elastic-linear', [ E nu], tz);
%  fsL_1 = fsL{1};
%  ksL_1 = KL{1};

%  % Non-Linear
%  materials(1).modelName  = 'elastic-rotEngStr';
%  [fsNL,KNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case21 , 'elastic-rotEngStr', [ E nu], tz);
%  fsNL_1 = fsNL{1};
%  ksNL_1 = KNL{1};

%  % Difference
%  dif_f_case2 = fsL_1 ./ fsNL_1
%  dif_K_case2 = ksL_1 ./ ksNL_1
 
 
%  elemDisps_case22 = zeros(nnodes*6,1) ;
%  dof_ux_2 = 2*6-5 ;

%  ux_2 = 1e-4 ;

%  elemDisps_case22(dof_ux_2,1) = ux_2 ;

%  % Linear
%  [fsL,KL,~] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case22 , 'elastic-linear', [ E nu], tz);
%  fsL_2 = fsL{1};
%  ksL_2 = KL{1};

%  % Non-Linear
%  materials(1).modelName  = 'elastic-rotEngStr';
%  [fsNL,KNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case22 , 'elastic-rotEngStr', [ E nu], tz);
%  fsNL_2 = fsNL{1};
%  ksNL_2 = KNL{1};

%  % Difference
%  dif_f_case2 = fsL_2 ./ fsNL_2
%  dif_K_case2 = ksL_2 ./ ksNL_2
 
%  [ fsL_1 fsNL_1 fsL_2 fsNL_2 ]
 
%  stop

fprintf('====================================================================================\n')
fprintf('Third Case - Node 3 vertical disp z-dir\n')
fprintf('====================================================================================\n')



mesh.nodesCoords = [ -Lx/2  0       0  ; 
                      Lx/2  0       0  ; 
                      0     Ly      0 ];

% r1_g = mesh.nodesCoords(1,:)';
% r2_g = mesh.nodesCoords(2,:)';
% r3_g = mesh.nodesCoords(3,:)';
% rc_g = (r1_g+r2_g+r3_g) / 3 ;

% mesh.nodesCoords = [ r1_g-rc_g  ; 
%                      r2_g-rc_g  ; 
%                      r3_g-rc_g  ];

elemDisps_case3 = zeros(nnodes*6,1) ;
dof_rx_1 = 1*6-4 ;
dof_rx_2 = 2*6-4 ;
dof_rx_3 = 3*6-4 ;
dof_uz_3 = 3*6-1 ;
dof_uz_2 = 2*6-1 ;
dof_uz_1 = 1*6-1 ;

uz_3 = 1e-6 ;
uz_2 = 1e-6 ;
uz_1 = 1e-6 ;

tan_theta_1 = uz_3 / Ly ;
theta_1 = atan(tan_theta_1) ;

rx_1 = theta_1 ;
rx_2 = theta_1 ;
rx_3 = theta_1 ;

elemDisps_case3(dof_rx_1,1) = rx_1 ;
% elemDisps_case3(dof_rx_2,1) = rx_2 ;
% elemDisps_case3(dof_rx_3,1) = rx_3 ;
% elemDisps_case3(dof_uz_3,1) = uz_3 ;
% elemDisps_case3(dof_uz_2,1) = uz_2*5 ;
elemDisps_case3(dof_uz_1,1) = uz_1*5 ;

% Linear
[fsL,KL,~] = internalForcesLinearShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case3 , 'elastic-linear', [ E nu], tz);
fsL = fsL{1};
ksL = KL{1};

% Non-Linear
materials(1).modelName  = 'elastic-rotEngStr';

[fsNL,KNL,~] = internalForcesShellTriangle(reshape( mesh.nodesCoords', 1,9 ), elemDisps_case3 , 'elastic-rotEngStr', [ E nu], tz, []);
fsNL = fsNL{1};
ksNL = KNL{1};

% % Difference
dif_f_case3 = fsL ./ fsNL 
dif_K_case3 = ksL ./ ksNL 

f_NL = ksNL * elemDisps_case3;

% a = [ 'ux' ; 'rx' ; 'uy' ; 'ry' ; 'uz' ; 'rz' ] ;
maxDif = max(max(dif_K_case3));
pos = find(maxDif==max(dif_K_case3));

poss_max = [ pos maxDif ]

[ fsL fsNL f_NL ]

% Rr_ana = [  1   0               0               ;
%             0   cos(theta_1)    -sin(theta_1)    ;
%             0   sin(theta_1)    cos(theta_1)    ] ;

% Ly_def = sqrt(Ly^2 + uz_3^2) ;
% uy_3_bar = Ly_def - Ly ;

% r1_g = mesh.nodesCoords(1,:)';
% r2_g = mesh.nodesCoords(2,:)';
% r3_g = mesh.nodesCoords(3,:)';
% rc_g = (r1_g+r2_g+r3_g) / 3 ;

% u1_g = zeros(3,1) ;
% u2_g = zeros(3,1) ;
% u3_g = zeros(3,1) ;
% u3_g(3) = uz_3 ;

% p1_g = r1_g + u1_g ;
% p2_g = r2_g + u2_g ;
% p3_g = r3_g + u3_g ;
% pc_g = (p1_g+p2_g+p3_g) / 3 ;

% To = eye(3) ;

% r1_o = To * (r1_g-rc_g) ;
% r2_o = To * (r2_g-rc_g) ;
% r3_o = To * (r3_g-rc_g) ;
% [ r1_o r2_o r3_o ];

% u1_def = Rr_ana'*(p1_g-pc_g) - r1_o ;
% stop








# Nada que ver

% A=3;
% l=2;
% E=5;
% I=7;

% Ka = [   E*A/l      -E*A/l ;
%         -E*A/l      E*A/l ] ;


% Kb = E*I/l^3 * [    12      6*l     -12     6*l     ;
%                     6*l     4*l^2   -6*l    2*l^2   ;
%                     -12     -6*l    12      -6*l    ;
%                     6*l     2*l^2   -6*l    4*l^2  ] ;

% ia = [ 1 4 ] ;
% ib = [ 2 3 5 6] ;

% Kl = zeros(6,6) ;

% Kl(ia,ia) = Ka ;
% Kl(ib,ib) = Kb ;


% % Haugen
%   u = u1_g
%   c = pog-rc_g
%   c_aux = ucg

%   I = eye(3) ;
%   Ro = Rr ;
%   xo = r1g ;
%   a = rc_g ; 

%   ud = u-c+(I-Ro)*(xo-a)
%   ud_def = Rr'*ud
%   stop
