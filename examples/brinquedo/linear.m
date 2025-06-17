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

close all; clear all; clc
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

% Mesh struct
mesh = struct();

mesh.nodesCoords = [ -Lx/2      0       0  ; 
                      Lx/2      0       0  ; 
                      0         Ly      0  ];

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
otherParams.problemName  = 'cantilever_linear';
otherParams.plots_format = 'vtk';

% general
P = 1 ;

dofs_3 = ((6*3-5):6*3) ;
% Pz+ case
boundaryConds(2).loadsBaseVals = [0 0 0 0 P 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_zp = matUs( (6*3-5):(6*3) , end) ;

% Pz- case
boundaryConds(2).loadsBaseVals = [0 0 0 0 -P 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_zm = matUs( (6*3-5):(6*3) , end) ;

[u_3_zp u_3_zm]

% Px+ case

boundaryConds(2).loadsBaseVals = [P 0 0 0 0 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_xp = matUs( (6*3-5):(6*3) , end) ;

% Px- case
boundaryConds(2).loadsBaseVals = [-P 0 0 0 0 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_xm = matUs( (6*3-5):(6*3) , end) ;

[u_3_xp u_3_xm]

% Py+ case
boundaryConds(2).loadsBaseVals = [0 0 P 0 0 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_yp = matUs( (6*3-5):(6*3) , end) ;

% Py- case
boundaryConds(2).loadsBaseVals = [0 0 -P 0 0 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_ym = matUs( (6*3-5):(6*3) , end) ;

[u_3_yp u_3_ym]

% Mx+ case
boundaryConds(2).loadsBaseVals = [0 P 0 0 0 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_mxp = matUs( (6*3-5):(6*3) , end) ;

% Mx- case
boundaryConds(2).loadsBaseVals = [0 -P 0 0 0 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_mxm = matUs( (6*3-5):(6*3) , end) ;

[u_3_mxp u_3_mxm]

% My case
boundaryConds(2).loadsBaseVals = [0 0 0 P 0 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_myp = matUs( (6*3-5):(6*3) , end) ;

% My- case
boundaryConds(2).loadsBaseVals = [0 0 0 -P 0 0];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_mym = matUs( (6*3-5):(6*3) , end) ;

[u_3_myp u_3_mym]

% Mz case
boundaryConds(2).loadsBaseVals = [0 0 0 0 0 P];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_mzp = matUs( (6*3-5):(6*3) , end) ;

% Mz- case
boundaryConds(2).loadsBaseVals = [0 0 0 0 0 -P];
[modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
u_3_mzm = matUs( (6*3-5):(6*3) , end) ;

[u_3_mzp u_3_mzm]

% =====================================
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

% close all; clear all; clc
% addpath(genpath([pwd '/../../src']));
% % md
% % md## Scalars
% E = 2100000000;
% nu = 0.0;
% tz = .1;

% Ly = 2;
% Lx = 2;


% materials                   = struct();
% materials(1).modelName      = 'elastic-linear';
% materials(1).modelParams    = [E nu];

% elements             = struct();
% elements(1).elemType = 'node';
% elements(2).elemType = 'triangle-shell';
% elements(2).elemCrossSecParams = {'thickness', tz };

% boundaryConds                  = struct();
% boundaryConds(1).imposDispDofs =  [1 2 3 4 5 6];
% boundaryConds(1).imposDispVals =  [0 0 0 0 0 0];

% boundaryConds(2).loadsCoordSys = 'global';
% boundaryConds(2).loadsTimeFact = @(t) t;

% % Mesh struct
% mesh = struct();

% mesh.nodesCoords = [ -Lx/2      0       0  ; 
%                       Lx/2      0       0  ; 
%                      -Lx/2     Ly       0  ;
%                       Lx/2     Ly       0  ];

% mesh.conecCell = {} ;
% mesh.conecCell{1,1} = [ 0 1 1  1  ];
% mesh.conecCell{2,1} = [ 0 1 1  2  ];
% mesh.conecCell{3,1} = [ 0 1 2  3  ];
% mesh.conecCell{4,1} = [ 0 1 2  4  ];
% mesh.conecCell{5,1} = [ 1 2 0  1 2 3  ];
% mesh.conecCell{6,1} = [ 1 2 0  2 4 3  ];

% % Initial conditions
% initialConds                  = struct();

% % Analysis settings

% analysisSettings               = struct();
% analysisSettings.methodName    = 'newtonRaphson';
% analysisSettings.deltaT        =   1;
% analysisSettings.finalTime     =   1;
% analysisSettings.stopTolDeltau =   1e-10;
% analysisSettings.stopTolForces =   1e-10;
% analysisSettings.stopTolIts    =   10;

% otherParams                  = struct();
% otherParams.problemName  = 'cantilever_linear';
% otherParams.plots_format = 'vtk';

% % general
% P = 1 ;

% dofs_3 = ((6*3-5):6*3) ;
% dofs_4 = ((6*4-5):6*4) ;
% % Pz+ case
% boundaryConds(2).loadsBaseVals = [0 0 0 0 P 0];
% [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% u_3_zp = matUs( (6*3-5):(6*3) , end) ;
% u_4_zp = matUs( (6*4-5):(6*4) , end) ;
% % modelSolutions = modelSolutions{1};
% % loadFactorsMat
% % modelSolutions.localInternalForces
% % modelSolutions.currLoadFactorsVals


% % Pz- case
% boundaryConds(2).loadsBaseVals = [0 0 0 0 -P 0];
% [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% u_3_zm = matUs( (6*3-5):(6*3) , end) ;
% u_4_zm = matUs( (6*4-5):(6*4) , end) ;
% % modelSolutions = modelSolutions{1};
% % loadFactorsMat
% % modelSolutions.localInternalForces
% % modelSolutions.currLoadFactorsVals
% [u_3_zp u_3_zm u_4_zp u_4_zm]

% % ================
% % Pz+- case
% mesh.conecCell{3,1} = [ 0 1 2  3  ];
% mesh.conecCell{4,1} = [ 0 1 3  4  ];
% boundaryConds(2).loadsBaseVals = [0 0 0 0 P 0];
% boundaryConds(3).loadsBaseVals = [0 0 0 0 -P 0];
% [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% u_3_zpm = matUs( (6*3-5):(6*3) , end) ;
% u_4_zpm = matUs( (6*4-5):(6*4) , end) ;

% % Pz-+ case
% boundaryConds(2).loadsBaseVals = [0 0 0 0 -P 0];
% boundaryConds(3).loadsBaseVals = [0 0 0 0 P 0];
% [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% u_3_zmp = matUs( (6*3-5):(6*3) , end) ;
% u_4_zmp = matUs( (6*4-5):(6*4) , end) ;

% [u_3_zpm u_4_zmp u_3_zmp u_4_zpm]


% % Px+ case

% boundaryConds(2).loadsBaseVals = [P 0 0 0 0 0];
% [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% u_3_xp = matUs( (6*3-5):(6*3) , end) ;
% u_4_xp = matUs( (6*4-5):(6*4) , end) ;

% % Px- case
% boundaryConds(2).loadsBaseVals = [-P 0 0 0 0 0];
% [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% u_3_xm = matUs( (6*3-5):(6*3) , end) ;
% u_4_xm = matUs( (6*4-5):(6*4) , end) ;

% [u_3_xp u_3_xm u_4_xp u_4_xm]

% % Py+ case
% boundaryConds(2).loadsBaseVals = [0 0 P 0 0 0];
% [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% u_3_yp = matUs( (6*3-5):(6*3) , end) ;
% u_4_yp = matUs( (6*4-5):(6*4) , end) ;

% % Py- case
% boundaryConds(2).loadsBaseVals = [0 0 -P 0 0 0];
% [modelInitSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
% [matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelInitSol, modelProperties, BCsData);
% u_3_ym = matUs( (6*3-5):(6*3) , end) ;
% u_4_ym = matUs( (6*4-5):(6*4) , end) ;

% [u_3_yp u_3_ym u_4_yp u_4_ym]