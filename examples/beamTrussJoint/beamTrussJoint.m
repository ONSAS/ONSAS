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
% md# Beam truss joint example
close all;
if ~strcmp(getenv('TESTS_RUN'), 'yes')
  clear all;
end
addpath(genpath([pwd '/../../src']));
% mdThe goal of this example is to provide a minimal validation of the integration between truss and frame elements in the same model. The structure considered is formed by two elements (one truss (t) and one beam (b)) and small displacements are considered.
% md
% md Truss geometrical and material properties are:
Et = 1e9;
nu = 3;
dt = .05;
At = pi * dt^2 / 4;
lt = 1;
nut = 0.3;
% md and frame geometrical and material properties are:
Eb = Et / 3;
db = 5 * dt;
Ab = pi * db^2 / 4;
lb = .5;
nub = 0.3;
Ib = pi * db^4 / 64;
% md##Numerical solution
% md### MEB parameters
% md### materials
% mdSince the example contains two different type of materials the fields of the `materials` struct will have two entries. Although the structure develops small displacements a Rotated Engineering strain material constitutive behavior is considered.
materials            = struct();
materials(1).modelName  = 'elastic-linear';
materials(1).modelParams = [Et nu];
%
materials(2).modelName  = 'elastic-linear';
materials(2).modelParams = [Eb nu];

% md### elements
% md
elements            = struct();
% mdThree different types of elements are considered: node, frame and truss, defined as follows:
elements(1).elemType = 'node';
elements(2).elemType = 'truss';
elements(3).elemType = 'frame';
% mdgeometrical properties are only assigned to the truss and beam elements:
elements(2).elemCrossSecParams{1, 1} = 'circle';
elements(2).elemCrossSecParams{2, 1} = dt;
elements(3).elemCrossSecParams{1, 1} = 'circle';
elements(3).elemCrossSecParams{2, 1} = db;
% mdTruss number of elements
%
numElemT  = 1;
numNodesT = numElemT + 1;
% md and beam:
numElemB  = 10;
numNodesB = numElemB + 1;

% md
% md### boundaryConds
% md
% mdThe fixed frame BC:
boundaryConds                  = struct();
boundaryConds(1).imposDispDofs = [1 2 3 4 5 6];
boundaryConds(1).imposDispVals = [0 0 0 0 0 0];
% mdloaded BC:
boundaryConds(2).loadsCoordSys = 'global';
boundaryConds(2).loadsTimeFact = @(t) 1e4 * t;
boundaryConds(2).loadsBaseVals = [0 0 0 0 1 0];
% mdsupport BC for node at the base of the truss:
boundaryConds(3).imposDispDofs = [1 3 5];
boundaryConds(3).imposDispVals = [0 0 0];
% md
% md### mesh parameters
% md
% mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh = struct();
mesh.nodesCoords = [(0:(numElemB))' * lb / numElemB zeros(numElemB + 1, 1) zeros(numElemB + 1, 1)
                    lb * ones(numElemT, 1)       zeros(numElemT, 1)   -(1:(numElemT))' * lt / numElemT];
% md where the columns 1,2 and 3 correspond to $x$, $y$ and $z$ coordinates, respectively, and the row $i$-th corresponds to the coordinates of node $i$.

% md The conectivity struct using MEB nomenclature is defined by using the following auxiliar Element and Nodes matrix:
% md Then the entry of node $1$ is introduced:
mesh.conecCell{1, 1} = [0 1 1  1];
% md the first MEB parameter (Material) is set as _zero_ (since nodes dont have material). The second parameter corresponds to the Element, and a _1_ is set since `node` is the first entry of the  `elements.elemType` cell. For the BC index, we consider that node $1$ is fixed, then the first index of the `boundaryConds` struct is used.
% md A similar approach is used for node where the beam ends $numNodesB$,
mesh.conecCell{2, 1} = [0 1 2  numNodesB];
% md analogosly for node $numNodesB + numElemT$ only the boundary condition is changed:
mesh.conecCell{3, 1} = [0 1 3  numNodesB + numElemT];
% md
% mdTo define the conecCell of elements a auxiliar auxConecElem matrix is defined using MEB nomenclature:
% mdMEB frame elements
frameElements = [ones(numElemB, 1) * 2, ones(numElemB, 1) * 3, zeros(numElemB, 1), ...
                 (1:numElemB)', (2:numElemB + 1)'];
% mdMEB truss elements
trussElements = [ones(numElemT, 1) * 1, ones(numElemT, 1) * 2, zeros(numElemT, 1), ...
                 (numElemB + 1:numElemB + numElemT)', (numElemB + 2:numElemB + numElemT + 1)'];
% mdCombine both element sets
auxConecElem = [frameElements; trussElements];
for i =  1:numElemB + numElemT
  mesh.conecCell{3 + i, 1} = auxConecElem(i, :);
end
% md### initial Conditions
% md homogeneous initial conditions are considered, then an empty struct is set:
initialConds                = struct();
% md### analysisSettings
% md The method used in the analysis is the Newton-Raphson, then the field `methodName` must be introduced as:
analysisSettings            = struct();
analysisSettings.methodName    = 'newtonRaphson';
% md and the following parameters correspond to the iterative numerical analysis settings
analysisSettings.deltaT        =   0.5;
analysisSettings.finalTime     =   1;
analysisSettings.stopTolDeltau =   1e-6;
analysisSettings.stopTolForces =   1e-6;
analysisSettings.stopTolIts    =   10;
% md
% md### otherParams
otherParams = struct();
otherParams.problemName = 'beamTrussJoint';
otherParams.plots_format = 'vtk';
% md
% md In order to validate this example the ONSAS code is run and the solution degree of freedom selected is the $uz$ displacement at the joint.
% mdFirst the input structs are converted to structs with the model information
[modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);
%
% mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions] = solveONSAS(modelCurrSol, modelProperties, BCsData);
% md

modelSolutions{3}.localInternalForces(11);

% md Its important to oultine analytical solution considering large dispalcements is not aviable. Thus the load applied is sleceted to produce amplitude displacements ($<0.05lt$). Consequently the small displacment solution of $uz$ is given by:
analyticFunc            = @(w)(Et * At / lt + 3 * Eb * Ib / lb^3) * w;
beamTruss_stiffRatio    = Et * At / lt / (3 * Eb * Ib / lb^3);
% mdThe numerical result is:
controlDof  = (numNodesB) * 6 - 1;
dispZnum    = matUs(controlDof, :);
% md and the verification boolean is computed as follows:
difLoadEngRot   = loadFactorsMat(:, 2)' - analyticFunc(dispZnum);
verifBoolean    = ((norm(difLoadEngRot) / norm(loadFactorsMat(:, 2))) <  1e-4);
% md
% md### Plots
% mdOutput diplacments and load factor function are plotted:
figure;
hold on;
grid on;
lw = 2.0;
lw2 = 1.0;
ms = 11;
plotfontsize = 18;
plot(dispZnum, analyticFunc(dispZnum), 'b-x', 'linewidth', lw, 'markersize', ms);
plot(dispZnum, loadFactorsMat(:, 2),    'r-s', 'linewidth', lw, 'markersize', ms);
labx = xlabel('Displacement u_z(t)');
laby = ylabel('\lambda(t)');
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize);
set(labx, 'FontSize', plotfontsize);
set(laby, 'FontSize', plotfontsize);
