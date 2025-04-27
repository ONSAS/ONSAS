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

function [modelCurrSol, modelProperties, BCsData] = initONSAS(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams)

  % checks if the fields defined are correct or not
  checkONSASFields(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams);

  % md set defaults
  [materials, elements, boundaryConds, analysisSettings, otherParams] = setDefaults(materials, elements, boundaryConds, analysisSettings, otherParams);

  % md sets the current version and welcomes user
  ONSASversion = '0.3.3';
  welcomeMessage(ONSASversion, otherParams);

  % creates outputdir in current location
  createOutputDir(otherParams);

  % =================================================================
  % md process boundary conds information and construct BCsData struct
  [Conec, Nodes, factorLoadsFextCell, loadFactorsFuncCell, diriDofs, neumDofs, KS, userLoadsFilename] ...
    = boundaryCondsProcessing(mesh, materials, elements, boundaryConds, analysisSettings);

  BCsData = constructBCsData(factorLoadsFextCell, loadFactorsFuncCell, neumDofs, KS, userLoadsFilename);
  % =================================================================

  nTimes = round(analysisSettings.finalTime / analysisSettings.deltaT) + 1; % number of times (including t=0)
  % if length( otherParams.plotParamsVector ) > 1
  %   nplots = min( [ nTimes otherParams.plotParamsVector(2) ] ) ;
  % else
  %   % default value: all
  nplots = nTimes;
  % end
  timesPlotsVec = round(linspace(1, nTimes, nplots)');

  % =================================================================
  % md construct modelProperties struct
  modelProperties = constructModelProperties(Nodes, Conec, materials, elements, analysisSettings, otherParams, timesPlotsVec);
  % =================================================================

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TEMPORARY
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  global spitMatrices
  if spitMatrices == true
    save('-mat', 'output/loads.mat', 'factorLoadsFextCell');
  end

  % Initialize VIV-related vectors if VIV is enabled
  if any(analysisSettings.crossFlowVIVBool) || any(analysisSettings.inLineVIVBool)
    numElements = size(modelProperties.Conec, 1);
    dofsPerElement = 2;
    numTimeSteps = round(analysisSettings.finalTime / analysisSettings.deltaT) + 1;
    if analysisSettings.crossFlowVIVBool
      global qvect
      qvect = zeros(numElements * dofsPerElement, numTimeSteps);
      qvect(:, 1) = initialConds.Q0;
    end
    if analysisSettings.inLineVIVBool
      global pvect
      pvect = zeros(numElements * dofsPerElement, numTimeSteps);
      pvect(:, 1) = initialConds.P0;
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % =================================================================
  % md process initial conds and construct modelSol struct
  currTime  = 0;
  timeIndex = 1;

  timeStepIters    = 0;
  timeStepStopCrit = 0;

  % process initial conditions
  nNodes = size(mesh.nodesCoords, 1);
  [U, Udot, Udotdot] = initialCondsProcessing(initialConds, nNodes);

  convDeltau   = zeros(size(U));

  % ~ previousStateCell = zeros( size(Conec,1), 3 ) ; % assumed only for trusses: scalar per element
  previousStateCell = cell(size(Conec, 1), 3);
  previousStateCell(:, 1) = {zeros(1, 3)};
  previousStateCell(:, 2) = {zeros(1, 3)};
  previousStateCell(:, 3) = {0};

  % comput internal forces and stresses
  [~, Stress, ~, localInternalForces, strain_vec, acum_plas_strain_vec] = assembler(modelProperties.Conec, modelProperties.elements, modelProperties.Nodes, modelProperties.materials, BCsData.KS, U, Udot, Udotdot, modelProperties.analysisSettings, [0 1 0 1], modelProperties.nodalDispDamping, currTime, previousStateCell);

  [FextG, currLoadFactorsVals]  = computeFext(modelProperties, BCsData, 0, length(U), [], {U, Udot, Udotdot});

  nextTime = currTime + analysisSettings.deltaT;

  % md call assembler
  [systemDeltauMatrix, systemDeltauRHS, ~, ~, ~, ~, modelProperties.exportFirstMatrices] = systemAssembler(modelProperties, BCsData, U, Udot, Udotdot, U, Udot, Udotdot, nextTime, [], previousStateCell);

  modelCurrSol = constructModelSol(timeIndex, currTime, U, Udot, Udotdot, Stress, convDeltau, ...
                                   currLoadFactorsVals, systemDeltauMatrix, systemDeltauRHS, timeStepStopCrit, timeStepIters, localInternalForces, previousStateCell);
  % =================================================================

  % md prints headers for solver output file
  printSolverOutput(otherParams.outputDir, otherParams.problemName, 0, []);
  printSolverOutput(otherParams.outputDir, otherParams.problemName, [2 timeIndex currTime 0 0], []);

  % md writes vtk file
  if strcmp(modelProperties.plots_format, 'vtk')
    vtkMainWriter(modelCurrSol, modelProperties);
  end

  if exist('controlDofs') == 0
    controlDofs = [];
    controlDofsAndFactors = [];
  end

  if length(controlDofs) > 0
    controlDofsAndFactors = zeros(size(controlDofs, 1), 2);

    % control dof info
    for i = 1:size(controlDofs, 1)
      aux                = nodes2dofs(controlDofs(i, 1), 6);
      controlDofsAndFactors(i, :) = [aux(controlDofs(i, 2)) controlDofs(i, 3)];
    end
  end

  % =========================================
  % function for creation of output directory
  % -----------------------------------------
function createOutputDir(otherParams)

  outputDir = otherParams.outputDir;

  if exist('./output/') ~= 7
    if otherParams.screenOutputBool
      fprintf('  - Creating directory ./output/ ...');
    end

    mkdir('./', './output/');

    if otherParams.screenOutputBool
      fprintf(' done. \n');
    end
  end

  % -----------------
  if exist(outputDir) == 7 % problemName is a directory
    % the content is erased
    if otherParams.screenOutputBool
      fprintf(['|  - Cleaning output directory ...']);
    end
    if isThisOctave
      confirm_recursive_rmdir(0);
    end

    % delete
    [aux, msg] = rmdir(outputDir, 's');

    % create empty
    mkdir(outputDir);

  elseif exist(['./' otherParams.problemName '/']) ~= 7 % problemName is not a directory
    % it is created
    if otherParams.screenOutputBool
      fprintf(['|  - Creating output directory ...']);
    end
    mkdir(outputDir);
  end
  if otherParams.screenOutputBool
    fprintf(' done. \n');
  end

  % =========================================
  % function for welcome message
  % -----------------------------------------
function welcomeMessage(ONSASversion, otherParams)
  if otherParams.screenOutputBool
    fprintf(['\n' ...
             '|=================================================|\n' ...
             '|         _ _             _ _     _ _     _ _     |\n' ...
             '|       /    /  /|   /  /       /    /  /         |\n' ...
             '|      /    /  / |  /  /_ _    /_ _ /  /_ _       |\n' ...
             '|     /    /  /  | /       /  /    /       /      |\n' ...
             '|    /_ _ /  /   |/   _ _ /  /    /   _ _ /       |\n' ...
             '|                                                 |\n' ...
             '|-------------------------------------------------|\n']);
    fprintf(['| Welcome to ONSAS v' ONSASversion '.                        |\n' ...
             '| This program comes with ABSOLUTELY NO WARRANTY. |\n' ...
             '| Please read the COPYING.txt and README.md files |\n' ...
             '|-------------------------------------------------|\n']);
    fprintf(['| Solving problem:  ' otherParams.problemName '\n']);
  end
