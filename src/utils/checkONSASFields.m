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
function checkONSASFields(materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams)

  checkFields(materials, {'modelName', 'modelParams', 'density', 'nodalMass'});

  checkFields(elements, {'elemType', 'elemTypeParams', 'massMatType', ...
                         'elemCrossSecParams', 'aeroNumericalParams', 'aeroCoefFunctions', 'chordVector'});

  checkFields(boundaryConds, {'loadsCoordSys', 'loadsTimeFact', 'loadsBaseVals', ...
                              'userLoadsFilename', 'imposDispDofs', 'imposDispVals', 'springDofs', 'springVals'});

  if length(initialConds) > 0
    if isfield(analysisSettings, 'crossFlowVIVBool') && isfield(analysisSettings, 'inLineVIVBool')
      checkFields(initialConds, {'U', 'Udot', 'Udotdot', 'Q0', 'P0'});
    elseif isfield(analysisSettings, 'inLineVIVBool')
      checkFields(initialConds, {'U', 'Udot', 'Udotdot', 'P0'});
    elseif isfield(analysisSettings, 'crossFlowVIVBool')
      checkFields(initialConds, {'U', 'Udot', 'Udotdot', 'Q0'});
    else
      checkFields(initialConds, {'U', 'Udot', 'Udotdot'});
    end
  end

  analysisSettingsDefaultFields = ...
   {'geometricNonLinearAero', 'numGaussPointsAeroForce', 'computeAeroStiffnessMatrix', ...
    'fluidProps', 'addedMassBool', 'booleanSelfWeight', 'methodName', 'deltaT', 'finalTime', 'incremArcLen', 'iniDeltaLamb', ...
    'stopTolDeltau', 'stopTolForces', 'stopTolIts', 'stabilityAnalysisFlag', 'modalAnalysisBoolean', 'posVariableLoadBC', ...
    'ALdominantDOF', 'crossFlowVIVBool', 'inLineVIVBool', 'constantLiftDir'};

  if isfield(analysisSettings, 'methodName')
    if strcmp(analysisSettings.methodName, 'newmark')
      analysisSettingsDefaultFields{end +1} = 'alphaNM';
      analysisSettingsDefaultFields{end +1} = 'deltaNM';
    elseif strcmp(analysisSettings.methodName, 'alphaHHT') || isfield(analysisSettings, 'alphaHHT')
      analysisSettingsDefaultFields{end +1} = 'alphaHHT';
    end
  end

  if isfield(analysisSettings, 'alphaHHT')
    analysisSettingsDefaultFields{end +1} = 'alphaHHT';
  end

  checkFields(analysisSettings, analysisSettingsDefaultFields);

  checkFields(otherParams, {'problemName', 'plots_format', ...
                            'plots_deltaTs_separation', 'controlDofs', 'storeBoolean', 'nodalDispDamping', 'exportFirstMatrices'});

function checkFields(mystruct, expectedFields)

  myfields = fieldnames(mystruct);
  len = length(expectedFields); % cell with expected fields
  for i = 1:length(myfields)
    j = 1;
    while (j <= len) && (strcmp(myfields{i}, expectedFields{j}) ~= 1)
      j = j + 1;
    end
    if j == (len + 1)
      error('The field %s is not a correct field name', myfields{i});
    end
  end
