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
% md set optional fields defaults
function [materials, elements, boundaryConds, analysisSettings, otherParams] = setDefaults(materials, elements, boundaryConds, analysisSettings, otherParams)

  % =========================================
  % materials
  materials         = checkOrSetDefault (materials, 'density', 0);

  % =========================================
  % elements
  elements = checkOrSetDefault(elements, 'massMatType', 'lumped');
  elements = checkOrSetDefault(elements, 'elemTypeParams', []);
  elements = checkOrSetDefault(elements, 'elemCrossSecParams', []);
  elements = checkOrSetDefault(elements, 'chordVector', []);
  elements = checkOrSetDefault(elements, 'aeroCoefFunctions', {[], [], []});

  elements = setDefauAeroCrossSecProps(elements);

  % =========================================
  % boundaryConds
  boundaryConds    =  checkOrSetDefault (boundaryConds, 'loadsTimeFact', []);
  boundaryConds    =  checkOrSetDefault (boundaryConds, 'loadsCoordSys', []);
  boundaryConds    =  checkOrSetDefault (boundaryConds, 'springDofs', []);

  % =========================================
  % analysisSettings
  analysisSettings  = checkOrSetDefault (analysisSettings, 'geometricNonLinearAero', true);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'numGaussPointsAeroForce', 4);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'computeAeroStiffnessMatrix', false);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'fluidProps', []);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'addedMassBool', false);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'booleanSelfWeight', false);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'Utp10', []);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'methodName', 'newtonRaphson');
  analysisSettings  = checkOrSetDefault (analysisSettings, 'deltaT', 1);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'finalTime', 1);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'stopTolDeltau', 1e-6);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'stopTolForces', 1e-6);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'stopTolIts', 15);
  if strcmp(analysisSettings.methodName, 'newmark')
    analysisSettings = checkOrSetDefault(analysisSettings, 'alphaNM', 0.25);
    analysisSettings = checkOrSetDefault(analysisSettings, 'deltaNM', 0.50);
  end
  if strcmp(analysisSettings.methodName, 'alphaHHT')
    analysisSettings = checkOrSetDefault(analysisSettings, 'alphaHHT', -0.05);
  end
  analysisSettings  = checkOrSetDefault (analysisSettings, 'stabilityAnalysisFlag', 0);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'modalAnalysisBoolean', 0);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'ALdominantDOF', []);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'incremArcLen', 1.5);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'iniDeltaLamb', 1e-3);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'posVariableLoadBC', []);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'crossFlowVIVBool', []);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'inLineVIVBool', []);
  analysisSettings  = checkOrSetDefault (analysisSettings, 'constantLiftDir', []);
  % =========================================
  % otherParams
  otherParams       = checkOrSetDefault(otherParams, 'screenOutputBool', 1);
  otherParams       = checkOrSetDefault(otherParams, 'plots_format', []);
  otherParams       = checkOrSetDefault(otherParams, 'plots_deltaTs_separation', 1);
  otherParams       = checkOrSetDefault(otherParams, 'nodalDispDamping', 0);
  otherParams       = checkOrSetDefault(otherParams, 'outputDir', ['./output/' otherParams.problemName '/']);
  otherParams       = checkOrSetDefault(otherParams, 'exportFirstMatrices', false);

  % md function that checks if a field is defined in a (scalar or array) struct
  % md and sets a default value if it is not defined.
function structName = checkOrSetDefault(structName, fieldName, default)

  if ~isfield(structName, fieldName)
    for i = 1:length(structName)
      aux(i)  = setfield(structName(i), fieldName, default);
    end
    structName = aux;
  end
