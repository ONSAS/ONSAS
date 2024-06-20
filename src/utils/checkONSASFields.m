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
function checkONSASFields( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams )

  checkFields(materials, {'modelName', 'modelParams', 'density', 'nodalMass'})

  checkFields(elements, {'elemType', 'elemTypeParams','massMatType',...
                            'elemCrossSecParams','aeroNumericalParams','aeroCoefFunctions','chordVector'});

   analysisSettingsDefaultFields = {'geometricNonLinearAero', 'numGaussPointsAeroForce','computeAeroStiffnessMatrix',...
    'fluidProps','addedMassBool','booleanSelfWeight','methodName','deltaT', 'finalTime', 'incremArcLen', 'iniDeltaLamb',...
     'stopTolDeltau', 'stopTolForces', 'stopTolIts','stabilityAnalysisFlag', 'modalAnalysisBoolean','posVariableLoadBC',...
      'ALdominantDOF'};
    if isfield(analysisSettings, 'methodName')
        if strcmp( analysisSettings.methodName, 'newmark' )
            analysisSettingsDefaultFields{end +1} = 'alphaNM' ; 
            analysisSettingsDefaultFields{end +1} = 'deltaNM' ; 
        end
        if strcmp( analysisSettings.methodName, 'alphaHHT' ) || isfield(analysisSettings, 'alphaHHT')
            analysisSettingsDefaultFields{end +1} = 'alphaHHT' ;
        end                      
    end
    analysisSettings
    if isfield(analysisSettings, 'alphaHHT')
        analysisSettingsDefaultFields{end +1} = 'alphaHHT'
    end
    analysisSettingsDefaultFields

  checkFields(analysisSettings, analysisSettingsDefaultFields)


function checkFields(mystruct, expectedFields)

    myfields = fieldnames(mystruct);
    len = length(expectedFields) ; % cell with expected fields
    for i = 1:length(myfields)
        j = 1;
        while (j <= len) && (strcmp(myfields{i},expectedFields{j}) ~= 1)
            j = j+1;
        end
        if j == (len+1)
            error('The field %s is not a correct field name', myfields{i}) 
        end     
    end
