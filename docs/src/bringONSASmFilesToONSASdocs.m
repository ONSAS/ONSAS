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

% md function used to convert .m files in the ONSAS repo to .md files in this repo
% md this function is executed from the makeAndPreview.sh file or can also be
% md executed from the bash terminal: `octave --eval "bringONSASmFilesToONSAS_docs('$ONSAS_PATH')"`
% md where \$ONSAS_PATH is the environment variable with the directory of ONSAS

function bringONSASmFilesToONSASdocs
  disp('running bringONSASmFilesToONSASdocs script...');
  dirONSASdocs = './examples/';
  dirONSASm     = '../../examples/';

  addpath(genpath('../../src/examples/'));

  ONSASmFiles = {   'beamLinearVibration/beamLinearVibration.m'
                 % ; 'linearAerodynamics/linearAerodynamics.m' ...
                 'dragBeamReconfiguration/dragBeamReconfiguration.m' ...
                 ; 'ringPlaneStrain/ringPlaneStrain.m'
                 % ; 'nonLinearAerodynamics/nonLinearAerodynamics.m' ...
                 'uniaxialExtension/uniaxialExtension.m' ...
                 ; 'uniaxialCompression/uniaxialCompression.m' ...
                 ; 'uniformCurvatureCantilever/uniformCurvatureCantilever.m' ...
                 ; 'simplePropeller/simplePropeller.m' ...
                 ; 'springMass/springMass.m' ...
                 ; 'staticVonMisesTruss/staticVonMisesTruss.m' ...
                };

  MDFiles    = {   'beamLinearVibration.md'
                %  ; 'linearAerodynamics.md' ...
                'dragBeamReconfiguration.md' ...
                ; 'ringPlaneStrain.md'
                % ; 'nonLinearAerodynamics.md' ...
                'uniaxialExtension.md' ...
                ; 'uniaxialCompression.md' ...
                ; 'cantileverBeam.md' ...
                ; 'simplePropeller.md' ...
                ; 'springMass.md' ...
                ; 'staticVonMisesTruss.md' ...
               };

  if exist(dirONSASdocs) ~= 7
    fprintf('creating examples dir...\n');
    mkdir('./examples/');
  end

  fprintf('converting:\n');
  for i = 1:length(ONSASmFiles)
    fprintf(['  - ' ONSASmFiles{i} '\n']);
    m2md([dirONSASm ONSASmFiles{i}], [dirONSASdocs MDFiles{i}], 1, 1);
  end

  if exist('./assets/generated/') ~= 7
    fprintf('creating generated figures dir...\n');
    mkdir('./assets/generated/');
  end

  movefile('../../examples/staticVonMisesTruss/output/vonMisesTrussCheck.png', './assets/generated/');
  movefile('../../examples/simplePropeller/output/verifPropeller.png', './assets/generated/');
  movefile('../../examples/springMass/output/springMassCheckU.png', './assets/generated/');
  movefile('../../examples/uniaxialCompression/output/verifCompression.png', './assets/generated/');
  movefile('../../examples/uniaxialExtension/output/verifUniaxial.png', './assets/generated/');
  movefile('../../examples/dragBeamReconfiguration/output/RvsCyCd.png', './assets/generated/');
  movefile('../../examples/dragBeamReconfiguration/output/defPlots.png', './assets/generated/');
  movefile('../../examples/dragBeamReconfiguration/output/zDisplacementVIV.png', './assets/generated/');
end
