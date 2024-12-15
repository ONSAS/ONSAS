% md function used to conver .m files in the ONSAS repo to .md files in this repo
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

end
