%md function used to conver .m files in the ONSAS repo to .md files in this repo
%md this function is executed from the makeAndPreview.sh file or can also be
%md executed from the bash terminal: `octave --eval "bringONSASmFilesToONSAS_docs('$ONSAS_PATH')"`
%md where \$ONSAS_PATH is the environment variable with the directory of ONSAS.m

function bringONSASmFilesToONSASdocs
  disp('running bringONSASmFilesToONSASdocs script...')
  dirONSASdocs = './examples/' ;
  dirONSASm     = '../../examples/'          ;

  addpath( genpath( '../../src/examples/' ));

  ONSASmFiles = {   'staticVonMisesTruss/onsasExample_staticVonMisesTruss.m' ...
                  ; 'uniformCurvatureCantilever/onsasExample_uniformCurvatureCantilever.m' ...
                  ; 'uniaxialExtension/uniaxialExtension.m' ...
                  ; 'springMass/springMass.m' ...
                  ; 'semiSphereWithInclusion/semiSphereWithInclusion.m' ...
                  ; 'linearAerodynamics/linearAerodynamics.m' ...                
                  ; 'nonLinearAerodynamics/nonLinearAerodynamics.m' ...                
                } ;

  MDFiles    = {   'staticVonMisesTruss.md' ...
                 ; 'cantileverBeam.md' ...
                 ; 'uniaxialExtension.md' ...
                 ; 'springMass.md' ...
                 ; 'semiSphereWithInclusion.md' ...
                 ; 'linearAerodynamics.md' ...
                 ; 'nonLinearAerodynamics.md' ...
                } ;

  if exist( dirONSASdocs ) ~= 7
    fprintf('creating examples dir...\n')
    mkdir( './examples/' );
  end

  fprintf( 'converting:\n' )
  for i=1:length( ONSASmFiles )
    fprintf([ '  - ' ONSASmFiles{i} '\n' ])
    m2md( [ dirONSASm ONSASmFiles{i} ] , [ dirONSASdocs MDFiles{i} ] , 1, 1 ) ;
  end