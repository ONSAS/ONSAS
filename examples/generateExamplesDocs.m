
% llamada jorge:    generateExamplesDocs('~/.julia/dev/ONSAS_docs/docs/src/')

function generateExamplesDocs( dirONSAS_docs )

addpath(genpath('../src/'));

fprintf( 'converting:\n' )

ONSASFiles = {'./staticVonMisesTruss/onsasExample_staticVonMisesTruss.m', ...
              './uniformCurvatureCantilever/onsasExample_uniformCurvatureCantilever.m' , ...
              './uniaxialExtension/onsasExample_uniaxialExtension.m' ...
              } ;

MDFiles    = {  'staticVonMisesTruss.md' ...
              , 'cantileverBeam.md' ...
              , 'uniaxialExtension.md' ...
              } ;

for i=1:length(ONSASFiles)
  fprintf([ '  - ' ONSASFiles{i} '\n' ])
  m2md( ONSASFiles{i} , [ dirONSAS_docs MDFiles{i} ] , 1 )
end
