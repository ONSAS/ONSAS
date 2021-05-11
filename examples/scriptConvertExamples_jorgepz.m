
addpath(genpath('../src/'));

fprintf( 'converting:\n' )

ONSASFiles = {'./staticVonMisesTruss/onsasExample_staticVonMisesTruss.m', ...
              './uniformCurvatureCantilever/onsasExample_uniformCurvatureCantilever.m' } ;
MDFiles    = {'~/.julia/dev/ONSAS_docs/docs/src/tutorials/StaticVonMisesTruss/staticVonMisesTruss.md', ...
              '~/.julia/dev/ONSAS_docs/docs/src/tutorials/CantileverBeam/cantileverBeam.md'} ;

for i=1:length(ONSASFiles)
  fprintf([ '  - ' ONSASFiles{i} '\n' ])
  m2md( ONSASFiles{i} , MDFiles{i},1)
end

%~ uniaxialExtension/onsasExample_uniaxialExtension.m
