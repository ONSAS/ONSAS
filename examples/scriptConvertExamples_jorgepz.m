
addpath(genpath('../src/'));

cantileverBeamONSASFile = './uniformCurvatureCantilever/onsasExample_uniformCurvatureCantilever_test.m' ;

cantileverBeamMDFile = '~/.julia/dev/ONSAS_docs/docs/src/tutorials/CantileverBeam/cantileverBeam.md' ;


m2md( cantileverBeamONSASFile , cantileverBeamMDFile,1)
