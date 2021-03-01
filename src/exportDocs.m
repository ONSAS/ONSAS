
addpath( genpath( '../../../m2md_repo/' ) )

pathToDocs = '~/.julia/dev/ONSAS_docs/' ;

m2md( '../examples/staticVonMisesTruss/onsasExample_staticVonMisesTruss.m', [ pathToDocs 'docs/src/tutorials/StaticVonMisesTruss/staticVonMisesTruss.md' ], 1 )
