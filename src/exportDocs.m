
addpath( genpath( '../../../m2md_repo/' ) )

pathToDocs = '~/.julia/dev/ONSAS_docs/' ;

m2md( '../examples/staticVonMisesTruss/onsasExample_staticVonMisesTruss.m', [ pathToDocs 'docs/src/tutorials/StaticVonMisesTruss/staticVonMisesTruss.md' ], 1 )

m2md( '../examples/heatAnalytic/example_HeatAnalytic.m', [ pathToDocs 'docs/src/tutorials/HeatAnalytic/HeatAnalytic.md' ], 1 )

%~ cd( pathToDocs ) ;
%~ [~,out] = system( ['cd ' pathToDocs ' \\ ls \\ ./preview.sh' ] ) ;
%~ disp(out)
