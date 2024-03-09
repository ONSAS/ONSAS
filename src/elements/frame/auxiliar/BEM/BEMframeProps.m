function [IDsection, radius, chords, beta, nodes] = BEMframeProps(BEMdataList, xs)

Idx1      = find( all( BEMdataList(1:end, 1:3) == xs(1:3)', 2) ) ;
Idx2      = find( all( BEMdataList(1:end, 1:3) == xs(4:6)', 2) ) ;

nodes = [ Idx1, Idx2 ];

radius    = [ BEMdataList(Idx1, 4)', BEMdataList(Idx2, 4)' ];
beta      = [ [BEMdataList(Idx1, 5), 0, 0 ]'  , [BEMdataList(Idx2, 5), 0, 0]'   ];
chords    = [ [0, -BEMdataList(Idx1, 6), 0]'  , [0, -BEMdataList(Idx2, 6), 0]'  ];
IDsection = [ BEMdataList(Idx1, 7)', BEMdataList(Idx2, 7)' ];



end