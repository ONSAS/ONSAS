function [twist, chord, IDsection, radius, Idx1, Idx2] = uBEMframeProps(uBEMdataCoords, xs)

Idx1      = find(uBEMdataCoords(1:end, 1:3) == xs(1:3)', 1) ;
Idx2      = find(uBEMdataCoords(1:end, 1:3) == xs(4:6)', 1) ;

radius    = [ uBEMdataCoords(Idx1, 7), uBEMdataCoords(Idx2, 7) ];

twist     = [ 0 uBEMdataCoords(Idx1, 4)            0             0 uBEMdataCoords(Idx2, 4)              0           ]';
chord     = [ 0             0           uBEMdataCoords(Idx1, 5)  0             0            uBEMdataCoords(Idx2, 5) ]';
IDsection = [ uBEMdataCoords(Idx1, 6), uBEMdataCoords(Idx2, 6) ]; 

end