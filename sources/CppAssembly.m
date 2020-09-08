nelems = size(Conec,1);

% --------------------------------------------------------------------
% --------------------------------------------------------------------
% write to file the following variables  

% write files
varsInps = [ Ut; paramOut] ;

save -ascii 'varsInps.dat' varsInps;
save -ascii 'Conec.dat' Conec;
save -ascii 'coordsElemsMat.dat' coordsElemsMat;
save -ascii 'materialsParamsMat.dat' materialsParamsMat;

% --------------------------------------------------------------------
% --------------------------------------------------------------------

% run sts
[status, output] = system('./sources/onsasAssembler.lnx') ;

% read file
FintGt = load( '-ascii', 'FintGt.dat') ;

if paramOut==2
  
  indsIKT = load( '-ascii', 'indsIKT.dat') ;
  indsJKT = load( '-ascii', 'indsJKT.dat') ;
  valsKT  = load( '-ascii', 'valsIKT.dat') ;
  
  KT = ...
    sparse( indsIKT, indsJKT, valsKT, size(KS,1), size(KS,1) ) ...
    + KS ;
end

StrainVec   = sparse( nelems, 6 ) ;
StressVec   = sparse( nelems, 6 ) ;



