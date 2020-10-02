

nelems = size(Conec,1);

% --------------------------------------------------------------------
% --------------------------------------------------------------------
% write to file the following variables  

% write files
%~ varsInps = [ Ut; paramOut] ;

%~ save -ascii 'varsInps.dat' varsInps;
save -ascii 'Conec.dat' Conec;
save -ascii  'numericalMethodParams.dat' numericalMethodParams;

save  'systemDeltauMatrix.dat' systemDeltauMatrix ;
status = system('tail -n +7 systemDeltauMatrix.dat > aux.dat' ); 
status = system('mv aux.dat systemDeltauMatrix.dat' ) ; 

save -ascii 'U.dat'                  U ;
save -ascii 'materialsParamsMat.dat' materialsParamsMat;
save -ascii 'coordsElemsMat.dat'     coordsElemsMat;
save -ascii 'crossSecsParamsMat.dat' crossSecsParamsMat;
% --------------------------------------------------------------------
% --------------------------------------------------------------------

% run sts
%~ [status, output] = system('../sources/timeStepIteration.lnx') ;
%~ output

[status] = system('../sources/timeStepIteration.lnx',0) ;

stop
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



