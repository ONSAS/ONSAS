
nelems = size(Conec,1);

% --------------------------------------------------------------------
% --------------------------------------------------------------------
% write to file the following variables  

% write files
%~ varsInps = [ Ut; paramOut] ;

%~ save -ascii 'varsInps.dat' varsInps;
save -ascii 'Conec.dat' Conec ;
numericalMethodParamsT = numericalMethodParams' ;
save -ascii  'numericalMethodParams.dat' numericalMethodParamsT ;

save  'systemDeltauMatrix.dat' systemDeltauMatrix ;
status = system('tail -n +7 systemDeltauMatrix.dat > aux.dat' ); 
status = system('mv aux.dat systemDeltauMatrix.dat' ) ; 

save  'KS.dat' KS ;
status = system('tail -n +7 KS.dat > aux.dat' ); 
status = system('mv aux.dat KS.dat' ) ;

save -ascii 'U.dat'                  U ;
save -ascii 'coordsElemsMat.dat'     coordsElemsMat;
save -ascii 'materialsParamsMat.dat' materialsParamsMat ;
save -ascii 'elementsParamsMat.dat'  elementsParamsMat  ;
save -ascii 'crossSecsParamsMat.dat' crossSecsParamsMat ;
save -ascii 'coordsElemsMat.dat' coordsElemsMat ;

save -ascii 'constantFext.dat' constantFext;
save -ascii 'variableFext.dat' variableFext;
currTime
scalarParams = [ currLoadFactor nextLoadFactor nodalDispDamping currTime ]' ;

neumdofs

save -ascii 'scalarParams.dat' scalarParams;
save -ascii 'neumdofs.dat' neumdofs;
% --------------------------------------------------------------------
% --------------------------------------------------------------------

% run sts
%~ [status, output] = system('../sources/timeStepIteration.lnx') ;
%~ output

[status] = system('../sources/timeStepIteration.lnx',0) ;

% read file
Utp1 = load( '-ascii', 'Utp1.dat') ;


