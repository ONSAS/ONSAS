nelems = size(Conec,1);

% --------------------------------------------------------------------
% --------------------------------------------------------------------
% write to file the following variables  

% write files
%~ varsInps = [ Ut; paramOut] ;

tiemposCppInterface = [];

auxt = cputime();

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

size(U)
save -ascii 'U.dat'                  U ;
save -ascii 'coordsElemsMat.dat'     coordsElemsMat;
save -ascii 'materialsParamsMat.dat' materialsParamsMat ;
save -ascii 'elementsParamsMat.dat'  elementsParamsMat  ;
save -ascii 'crossSecsParamsMat.dat' crossSecsParamsMat ;
save -ascii 'coordsElemsMat.dat' coordsElemsMat ;

save -ascii 'constantFext.dat' constantFext;
save -ascii 'variableFext.dat' variableFext;

scalarParams = [ currLoadFactor nextLoadFactor nodalDispDamping currTime timeIndex ]' ;

save -ascii 'scalarParams.dat' scalarParams;
save -ascii 'neumdofs.dat' neumdofs;

sfid = fopen('strings.txt','w') ;
fprintf( sfid, [ outputDir '\n' problemName '\n'] ) ;
fclose( sfid );


tiemposCppInterface(1) = cputime() - auxt ;
auxt = cputime() ;

% --------------------------------------------------------------------
% --------------------------------------------------------------------

% run sts
%~ [status, output] = system('../sources/timeStepIteration.lnx') ;
%~ output
auxPathBin = file_in_loadpath('timeStepIteration.lnx')

[status] = system( auxPathBin ,0) ;

tiemposCppInterface(2) = cputime() - auxt ;
auxt = cputime() ;

% read file
Ut         = load( '-ascii', 'Ut.dat'        ) ;
Utp1       = load( '-ascii', 'Utp1.dat'      ) ;
Udottp1    = load( '-ascii', 'Udottp1.dat'   ) ;
Udotdottp1 = load( '-ascii', 'Udotdottp1.dat') ;

Stresstp1 = zeros( size( Stress));

auxOutValsVec = load( '-ascii', 'auxOutValsVec.dat') ;

nextTime    = auxOutValsVec(1) ;
stopCritPar = auxOutValsVec(2) ;
dispIters   = auxOutValsVec(3) ;
solutionMethod   = auxOutValsVec(4) ;

auxspmat = load('systemDeltauMatrixCpp.dat');
systemDeltauMatrix = sparse( auxspmat(:,1)+1, auxspmat(:,2)+1, auxspmat(:,3) ) ;


tiemposCppInterface(3) = cputime() - auxt 


