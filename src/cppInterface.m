
%mdThe goal of this function is to convert the input data to a text file format, run ONSAS.cpp and get the results

function [ matUs, loadFactorsMat ] = cppInterface( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams )

onsInputFile = [ otherParams.problemName '.ons' ] ;
fidTextInputFile = fopen( onsInputFile, 'w' ) ;

% --- write materials ---
numMaterials = length( materials.hyperElasModel );
fprintf( fidTextInputFile , 'materials\n')
fprintf( fidTextInputFile , 'numMaterials %03i\n', numMaterials )
fprintf( fidTextInputFile , 'materials.hyperElasModel\n')
for i=1:numMaterials
  fprintf( fidTextInputFile , '%s\n', materials.hyperElasModel{i} )
end
fprintf( fidTextInputFile , 'materials.hyperElasParams\n')
for i=1:numMaterials
  fprintf( fidTextInputFile , '%15.10e ', materials.hyperElasParams{i} )
end
fprintf( fidTextInputFile , '\n')

fprintf( fidTextInputFile , '----------------\n') % separator

% --- write elements ---
numElements = length( elements.elemType )
fprintf( fidTextInputFile , 'elements\n')
fprintf( fidTextInputFile , 'numElements %03i\n', numElements )
fprintf( fidTextInputFile , 'elements.elemType\n')
for i=1:numElements
  fprintf( fidTextInputFile , '%s\n', elements.elemType{i} )
end
fprintf( fidTextInputFile , 'elements.elemTypeParams\n')
for i=1:numElements
  if isempty( elements.elemTypeParams{i} )
    fprintf( fidTextInputFile , '\n' )
  end
end
fprintf( fidTextInputFile , 'elements.elemTypeGeometry\n')
for i=1:numElements
  if isempty( elements.elemTypeGeometry{i} )
    fprintf( fidTextInputFile , '\n' )
  end
end

%md close and delete the ons file
fclose( fidTextInputFile ) ;


%~ [status, output] = system('../sources/timeStepIteration.lnx') ;
%~ output
auxPathBin = file_in_loadpath('ONSAS.lnx');
[ auxPathBin ' ' onsInputFile ]
[ status ] = system( [ auxPathBin ' ' onsInputFile ] , 0 ) ;


% delete( [ otherParams.problemName '.ons' ] ) ;
stop





Conec = modelProperties.Conec

nelems = size(Conec,1) ;

% --------------------------------------------------------------------
% --------------------------------------------------------------------
% write to file the following variables

% write files
%~ varsInps = [ Ut; paramOut] ;

tiemposCppInterface = [];

auxt = cputime();

analy = modelProperties.analysisSettings
save -ascii 'testAnalysis.dat' analy

stop

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
