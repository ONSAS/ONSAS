  timer = time();
  nelems    = size(Conec,1);

  currDir = pwd ;
  cd( '~/work/repos/inomiso/onsaspp/src' )

  % --------------------------------------------------------------------
  % --------------------------------------------------------------------
  % write to file the following variables  
  
  % write files
  varsInps = [ Ut; paramOut] ;
  
  save -ascii 'varsInps.dat' varsInps;
  save -ascii 'Conec.dat' Conec;
  save -ascii 'coordsElemsMat.dat' coordsElemsMat;
  save -ascii 'hyperElasParamsMat.dat' hyperElasParamsMat;

  timeEscritaArchivos = time() - timer

  % --------------------------------------------------------------------
  % --------------------------------------------------------------------
  
  % run sts
  timer = time();
  [status, output] = system('./onsasAssembler')
  %~ system('g++ onsasAssembler.cpp -larmadillo -o onsasAssembler')
  %~ system('make')
  %~ system('./onsasAssembler > salidaaa.txt');
  timeEnsamblado = time() - timer
  
  
  % read file
  FintGt = load( '-ascii', 'FintGt.dat') ;
  
  if paramOut==2

  timer = time();
    
    indsIKT = load( '-ascii', 'indsIKT.dat') ;
    indsJKT = load( '-ascii', 'indsJKT.dat') ;
    valsKT  = load( '-ascii', 'valsIKT.dat') ;
  timeLectura = time() - timer
    
    KT = ...
      sparse( indsIKT, indsJKT, valsKT, size(KS,1), size(KS,1) ) ...
      + KS ;
  end

  %~ issparse(KT)
  timeLLamadaCpp = time() - timer
  
  % ------------------
  %~ if stopit
    %~ cd([currDir '/examples']), stop
  %~ else
    cd([currDir ])
  %~ end

  StrainVec   = sparse( nelems, 6 ) ;
  StressVec   = sparse( nelems, 6 ) ;



