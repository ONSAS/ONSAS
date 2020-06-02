% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.


%function for assembly of tangent stiffness matrix and/or internal forces vector.
%
% Inputs:
%   paramOut: parameter used to set output: only internal forces vector (1) or only tangent matrices (2)  
%

function [FintGt, KT, StrainVec, StressVec ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut, bendStiff, paramOut )

booleanCppAssembler = 0 ;

% -----------------------------------------------
% -----------------------------------------------
% C++ assembler
if booleanCppAssembler

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




% -----------------------------------------------
% -----------------------------------------------
% --- octave assembler ---
else

  % -----------------------------------------------
  nelems    = size(Conec,1)   ;
  nnodes    = length( Ut) / 6 ;
  
  %~ profile clear, profile on
  
  % creates Fint vector
  FintGt = zeros( nnodes*6 , 1 ) ;
  
  % assumes maximum 4 nodes per element
  indsIKT = uint32( zeros( nelems*24*24, 1 ) ) ;
  indsJKT = uint32( zeros( nelems*24*24, 1 ) ) ;
  valsKT  =         zeros( nelems*24*24, 1 ) ;

  counterInds = 0 ;
  
  StrainVec   = zeros( nelems, 6 ) ;
  StressVec   = zeros( nelems, 6 ) ;
  % ----------------------------------------------
  
  contTiempoLlamadasIndexs = 0;
  contTiempoLlamadasAssembly = 0;
  contTiempoLlamadasAssemblyFint = 0;

  % ----------------------------------------------
  % loop for assembly
  for elem = 1:nelems
  
    switch Conec(elem,7)
  
    % -------------------------------------------
    case 1 % Co-rotational Truss with Engineering strain
  
      % obtains nodes and dofs of element
      nodeselem = Conec(elem,1:2)' ;
      dofselem  = nodes2dofs( nodeselem , 6 ) ;
      dispsElem = u2ElemDisps( Ut , dofselem ) ;
  
      dofselemRed = dofselem(1:2:end)  ;
      
      sizeTensor = 1 ;
  
      A  = secGeomProps(Conec(elem,6),1) ;
      hyperAux  = hyperElasParamsMat( Conec(elem,5),:) ;
      
      [ Finte, KTe, stress, dstressdeps, strain ] = elementTrussEngStr( coordsElemsMat(elem,1:12)', dispsElem, hyperAux , A, paramOut ) ;
      
      %~ size(Finte)
      %~ size(KTe)
      %~ stop  
    % -------------------------------------------
    case 2 % Co-rotational Frame element (bernoulli beam)
  
      % obtains nodes and dofs of element
      nodeselem = Conec(elem,1:2)' ;
      dofselem  = nodes2dofs( nodeselem , 6 ) ;
      dispsElem = u2ElemDisps( Ut , dofselem ) ;

      dofselemRed = dofselem  ;
  
      sizeTensor = 1 ;
  
      A   = secGeomProps(Conec(elem,6),1) ;
      Iyy = secGeomProps(Conec(elem,6),2) ;
      Izz = secGeomProps(Conec(elem,6),3) ;
      J   = secGeomProps(Conec(elem,6),4) ;
  
      xs = coordsElemsMat(elem,1:2:end)'        ;
      E  = hyperElasParamsMat( Conec(elem,5),2) ;
      nu = hyperElasParamsMat( Conec(elem,5),3) ;
      G  = E/(2*(1+nu)) ;
  
      [ Finte, KTe, strain, stress ]= elementBeam3DInternLoads( xs, dispsElem , [E G A Iyy Izz J] ) ;
  
      KL0e = KTe;
  
  
    % -------------------------------------------
    case 3 % linear solid element
      
      
      % obtains nodes and dofs of element
      nodeselem = Conec(elem,1:4)' ;
      dofselem  = nodes2dofs( nodeselem , 6 ) ;
      dispsElem = u2ElemDisps( Ut , dofstet ) ;
     
      dofselemRed = dofselem(1:2:length(dofselem)) ;
      
      %~ dofselem = dofstet ;
   
      tetcoordmat        = zeros(3,4) ;
      tetcoordmat(1,1:4) = coordsElemsMat(elem,1:6:end) ;
      tetcoordmat(2,1:4) = coordsElemsMat(elem,3:6:end) ;
      tetcoordmat(3,1:4) = coordsElemsMat(elem,5:6:end) ;
      
      sizeTensor = 6 ;
  
      if hyperElasParamsMat( Conec(elem,5), 1 ) == 6
        E  = hyperElasParamsMat( Conec(elem,5),2) ;
        nu = hyperElasParamsMat( Conec(elem,5),3) ;
  
        if paramOut==1
        
          iniAss = time() ;
  
          [ Finte ] = elementTetraSVKSolidInternLoadsTangMat( tetcoordmat, dispsElem , [E nu], paramOut ) ; 
          
          contTiempoLlamadasAssemblyFint = contTiempoLlamadasAssemblyFint + ( time() - iniAss) ;
          
          strain = zeros(6,1) ;
          stress = zeros(6,1) ;
  
        elseif paramOut == 2
  
          iniAss = time() ;
  
          [ Finte, KTe, strain, stress ] = elementTetraSVKSolidInternLoadsTangMat ( tetcoordmat, dispsElem , [E nu] , paramOut) ;
  
          contTiempoLlamadasAssembly = contTiempoLlamadasAssembly + ( time() - iniAss) ;
  
        end
        
      else
        E  = hyperElasParamsMat( Conec(elem,5),2) ;
        nu = hyperElasParamsMat( Conec(elem,5),3) ;
        [ Finte, KTe, strain, stress ]= elementTetraSolidInternLoadsTangMat ( tetcoordmat, dispsElem , [E nu], paramOut ) ;
      end
  
    end   % case tipo elemento
    % -------------------------------------------
  
  
  
    % -------------------------------------------
    if paramOut == 1
      % internal loads vector assembly
      FintGt ( dofselemRed ) = FintGt( dofselemRed ) + Finte ;
      %~ FintGt ( dofstet ) = FintGt( dofstet ) + Finte ;
    
      StrainVec(elem,(1:sizeTensor) ) = strain ;
      StressVec(elem,(1:sizeTensor) ) = stress ;
    
    elseif paramOut == 2
      % matrices assembly
      %~ KT  (dofselem,dofselem) = KT(dofselem,dofselem) + KTe     ;
    %~ else
      for indRow = 1:length( dofselemRed )
  
        %~ indVec = (indRow+1)/2 ;
      
        %~ entriesSparseStorVecs = (elem-1)*24*24 + (indRow-1) * 24 + (1:24) ;
        entriesSparseStorVecs = counterInds + (1:length( dofselemRed) ) ;        
        
        indsIKT ( entriesSparseStorVecs  ) = dofselemRed( indRow )     ;
        %~ indsJKT ( entriesSparseStorVecs ) = dofselem            ;
        %~ valsKT  ( entriesSparseStorVecs ) = KTe( indRow, : ) ;
  
        indsJKT ( entriesSparseStorVecs ) = dofselemRed       ;
        valsKT  ( entriesSparseStorVecs ) = KTe( indRow, : )' ;
        
        counterInds = counterInds + length( dofselemRed ) ;

      end
    
    end % if paramout
    % ---------------------
   
  end % for elements ----
  
  %~ loopelemtie = time() - chetiem ;
  
  indsIKT = indsIKT(1:counterInds) ;
  indsJKT = indsJKT(1:counterInds) ;
  valsKT  = valsKT (1:counterInds) ;
  
  if paramOut == 2
    KT     = sparse( indsIKT, indsJKT, valsKT, size(KS,1), size(KS,1) ) + KS ;
  end
  
  FintGt = FintGt + KS * Ut ;
  
  if length(bendStiff) >0
  
    Nodes = conv ( Conec, coordsElemsMat+dispsElemsMat ) ;
  
    [ ~, KTAngSpr ] = loadsAngleSpring( Nodes, Conec, bendStiff ) ;
  
    fextAngSpr = KTAngSpr*Ut ;
  
    KT     += sparse(KTAngSpr)   ;
    FintGt += fextAngSpr ;
  
  end % if bend stiff ----
  
  %~ fintiem = time() - chetiem;
  % ------------------------------------

end % if booleanCppAssembler
% ----------------------------------------




%~ KTsparse = sparse( indsIKT, indsJKT, valsKT ) ;

%~ KT     = KT  + KS ;
%~ FintGt = FintGt + KS*Ut ;


%~ if length(bendStiff) >0

  %~ Nodes = conv ( Conec, coordsElemsMat+dispsElemsMat ) ;

  %~ [ ~, KTAngSpr ] = loadsAngleSpring( Nodes, Conec, bendStiff ) ;

  %~ fextAngSpr = KTAngSpr*Ut ;

  %~ KT     = KT     + KTAngSpr   ;
  %~ FintGt = FintGt + fextAngSpr ;

%~ end

%~ % ------------------------------------
%~ >>>>>>> 85df745bf8cc84eb567a52786a477589e9d8673e





function nodesmat = conv ( conec, coordsElemsMat ) 
nodesmat  = [] ;
nodesread = [] ;

for i=1:size(conec,1)
  for j=1:2
    if length( find( nodesread == conec(i,j) ) ) == 0
      nodesmat( conec(i,j),:) = coordsElemsMat( i, (j-1)*6+(1:2:5) ) ;
    end
  end
end

% ---------------------
% function to convert vector of displacements into displacements of element.
% ---------------------
function elemDisps = u2ElemDisps( U, dofselem)

elemDisps = U(dofselem);
