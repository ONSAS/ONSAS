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


%function for assembly of tangent matrices and/or forces vectors.
%
% Inputs:
%   paramOut: parameter used to set output:
%     only internal forces vector (1) or only tangent matrices (2)  
%

function [ Assembled, StressVec ] = assembler ( Conec, crossSecsParams, coordsElemsMat, materialsParamsMat, KS, Ut, paramOut, Udotdott, booleanConsistentMassMat )

booleanCppAssembler = 0 ;

% -----------------------------------------------
% ---     C++ assembler    ---
% -----------------------------------------------
if booleanCppAssembler
  %~ CppAssembly
  
% -----------------------------------------------
% ---    octave assembler    ---
% -----------------------------------------------
else
  % -----------------------------------------------
  nelems  = size(Conec,1) ;   nnodes  = length( Ut) / 6 ;
  
  %~ profile clear, profile on

  Assembled = cell( 2, 1 ) ;
    
  % --- creates Fint vector ---
  if paramOut == 1
    Fintt = zeros( nnodes*6 , 1 ) ;
  
    if nargout == 2
      StressVec   = zeros( nelems, 6 ) ;
    end
    
    Fmast = zeros( nnodes*6 , 1 ) ;
  
  elseif paramOut == 2
  
    % assumes maximum 4 nodes per element
    indsIKT = uint32( zeros( nelems*24*24, 1 ) ) ;
    indsJKT = uint32( zeros( nelems*24*24, 1 ) ) ;
    valsKT  =         zeros( nelems*24*24, 1 )   ;
    valsMT  =         zeros( nelems*24*24, 1 )   ;
    
    counterInds = 0 ;
  
  end
  
  % ----------------------------------------------
  
  contTiempoLlamadasIndexs       = 0;
  contTiempoLlamadasAssembly     = 0;
  contTiempoLlamadasAssemblyFint = 0;


  %~ if solutionMethod == 3 || solutionMethod == 4
  %~ Udotdottp1
    %~ Fine    = massMat * Udotdottp1 ;
    %~ Finered = Fine( neumdofs ) ;
  %~ else, Finered    = [] ; end


  % ----------------------------------------------
  % loop for assembly
  for elem = 1:nelems

    elemCrossSecParams = crossSecsParams( Conec( elem, 6 ) , : ) ;

    elemMaterialParams     = materialsParamsMat( Conec( elem, 5), : ) ;
    
    elemrho                = elemMaterialParams( 1     )              ;
    elemConstitutiveParams = elemMaterialParams( 2:end )              ;

    switch Conec(elem,7)
  
    % -------------------------------------------
    case 1 % Co-rotational Truss with Engineering strain
  
      % obtains nodes and dofs of element
      nodeselem = Conec(elem,1:2)'             ;
      dofselem  = nodes2dofs( nodeselem , 6 )  ;
      dispsElem = u2ElemDisps( Ut , dofselem ) ;
  
      dofselemRed = dofselem(1:2:end)  ;
        
      A = elemCrossSecParams(1) ;
      
      sizeTensor = 1 ;
      
      [ Finte, KTe, stress, dstressdeps, strain ] = elementTrussInternForce( coordsElemsMat(elem,1:12)', dispsElem, elemConstitutiveParams, A, paramOut ) ;
       
      if elemrho > 0
        dotdotdispsElem  = u2ElemDisps( Udotdott , dofselem ) ;
        
        [ Fmase, Mmase ] = elementTrussMassForce( coordsElemsMat(elem,1:12)', elemrho, A, booleanConsistentMassMat, paramOut, dotdotdispsElem  );
      end
      
    % -------------------------------------------
    case 2 % Co-rotational Frame element (bernoulli beam)
  
      % obtains nodes and dofs of element
      nodeselem = Conec(elem,1:2)'             ;  
      dofselem  = nodes2dofs( nodeselem , 6 )  ;
      dispsElem = u2ElemDisps( Ut , dofselem ) ;

      dofselemRed = dofselem  ;
  
      sizeTensor = 1 ;
  
      A   = elemCrossSecParams(Conec(elem,6),1) ;
      Iyy = elemCrossSecParams(Conec(elem,6),2) ;
      Izz = elemCrossSecParams(Conec(elem,6),3) ;
      J   = elemCrossSecParams(Conec(elem,6),4) ;
  
      xs = coordsElemsMat(elem,1:2:end)'        ;
      E  = elemConstitutiveParams(2) ;
      nu = elemConstitutiveParams(3) ;
      G  = E/(2*(1+nu)) ;
      
      [ Finte, KTe, strain, stress ]= elementBeamInternLoads( xs, dispsElem , [E G A Iyy Izz J] ) ;
  
      if elemrho > 0
        [Fmase,~,~] = elementFuerzaInercial(xs, Dte, Ddote, Ddotdote, params,Jrho );
      end

  
    % -------------------------------------------
    case 3 % linear solid element
      
      
      % obtains nodes and dofs of element
      nodeselem = Conec(elem,1:4)' ;
      dofselem  = nodes2dofs( nodeselem , 6 ) ;
      dofselemRed = dofselem(1:2:end) ;

      dispsElem = u2ElemDisps( Ut , dofselemRed ) ;
         
      tetcoordmat        = zeros(3,4) ;
      tetcoordmat(1,1:4) = coordsElemsMat(elem,1:6:end) ;
      tetcoordmat(2,1:4) = coordsElemsMat(elem,3:6:end) ;
      tetcoordmat(3,1:4) = coordsElemsMat(elem,5:6:end) ;
      
      sizeTensor = 6 ;
  
      if elemMaterialParams( 1 ) == 6
        E  = elemMaterialParams( 2) ;
        nu = elemMaterialParams( 3) ;
  
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
        E  = elemMaterialParams( 2) ;
        nu = elemMaterialParams( 3) ;
        [ Finte, KTe, strain, stress ]= elementTetraSolidInternLoadsTangMat ( tetcoordmat, dispsElem , [E nu], paramOut ) ;
      end
  
    end   % case tipo elemento
    % -------------------------------------------


  


  
  
    % -------------------------------------------
    % ---   assemble   ----
    % -------------------------------------------
    if paramOut == 1
      % internal loads vector assembly
      Fintt ( dofselemRed ) = Fintt( dofselemRed ) + Finte ;
      %~ FintGt ( dofstet ) = FintGt( dofstet ) + Finte ;
      
      if elemrho > 0
      elemrho
        Fmast ( dofselemRed ) = Fmast( dofselemRed ) + Fmase ;
      end
        
      if nargout == 3
        %~ StrainVec(elem,(1:sizeTensor) ) = strain ;
        StressVec(elem,(1:sizeTensor) ) = stress ;
      end
      
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

        if elemrho > 0
          valsMT( entriesSparseStorVecs ) = Mmase( indRow, : )' ;
        end
        
        counterInds = counterInds + length( dofselemRed ) ;
      end
    
    end % if paramOut

  end % for elements ----
    
  
  if paramOut == 1,

    Fintt = Fintt + KS * Ut ;

    Assembled{1} = Fintt ;
    Assembled{2} = Fmast ;

  elseif paramOut == 2,
  
    indsIKT = indsIKT(1:counterInds) ;
    indsJKT = indsJKT(1:counterInds) ;
    valsKT  = valsKT (1:counterInds) ;
    KT      = sparse( indsIKT, indsJKT, valsKT, size(KS,1), size(KS,1) ) + KS ;

    Assembled{1} = KT ;
    
    valsMT  = valsMT (1:counterInds) ;
    MT      = sparse( indsIKT, indsJKT, valsMT, size(KS,1), size(KS,1) )      ;
    Assembled{2} = MT ;
    
  end
  
end % if booleanCppAssembler
% ----------------------------------------


% ==============================================================================
%
%
% ==============================================================================

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

% ==============================================================================
%
% function to convert vector of displacements into displacements of element.
%
% ==============================================================================
% _____&&&&&&&&&&&&&& GENERALIZAR PARA RELEASES &&&&&&&&&&&&&&&
function elemDisps = u2ElemDisps( U, dofselem)

elemDisps = U(dofselem);
