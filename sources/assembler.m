% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro  
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

%
% function for assembly of:
%  - residual force vectors    (paramOut = 1 )
%  - residual tangent matrices (paramOut = 2)
%  - matrix with stresses      (paramOut = 3).
%

function Assembled = assembler ( Conec, crossSecsParamsMat, coordsElemsMat, ...
  materialsParamsMat, KS, Ut, paramOut, Udott, Udotdott, nodalDispDamping, ...
  solutionMethod, elementsParamsMat )

booleanCppAssembler = 0 ;

% -----------------------------------------------
% ---     C++ assembler                       ---
% -----------------------------------------------
if booleanCppAssembler
  CppAssembly

% -----------------------------------------------
% ---    octave/matlab assembler    ---
% -----------------------------------------------
else
  % -----------------------------------------------
  nElems  = size(Conec, 1) ;   nNodes  = length( Ut ) / 6 ;


  % ====================================================================
  %  --- 1 declarations ---
  % ====================================================================

  switch paramOut
  % -------  residual forces vector ------------------------------------
  case 1
    % --- creates Fint vector ---
    Fint = zeros( nNodes*6 , 1 ) ;
    Fmas = zeros( nNodes*6 , 1 ) ;
    Fvis = zeros( nNodes*6 , 1 ) ;

  % -------  tangent matrix        -------------------------------------
  case 2
    %~ if octaveBoolean
      %~ indsIK = uint32( zeros( nElems*24*24, 1 ) ) ;
      %~ indsJK = uint32( zeros( nElems*24*24, 1 ) ) ;
    indsIK =         zeros( nElems*24*24, 1 )   ;
    indsJK =         zeros( nElems*24*24, 1 )   ;
    valsK  =         zeros( nElems*24*24, 1 )   ;

    valsC  =         zeros( nElems*24*24, 1 )   ;
    valsM  =         zeros( nElems*24*24, 1 )   ;

    counterInds = 0 ; % counter non-zero indexes

  % -------  matrix with stress per element ----------------------------
  case 3
    stressMat = zeros( nElems, 6 ) ;
  end
  % ====================================================================


  % ====================================================================
  %  --- 2 loop assembly ---
  % ====================================================================

  for elem = 1:nElems


    % extract element properties
    elemMaterialParams     = materialsParamsMat( Conec( elem, 5) , : ) ;
    elemrho                = elemMaterialParams( 1     )              ;
    elemConstitutiveParams = elemMaterialParams( 2:end )             ;

    elemElementParams      = elementsParamsMat( Conec( elem, 6),: ) ;

    typeElem = elemElementParams(1) ;

    [numNodes, dofsStep] = elementTypeInfo ( typeElem ) ;

    % obtains nodes and dofs of element
    nodeselem   = Conec( elem, 1:numNodes )'      ;
    dofselem    = nodes2dofs( nodeselem , 6 )  ;
    dofselemRed = dofselem ( 1 : dofsStep : end ) ;

    elemDisps   = u2ElemDisps( Ut , dofselemRed ) ;

    elemCoords  = coordsElemsMat( elem, 1:dofsStep:( numNodes*6 ) ) ;

    if typeElem == 2 || typeElem == 3
      elemCrossSecParams     = crossSecsParamsMat ( Conec( elem, 4+4 ) , : ) ;
    end

    stress = [] ;

    switch typeElem

    % -----------   truss element   ------------------------------------
    case 2
      
      if elemCrossSecParams(1) == 1
        A = elemCrossSecParams(2) ;

      elseif elemCrossSecParams(1) == 2
        A = elemCrossSecParams(2)*elemCrossSecParams(3) ;

      elseif elemCrossSecParams(1) == 3
        A = pi*elemCrossSecParams(2)^2 / 4.0 ;

      else
        error('missing cross section case')
      end
      
      [ fs, ks, stress ] = elementTrussInternForce( elemCoords, elemDisps, elemConstitutiveParams, A, paramOut ) ;

      Finte = fs{1} ;
      Ke    = ks{1} ;

      if solutionMethod > 2
        booleanConsistentMassMat = elemElementParams(2) ;

        dotdotdispsElem  = u2ElemDisps( Udotdott , dofselem ) ;
        [ Fmase, Mmase ] = elementTrussMassForce( elemCoords, elemrho, A, elemElementParams(2), paramOut, dotdotdispsElem ) ;
        %
        Ce = zeros(size(Mmase));%+4.6e-4*Ke;
        %~ Fmase = fs{3} ;   Ce    = ks{2} ;   Mmase = ks{3} ;
        %~ Fmase = fs{3} ;   Ce    = ks{2} ;   Mmase = ks{3} ;
      end


    % -----------   frame element   ------------------------------------
    case 3

      [ fs, ks, stress ] = elementBeamForces( elemCoords, elemCrossSecParams, elemConstitutiveParams, solutionMethod,  u2ElemDisps( Ut       , dofselem ) , ...
                                               u2ElemDisps( Udott    , dofselem ) , ...
                                               u2ElemDisps( Udotdott , dofselem ), elemrho ) ;
      Finte = fs{1} ;  Ke    = ks{1} ;

      if solutionMethod > 2
        Fmase = fs{3} ;  Ce    = ks{2} ;   Mmase = ks{3} ;
      end


    % ---------  tetrahedron solid element -----------------------------
    case 4
      [ Finte, Ke, stress ] = elementTetraSolid( elemCoords, elemDisps, ...
                              elemConstitutiveParams, paramOut, elemElementParams(2)) ;

    end   % case tipo elemento
    % -------------------------------------------


    % -------------------------------------------
    % ---   assemble   ----
    % -------------------------------------------
    switch paramOut
    case 1
      % internal loads vector assembly
      Fint ( dofselemRed ) = Fint( dofselemRed ) + Finte ;
      if solutionMethod > 2
        Fmas ( dofselemRed ) = Fmas( dofselemRed ) + Fmase ;
      end

    case 2

      for indRow = 1:length( dofselemRed )

        entriesSparseStorVecs = counterInds + (1:length( dofselemRed) ) ;

        indsIK ( entriesSparseStorVecs  ) = dofselemRed( indRow )     ;
        indsJK ( entriesSparseStorVecs )  = dofselemRed       ;
        valsK  ( entriesSparseStorVecs )  = Ke( indRow, : )' ;

        if solutionMethod > 2
          valsM( entriesSparseStorVecs ) = Mmase( indRow, : )' ;
          if exist('Ce')~=0
            valsC( entriesSparseStorVecs ) = Ce   ( indRow, : )' ;
          end
        end

        counterInds = counterInds + length( dofselemRed ) ;
      end

    case 3
      StressVec( elem, (1:length(stress) ) ) = stress ;

    end % case paramOut ---

  end % for elements ----

  % ============================================================================




  % ============================================================================
  %  --- 3 global additions and output ---
  % ============================================================================

  Assembled = cell( 1 + ( ( solutionMethod > 2) && ( paramOut <= 2 ) ), 1 ) ;

  if solutionMethod > 2
    dampingMat          = sparse( nNodes*6, nNodes*6 ) ;
    dampingMat(1:2:end) = nodalDispDamping             ;
    dampingMat(2:2:end) = nodalDispDamping * 0.01      ;
  end

  switch paramOut

  case 1,

    Fint = Fint + KS * Ut ;

    Assembled{1} = Fint ;

    if solutionMethod > 2,
      Fvis = dampingMat * Udott ;
    end

    Assembled{2} = Fvis ;
    Assembled{3} = Fmas ;

  case 2,

    indsIK = indsIK(1:counterInds) ;
    indsJK = indsJK(1:counterInds) ;
    valsK  = valsK (1:counterInds) ;
    K      = sparse( indsIK, indsJK, valsK, size(KS,1), size(KS,1) ) + KS ;

    Assembled{1} = K ;

    if solutionMethod > 2
      valsM = valsM (1:counterInds) ;
      valsC = valsC (1:counterInds) ;
      M     = sparse( indsIK, indsJK, valsM , size(KS,1), size(KS,1) )  ;
      C     = sparse( indsIK, indsJK, valsC , size(KS,1), size(KS,1) ) + dampingMat ;
    else
      M = sparse(size( K ) ) ;
      C = sparse(size( K ) ) ;
    end

    Assembled{2} = C ;
    Assembled{3} = M ;

  case 3
    Assembled{1} = StressVec ;

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
