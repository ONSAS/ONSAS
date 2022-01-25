% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
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

%mdThis function computes the assembled force vectors, tangent matrices and stress matrices.
function [ fsCell, stressMat, tangMatsCell ] = assembler ( Conec, elements, Nodes,
                                                           materials, KS, Ut, Udott, Udotdott,
                                                           analysisSettings, outputBooleans, nodalDispDamping,
                                                           timeVar )

fsBool     = outputBooleans(1) ; stressBool = outputBooleans(2) ; tangBool   = outputBooleans(3) ;

nElems     = size(Conec, 1) ;
nNodes     = size(Nodes, 1) ;
% ====================================================================
%  --- 1 declarations ---
% ====================================================================

% -------  residual forces vector ------------------------------------
if fsBool
  % --- creates Fint vector ---
  Fint  = zeros( nNodes*6 , 1 ) ;
  Fmas  = zeros( nNodes*6 , 1 ) ;
  Fvis  = zeros( nNodes*6 , 1 ) ;
  Faero = zeros( nNodes*6 , 1 ) ;
end

% -------  tangent matrix        -------------------------------------
if tangBool
  indsIK =         zeros( nElems*24*24, 1 )   ;
  indsJK =         zeros( nElems*24*24, 1 )   ;
  valsK  =         zeros( nElems*24*24, 1 )   ;

  valsC  =         zeros( nElems*24*24, 1 )   ;
  valsM  =         zeros( nElems*24*24, 1 )   ;

  counterInds = 0 ; % counter non-zero indexes
end

% -------  matrix with stress per element ----------------------------
if stressBool
  stressMat = zeros( nElems, 6 ) ;
else
  stressMat = [] ;
end
% ====================================================================


dynamicProblemBool = strcmp( analysisSettings.methodName, 'newmark' ) || strcmp( analysisSettings.methodName, 'alphaHHT' ) ;

% ====================================================================
%  --- 2 loop assembly ---
% ====================================================================

for elem = 1:nElems
  mebiVec = Conec( elem, 1:4) ;

  %md extract element properties
  hyperElasModel   = materials( mebiVec( 1 ) ).hyperElasModel  ;
  hyperElasParams  = materials( mebiVec( 1 ) ).hyperElasParams ;
  density          = materials( mebiVec( 1 ) ).density         ;

  elemType         = elements( mebiVec( 2 ) ).elemType         ;
  elemTypeParams   = elements( mebiVec( 2 ) ).elemTypeParams   ;
  elemTypeGeometry = elements( mebiVec( 2 ) ).elemTypeGeometry ;

  %md extract aerodinamic properties
  elemTypeAero     = elements( mebiVec( 2 ) ).elemTypeAero     ;
  userDragCoef     = elements( mebiVec( 2 ) ).userDragCoef     ;
  userLiftCoef     = elements( mebiVec( 2 ) ).userLiftCoef     ;
  userMomentCoef   = elements( mebiVec( 2 ) ).userMomentCoef   ;

  %md compute aerodynamic compute force boolean
  AeroCoefficentsBool = ~isempty( userDragCoef ) || ~isempty( userMomentCoef ) || ~isempty( userLiftCoef ) ;
  aeroBool = ~isempty( elemTypeAero ) && AeroCoefficentsBool ;
  %md chcek unless one coefficient is defined

  %md obtain elemeny info
  [numNodes, dofsStep] = elementTypeInfo ( elemType ) ;

  %md obtains nodes and dofs of element
  nodeselem   = Conec( elem, (4+1):(4+numNodes) )' ;
  dofselem    = nodes2dofs( nodeselem , 6 )        ;
  dofselemRed = dofselem ( 1 : dofsStep : end )    ;

  %md elemDisps contains the displacements corresponding to the dofs of the element
  elemDisps   = u2ElemDisps( Ut , dofselemRed ) ;

  elemNodesxyzRefCoords  = reshape( Nodes( Conec( elem, (4+1):(4+numNodes) )' , : )',1,3*numNodes) ;

  stressElem = [] ;

  % -----------   truss element   ------------------------------
  if strcmp( elemType, 'truss')

    A  = crossSectionProps ( elemTypeGeometry, density ) ;

    [ fs, ks, stressElem ] = elementTrussInternForce( elemNodesxyzRefCoords, elemDisps, hyperElasModel, hyperElasParams, A ) ;

    Finte = fs{1} ;    Ke    = ks{1} ;

    if dynamicProblemBool
      booleanConsistentMassMat = elemTypeParams(1) ;

      dotdotdispsElem  = u2ElemDisps( Udotdott , dofselem ) ;
      [ Fmase, Mmase ] = elementTrussMassForce( elemNodesxyzRefCoords, density, A, booleanConsistentMassMat, dotdotdispsElem ) ;
      %
      Ce = zeros( size( Mmase ) ) ; % only global damping considered (assembled after elements loop)
    end


  % -----------   frame element   ------------------------------------
  elseif strcmp( elemType, 'frame')

		if strcmp(hyperElasModel, 'linearElastic')

			[ Finte, Ke ] = linearStiffMatBeam3D(elemNodesxyzRefCoords, elemTypeGeometry, density, hyperElasParams, elemDisps ) ;

		elseif strcmp( hyperElasModel, '1DrotEngStrain')

      [ fs, ks, stressElem ] = elementBeamForces( elemNodesxyzRefCoords, elemTypeGeometry, [ 1 hyperElasParams ], u2ElemDisps( Ut       , dofselem ) , ...
                                               u2ElemDisps( Udott    , dofselem ) , ...
                                               u2ElemDisps( Udotdott , dofselem ), density ) ;
      Finte = fs{1} ;  Ke = ks{1} ;

      if dynamicProblemBool
        Fmase = fs{3} ;Ce = ks{2} ; Mmase = ks{3} ;
      end
    else
      error('wrong hyperElasModel for frame element.')
   end
    if ~isempty(elemTypeAero) &&  aeroBool == 0
      error("Drag, Lift or Moment coefficients must be defined in elements struct\n ")
    end
    %md chcek wind velocity is defined
    if AeroCoefficentsBool &&  aeroBool == 0
      error("elemTypeAero chord vector must be defined in elements struct \n")
    end
    if aeroBool
      % extract wind function name
      userWindVel = analysisSettings.userWindVel ;
      % extract nonLinearity in aero force boolean
      geometricNonLinearAero = analysisSettings.geometricNonLinearAero ;
      numGaussPoints = 2 ;
      % read aero paramters of the element
      elemTypeAero      = elements( mebiVec( 2 ) ).elemTypeAero ;
      if ~isempty(userDragCoef)
        userDragCoef    = elements( mebiVec( 2 ) ).userDragCoef ;
      end
      if ~isempty(userLiftCoef)
        userLiftCoef    = elements( mebiVec( 2 ) ).userLiftCoef ;
      end
      if ~isempty(userMomentCoef)
        userMomentCoef  = elements( mebiVec( 2 ) ).userMomentCoef ;
      end
      elemTypeGeometry  = elements( mebiVec( 2 ) ).elemTypeGeometry;
      % compute force
      FaeroElem = aeroForce( elemNodesxyzRefCoords, elemTypeGeometry    ,
                             u2ElemDisps( Ut       , dofselem )         ,
                             u2ElemDisps( Udott    , dofselem )         ,
                             u2ElemDisps( Udotdott , dofselem )         ,
                             userDragCoef, userLiftCoef, userMomentCoef ,
                             elemTypeAero, userWindVel ,  geometricNonLinearAero,
                             timeVar ) ;

    end

   % ---------  triangle solid element -----------------------------
  elseif strcmp( elemType, 'triangle')

    thickness = elemTypeGeometry ;

    if strcmp( hyperElasModel, 'linearElastic' )

    % -----------   aerodynmamic force   ------------------------------------
      planeStateFlag = elemTypeParams ;
      dotdotdispsElem  = u2ElemDisps( Udotdott , dofselemRed ) ;

      [ fs, ks, stress ] = elementTriangSolid( elemNodesxyzRefCoords, elemDisps, ...
                            [1 hyperElasParams], 2, thickness, planeStateFlag, dotdotdispsElem, density ) ;
        %
        Finte = fs{1};
        Ke    = ks{1};
        Fmase = fs{3};
        Mmase = ks{3};
        Ce = zeros( size( Mmase ) ) ; % only global damping considered (assembled after elements loop)

    end

  % ---------  tetrahedron solid element -----------------------------
  elseif strcmp( elemType, 'tetrahedron')

    if strcmp( hyperElasModel, 'SVK' )
      auxMatNum = 2 ;

    elseif strcmp( hyperElasModel, 'NHC' )
      auxMatNum = 3 ;
    else
      hyperElasModel
      error('material not implemented yet! open an issue.')
    end

   if isempty(elemTypeParams)
     % (1 analytic 2 complex step)
     consMatFlag = 1 ; % default: 1
   else
     consMatFlag = elemTypeParams(1) ;
   end
   [ Finte, Ke, stress ] = elementTetraSolid( elemNodesxyzRefCoords, elemDisps, ...
                            [ auxMatNum hyperElasParams], 2, consMatFlag ) ;

  end   % case in typee of element ----
  % -------------------------------------------

  %md### Assembly
  %md
  if fsBool
    % internal loads vector assembly
    Fint ( dofselemRed ) = Fint( dofselemRed ) + Finte ;
    if dynamicProblemBool
      Fmas ( dofselemRed ) = Fmas( dofselemRed ) + Fmase ;
    end
    if aeroBool
      Faero( dofselemRed ) = Faero( dofselemRed ) + FaeroElem ;
    end
  end

  if tangBool
    for indRow = 1:length( dofselemRed )

      entriesSparseStorVecs = counterInds + (1:length( dofselemRed) ) ;

      indsIK ( entriesSparseStorVecs )  = dofselemRed( indRow ) ;
      indsJK ( entriesSparseStorVecs )  = dofselemRed ;
      valsK  ( entriesSparseStorVecs )  = Ke( indRow, : )' ;

      if dynamicProblemBool
        valsM( entriesSparseStorVecs ) = Mmase( indRow, : )' ;
        if exist('Ce')~=0
          valsC( entriesSparseStorVecs ) = Ce   ( indRow, : )' ;
        end
      end

      counterInds = counterInds + length( dofselemRed ) ;
    end
  end


  if stressBool
    stressMat( elem, (1:length(stressElem) ) ) = stressElem ;
  end % if stress

end % for elements ----


% ============================================================================
%  --- 3 global additions and output ---
% ============================================================================

fsCell       = cell( 3, 1 ) ;
tangMatsCell = cell( 3, 1 ) ;

if dynamicProblemBool
  dampingMat          = sparse( nNodes*6, nNodes*6 ) ;
  dampingMat(1:2:end) = nodalDispDamping             ;
  dampingMat(2:2:end) = nodalDispDamping * 0.01      ;
end

if fsBool
  Fint = Fint + KS * Ut ;

  fsCell{1} = Fint ;

  if dynamicProblemBool,
    Fvis = dampingMat * Udott ;
  end

  fsCell{2} = Fvis  ;
  fsCell{3} = Fmas  ;
  fsCell{4} = Faero ;
end

if tangBool

  indsIK = indsIK(1:counterInds) ;
  indsJK = indsJK(1:counterInds) ;
  valsK  = valsK (1:counterInds) ;
  K      = sparse( indsIK, indsJK, valsK, size(KS,1), size(KS,1) ) + KS ;

  tangMatsCell{1} = K ;

  if dynamicProblemBool
    valsM = valsM (1:counterInds) ;
    valsC = valsC (1:counterInds) ;
    M     = sparse( indsIK, indsJK, valsM , size(KS,1), size(KS,1) )  ;
    C     = sparse( indsIK, indsJK, valsC , size(KS,1), size(KS,1) ) + dampingMat ;
  else
    M = sparse(size( K ) ) ;
    C = sparse(size( K ) ) ;
  end

  tangMatsCell{2} = C ;
  tangMatsCell{3} = M ;
end

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

elemDisps = U( dofselem ) ;
