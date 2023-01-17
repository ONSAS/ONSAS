% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
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
function [ fsCell, stressMat, tangMatsCell, matFint, strain_vec, acum_plas_strain_vec ] = assembler( Conec, elements, Nodes,...
                                                           materials, KS, Ut, Udott, Udotdott,...
                                                           analysisSettings, outputBooleans, nodalDispDamping,...
                                                           timeVar, previous_state_mat )

fsBool     = outputBooleans(1) ; stressBool = outputBooleans(2) ; tangBool   = outputBooleans(3) ; matFintBool = outputBooleans(4) ;

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

  % "allocates" space for the bigest possible matrices (4 nodes per element)
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

% -------  matrix with internal forces per element -------------------
if matFintBool
	matFint = zeros( nElems, 6*4 ) ;
else
	matFint = [] ;
end

stress_n_vec           =  previous_state_mat(:,1) ;
strain_n_vec           = previous_state_mat(:,2)  ;
acum_plas_strain_n_vec =  previous_state_mat(:,3) ;

strain_vec = zeros(size( strain_n_vec )) ;
acum_plas_strain_vec = zeros(size(acum_plas_strain_n_vec)) ;


% ====================================================================


dynamicProblemBool = strcmp( analysisSettings.methodName, 'newmark' ) || strcmp( analysisSettings.methodName, 'alphaHHT' ) ;

% ====================================================================
%  --- 2 loop assembly ---
% ====================================================================

for elem = 1:nElems

  mebiVec = Conec( elem, 1:4) ;

  %md extract element properties
  hyperElasModel     = materials( mebiVec( 1 ) ).hyperElasModel   ;
  hyperElasParams    = materials( mebiVec( 1 ) ).hyperElasParams  ;
  density            = materials( mebiVec( 1 ) ).density          ;

  elemType           = elements( mebiVec( 2 ) ).elemType          ;
  elemTypeParams     = elements( mebiVec( 2 ) ).elemTypeParams    ;
  massMatType        = elements( mebiVec( 2 ) ).massMatType       ;
  elemCrossSecParams = elements( mebiVec( 2 ) ).elemCrossSecParams;

  %md extract aerodynamic properties
  elemTypeAero       = elements( mebiVec( 2 ) ).elemTypeAero      ;
  aeroCoefs          = elements( mebiVec( 2 ) ).aeroCoefs         ;

  %md compute aerodynamic compute force booleans
  aeroBool = ~isempty(analysisSettings.fluidProps) || ...
             ~isempty( elemTypeAero ) || ~isempty( aeroCoefs ) ;

  %md obtain element info
  [numNodes, dofsStep] = elementTypeInfo ( elemType ) ;

  %md obtains nodes and dofs of element
  nodeselem   = Conec( elem, (4+1):(4+numNodes) )' ;
  dofselem    = nodes2dofs( nodeselem , 6 )        ;
  dofselemRed = dofselem ( 1 : dofsStep : end )    ;

  %md elemDisps contains the displacements corresponding to the dofs of the element
  elemDisps   = u2ElemDisps( Ut , dofselemRed ) ;

  %md dotdotdispsElem contains the accelerations corresponding to the dofs of the element
  if dynamicProblemBool
    dotdotdispsElem  = u2ElemDisps( Udotdott , dofselemRed ) ;
  end

  elemNodesxyzRefCoords  = reshape( Nodes( Conec( elem, (4+1):(4+numNodes) )' , : )',1,3*numNodes) ;

  stressElem = [] ;

  % -----------   node element   ------------------------------
  if strcmp( elemType, 'node')
    nodalMass = materials( mebiVec( 1 ) ).nodalMass ;
    (iscolumn(nodalMass)) && ( nodalMass == nodalMass' ) ;

    Finte = zeros(3,1) ;    Ke    = zeros(3,3) ;

    if dynamicProblemBool
      % Mmase = spdiags( nodalMass', 0, 3, 3 ) ; % sparse option for someday...
      Mmase = diag( nodalMass ) ;
      Fmase = Mmase * dotdotdispsElem        ;
    end

  % -----------   truss element   ------------------------------
  elseif strcmp( elemType, 'truss')

    A  = crossSectionProps ( elemCrossSecParams, density ) ;
    previous_state = [ stress_n_vec(elem) strain_n_vec(elem) acum_plas_strain_n_vec(elem) ] ;

    [ fs, ks, stressElem, ~, strain, acum_plas_strain ] = elementTrussInternForce( elemNodesxyzRefCoords, elemDisps, hyperElasModel, hyperElasParams, A, previous_state ) ;

    Finte = fs{1} ;  Ke = ks{1} ;

    if dynamicProblemBool
      dotdotdispsElem  = u2ElemDisps( Udotdott , dofselem ) ;
      [ Fmase, Mmase ] = elementTrussMassForce( elemNodesxyzRefCoords, density, A, massMatType, dotdotdispsElem ) ;
      %
      Ce = zeros( size( Mmase ) ) ; % only global damping considered (assembled after elements loop)
    end


  % -----------   frame element   ------------------------------------
  elseif strcmp( elemType, 'frame')

		if strcmp(hyperElasModel, 'linearElastic')
			boolLinear = 1 ;
			boolMatNonLin = 0 ;
		elseif strcmp(hyperElasModel, 'biLinear') || strcmp(hyperElasModel, 'userFunc')
			boolLinear = 1 ;
			boolMatNonLin = 1 ;
		else
			boolLinear = 0 ;
		end

		if  boolLinear == 1

			[ fs, ks ] = linearStiffMatBeam3D(elemNodesxyzRefCoords, elemCrossSecParams, massMatType, density, hyperElasModel, hyperElasParams, u2ElemDisps( Ut, dofselem ), u2ElemDisps( Udotdott , dofselem ), tangBool, boolMatNonLin, matFintBool, elem ) ;

      Finte = fs{1} ;  Ke = ks{1} ;

      if dynamicProblemBool
        Fmase = fs{3} ; Mmase = ks{3} ;
      end

		elseif strcmp( hyperElasModel, '1DrotEngStrain')

      [ fs, ks, stress, rotData ] = frame_internal_force( elemNodesxyzRefCoords , ...
                                                             elemCrossSecParams    , ...
                                                             [ 1 hyperElasParams ] , ...
                                                             u2ElemDisps( Ut, dofselem ) ) ;
      Finte = fs{1} ;  Ke = ks{1} ;

      if dynamicProblemBool

        [ fs, ks  ] = frame_inertial_force( elemNodesxyzRefCoords               , ...
                                            elemCrossSecParams                  , ...
                                            [ 1 hyperElasParams ]               , ...
                                            u2ElemDisps( Ut, dofselem )         , ...
                                            u2ElemDisps( Udott    , dofselem )  , ...
                                            u2ElemDisps( Udotdott , dofselem )  , ...
                                            density, massMatType  ) ;


        Fmase = fs{3} ; Ce = ks{2} ; Mmase = ks{3} ;
      end
    else
      error('wrong hyperElasModel for frame element.')
    end

    %md compute fluid forces on the element
    if aeroBool && fsBool
      % extract boolean to compute aerodynamic displacements tangent matrix
      aeroTangBool = elements( mebiVec( 2 ) ).elemTypeAero(5) ;

      [FaeroElem, MataeroEelem] = frame_fluid_force( elemNodesxyzRefCoords, ...
                                     elemCrossSecParams                   , ...
                                     u2ElemDisps( Ut       , dofselem )   , ...
                                     u2ElemDisps( Udott    , dofselem )   , ...
                                     u2ElemDisps( Udotdott , dofselem )   , ...
                                     elements( mebiVec( 2 ) ).aeroCoefs, elements( mebiVec( 2 ) ).elemTypeAero,...
                                     analysisSettings, timeVar, elem, aeroTangBool ) ;
    end

  % ---------  triangle solid element -----------------------------
  elseif strcmp( elemType, 'triangle')

    thickness = elemCrossSecParams ;

    if strcmp( hyperElasModel, 'linearElastic' )

      planeStateFlag = elemTypeParams ;
      dotdotdispsElem  = u2ElemDisps( Udotdott , dofselemRed ) ;

      [ fs, ks, stress ] = elementTriangSolid( elemNodesxyzRefCoords, elemDisps, ...
                            [1 hyperElasParams], 2, thickness, planeStateFlag, dotdotdispsElem, density ) ;

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
   [ Finte, Ke, stressElem ] = elementTetraSolid( elemNodesxyzRefCoords, elemDisps, ...
                            [ auxMatNum hyperElasParams], 2, consMatFlag ) ;

  end   % case in typee of element ----
  % -------------------------------------------


  %md### Assembly
  %md
  if fsBool
    % internal loads vector assembly
    if norm( Finte ) > 0.0
      Fint ( dofselemRed ) = Fint( dofselemRed ) + Finte ;
    end
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

      if aeroBool && exist('aeroTangBool') && aeroTangBool
        % add displacements minus since is an external force
        valsK  ( entriesSparseStorVecs )  = Ke( indRow, : )' - MataeroEelem( indRow, : )' ;
      else
        valsK  ( entriesSparseStorVecs )  = Ke( indRow, : )' ;
      end

      if dynamicProblemBool
        valsM( entriesSparseStorVecs ) = Mmase( indRow, : )' ;
        if exist('Ce')~=0
          valsC( entriesSparseStorVecs ) = Ce( indRow, : )' ;
        end
      end

      counterInds = counterInds + length( dofselemRed ) ;
    end
  end


  if stressBool
    stressMat( elem, (1:length(stressElem) ) ) = stressElem ;

    if exist('strain')==1
      strain_vec( elem )           = strain ;
      acum_plas_strain_vec( elem,1 ) = acum_plas_strain ;
    end
  end % if stress

	if matFintBool
		matFint(elem,1:dofsStep:numNodes*6) = Finte' ;
	end

end % for elements ----

% ============================================================================
%  --- 3 global additions and output ---
% ============================================================================

fsCell       = cell( 3, 1 ) ;
tangMatsCell = cell( 3, 1 ) ;

if dynamicProblemBool
  dampingMat          = sparse( nNodes*6, nNodes*6 ) ;
  dampingMat(1:2:end,1:2:end) = nodalDispDamping*speye(nNodes*3,nNodes*3)             ;
  % dampingMat(2:2:end,1:2:end) = nodalDispDamping * 0.00      ;
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

  global globalReactionForces
  global glboalNodeReactionForces
  if ~isempty(globalReactionForces) && (round(timeVar) == timeVar) && (timeVar ~= 0)
    dofsRForces = (glboalNodeReactionForces - 1) * 6 + 1 : glboalNodeReactionForces * 6  ;
    globalReactionForces((timeVar -1)*6 + 1: (timeVar)*6) = Faero(dofsRForces) - Fint(dofsRForces) - Fmas(dofsRForces) - Fvis(dofsRForces) ;
  end

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
