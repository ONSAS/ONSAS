% Copyright 2024, ONSAS Authors (see documentation)
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
%mdThis function computes the assembled force vectors, tangent matrices and stress matrices.
function [ fsCell, stressMat, tangMatsCell, ... 
            localInternalForces, strain_vec, acum_plas_strain_vec ] = assembler( ...
             Conec, elements, Nodes,...
             materials, KS, Ut, Udott, Udotdott,...
             analysisSettings, outputBooleans, nodalDispDamping,...
             timeVar, previousStateCell )

% ====================================================================
%  --- 1 declarations ---
% ====================================================================

fsBool   = outputBooleans(1) ; stressBool = outputBooleans(2) ;
tangBool = outputBooleans(3) ; matFintBool = outputBooleans(4) ;

nElems   = size(Conec, 1) ;
nNodes   = size(Nodes, 1) ;

% -------  forces vectors -------------------------------------------
if fsBool
  % --- creates Fint vector ---
  Fint  = zeros( nNodes*6 , 1 ) ;
  Fmas  = zeros( nNodes*6 , 1 ) ;
  Fvis  = zeros( nNodes*6 , 1 ) ;
  Faero = zeros( nNodes*6 , 1 ) ;
  Fther = zeros( nNodes*6 , 1 ) ;
end

% -------  tangent matrices   -------------------------------------
if tangBool

  % "allocates" space for the bigest possible matrices (4 nodes per element)
  indsIK = zeros( nElems*24*24, 1 )   ;
  indsJK = zeros( nElems*24*24, 1 )   ;
  valsK  = zeros( nElems*24*24, 1 )   ;

  valsC  = zeros( nElems*24*24, 1 )   ;
  valsM  = zeros( nElems*24*24, 1 )   ;

  counterInds = 0 ; % counter non-zero indexes
end

% -------  matrix with stress per element ----------------------------
% due to octave constraints a matrix of zeros data is used,
% however a boolean first column was added to save 1 if the element stress/fint was computed of zero if not.
%
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

localInternalForces = struct();

% Previous state
stress_n_vec           =  previousStateCell(:,1) ;
strain_n_vec           =  previousStateCell(:,2) ;
acum_plas_strain_n_vec =  previousStateCell(:,3) ;

strain_vec = cell( size(strain_n_vec, 1), 1 ) ;
acum_plas_strain_vec = cell( size(acum_plas_strain_n_vec, 1), 1 ) ;

dynamicProblemBool = strcmp( analysisSettings.methodName, 'newmark' ) ...
                  || strcmp( analysisSettings.methodName, 'alphaHHT' ) ;
% ====================================================================


% ====================================================================
%  --- 2 loop assembly ---
% ====================================================================

for elem = 1:nElems

  mebVec = Conec( elem, 1:3) ;

  %md extract element properties
  modelName          = materials( mebVec( 1 ) ).modelName   ;
  modelParams        = materials( mebVec( 1 ) ).modelParams  ;
  density            = materials( mebVec( 1 ) ).density          ;

  elemType           = elements( mebVec( 2 ) ).elemType          ;
  elemTypeParams     = elements( mebVec( 2 ) ).elemTypeParams    ;
  massMatType        = elements( mebVec( 2 ) ).massMatType       ;
  elemCrossSecParams = elements( mebVec( 2 ) ).elemCrossSecParams;

  %md extract aerodynamic properties
  aeroCoefs           = elements( mebVec( 2 ) ).aeroCoefFunctions       ;
  chordVector         = elements( mebVec( 2 ) ).chordVector         ;
  aeroNumericalParams = elements( mebVec( 2 ) ).aeroNumericalParams ;

  %md compute aerodynamic force booleans
  aeroBool = ~isempty(analysisSettings.fluidProps) ;

  %md obtain element info
  [numNodes, nodalDofsEntries] = elementTypeDofs( elemType ) ;

  % obtains nodes and dofs of element
  nodeselem   = Conec( elem, (3+1):(3+numNodes) )' ;
  dofselem    = nodes2dofs( nodeselem , 6 )   ;     

  % construct vector of degrees of freedom of element
  auxA = repmat( nodalDofsEntries, length(dofselem)/6,1 )  ;
  auxB = repelem( (0:6:length(dofselem)-1)',length(nodalDofsEntries),1) ;
  dofselemRed = dofselem( auxA+auxB )   ;


  %md elemDisps contains the displacements corresponding to the dofs of the element
  elemDisps       = u2ElemDisps( Ut      , dofselemRed ) ;
  dotdispsElem    = u2ElemDisps( Udott   , dofselemRed ) ;
  dotdotdispsElem = u2ElemDisps( Udotdott, dofselemRed ) ;

  elemNodesxyzRefCoords  = reshape( Nodes( nodeselem, : )', 1, 3*numNodes ) ;

  stressElem = [] ;
  fintLocCoord = [] ;


  % -----------   node element   ------------------------------
  if strcmp( elemType, 'node')
    nodalMass = materials( mebVec( 1 ) ).nodalMass ;
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
    previous_state = { stress_n_vec{elem}; strain_n_vec{elem}; acum_plas_strain_n_vec{elem} } ;

    [ fs, ks, stressElem, ~, strain, acum_plas_strain ] = elementTrussInternForce( elemNodesxyzRefCoords, elemDisps, modelName, modelParams, A, previous_state ) ;

    Finte = fs{1} ;  Ke = ks{1} ;

    localInternalForces(elem).Nx = norm( Finte(1:3) ) ;

    if dynamicProblemBool
      [ Fmase, Mmase ] = elementTrussMassForce( elemNodesxyzRefCoords, density, A, massMatType, dotdotdispsElem ) ;
      %
      Ce = zeros( size( Mmase ) ) ; % only global damping considered (assembled after elements loop)
    end
  
  % -----------   frame element   ------------------------------------
  elseif strcmp( elemType, 'frame')
    
		if strcmp(modelName, 'elastic-linear')

			[ fs, ks, fintLocCoord ] = elementFrameLinear(elemNodesxyzRefCoords, elemCrossSecParams, massMatType, density, modelName, modelParams, elemDisps, dotdotdispsElem) ;
      
      Finte = fs{1} ;  Ke = ks{1} ;

      Nx = fintLocCoord(1);   My = fintLocCoord(4);   Mz = fintLocCoord(6);
      if dynamicProblemBool
        Fmase = fs{3} ; Mmase = ks{3} ;
      end

		elseif strcmp( modelName, 'elastic-rotEngStr')

      [ fs, ks, stress, rotData, fintLocCoord ] = frame_internal_force( elemNodesxyzRefCoords , ...
                                                             elemCrossSecParams    , ...
                                                             [ 1 modelParams ] , ...
                                                             elemDisps ) ;
      Finte = fs{1} ;  Ke = ks{1} ;

      Nx = fintLocCoord(1);   My = fintLocCoord(2);   Mz = fintLocCoord(3);

      if dynamicProblemBool
        [ fs, ks  ] = frame_inertial_force( elemNodesxyzRefCoords , elemCrossSecParams, ...
                                            [ 1 modelParams ], elemDisps, ...
                                            dotdispsElem, dotdotdispsElem  , ...
                                            density, massMatType, analysisSettings ) ;


        Fmase = fs{3} ; Ce = ks{2} ; Mmase = ks{3} ;
      end
    else
      error('wrong modelName for frame element.')
    end

    localInternalForces(elem).Nx = Nx ;
    localInternalForces(elem).My = My ;
    localInternalForces(elem).Mz = Mz ;

    %md compute fluid forces on the element
    if aeroBool && fsBool
      [FaeroElem, MataeroEelem] = frame_fluid_force( elemNodesxyzRefCoords,        ...
                                     elemCrossSecParams                   ,        ...
                                     elemDisps   ,        ...
                                     dotdispsElem   ,        ...
                                     dotdotdispsElem   ,        ...
                                     aeroCoefs, chordVector, aeroNumericalParams,  ...
                                     analysisSettings, timeVar, elem, ...
                                     aeroNumericalParams{2}  ) ;
                                     

    end

  % ---------  triangle solid element -----------------------------
  elseif strcmp( elemType, 'triangle')

    thickness = elemCrossSecParams ;
		planeStateFlag = elemTypeParams ;
		
		dotdotdispsElem  = u2ElemDisps( Udotdott , dofselemRed ) ;
		
		previous_state = { stress_n_vec{elem} ; strain_n_vec{elem} ; acum_plas_strain_n_vec{elem} } ;
		  
		[ fs, ks, stressElem, strain, acum_plas_strain ] = 	elementTriangSolid( elemNodesxyzRefCoords, elemDisps, ...
																										modelName, [1 modelParams], 2, thickness, planeStateFlag, ...
																										dotdotdispsElem, density, previous_state ) ;

    localInternalForces(elem).Mx  = 0 ;
    localInternalForces(elem).My  = 0 ;
    localInternalForces(elem).Mxy = 0 ;

    Finte = fs{1};
		Ke    = ks{1};
		
		if dynamicProblemBool
			Fmase = fs{3};
			Mmase = ks{3};
			Ce = zeros( size( Mmase ) ) ; % only global damping considered (assembled after elements loop)
		end

  elseif strcmp( elemType, 'triangle-plate')

    thickness = elemCrossSecParams{2};
    
    [ fs, ks, fintLocCoord ] = 	internal_forces_plate_triangle( elemNodesxyzRefCoords, elemDisps, modelName, ...
      modelParams, thickness ) ;

    localInternalForces(elem).Mx  = fintLocCoord(1) ;
    localInternalForces(elem).My  = fintLocCoord(2) ;
    localInternalForces(elem).Mxy = fintLocCoord(3) ;

    Finte = fs{1};
		Ke    = ks{1};

  % ---------  tetrahedron solid element -----------------------------
  elseif strcmp( elemType, 'tetrahedron')

    if strcmp( modelName, 'SVK' )
      auxMatNum = 2 ;

    elseif strcmp( modelName, 'NHC' )
      auxMatNum = 3 ;
    else
      modelName
      error('material not implemented yet! open an issue.')
    end

    localInternalForces(elem).Mx  = 0 ;
    localInternalForces(elem).My  = 0 ;
    localInternalForces(elem).Mxy = 0 ;

   if isempty(elemTypeParams)
     % (1 analytic 2 complex step)
     consMatFlag = 1 ; % default: 1
   else
     consMatFlag = elemTypeParams(1) ;
   end
   [ Finte, Ke, stressElem ] = elementTetraSolid( elemNodesxyzRefCoords, elemDisps, ...
                            [ auxMatNum modelParams], 2, consMatFlag ) ;

  end   % case in type of element ----
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
    if aeroBool && strcmp(elemType,'frame')
      Faero( dofselemRed ) = Faero( dofselemRed ) + FaeroElem ;
    end

    if exist('Fthere')==1 && ( norm( Fthere ) > 0.0 )
      Fther( dofselemRed ) = Fther( dofselemRed ) + Fthere ;
    end
  end

  if tangBool
    for indRow = 1:length( dofselemRed )

      entriesSparseStorVecs = counterInds + (1:length( dofselemRed) ) ;

      indsIK ( entriesSparseStorVecs )  = dofselemRed( indRow ) ;
      indsJK ( entriesSparseStorVecs )  = dofselemRed ;

      if aeroBool && strcmp(elemType,'frame') && aeroNumericalParams{2}
        % add displacements minus since is an external force
        valsK  ( entriesSparseStorVecs )  = Ke( indRow, : )' - MataeroEelem( indRow, : )' ;
      else
        valsK  ( entriesSparseStorVecs )  = Ke( indRow, : )' ;
      end

      if dynamicProblemBool
        valsM( entriesSparseStorVecs ) = Mmase( indRow, : )' ;
        if exist('Ce')==1
          valsC( entriesSparseStorVecs ) = Ce( indRow, : )' ;
        end
      end

      counterInds = counterInds + length( dofselemRed ) ;
    end
  end


  if stressBool
    stressMat( elem, (1:length(stressElem) ) ) = stressElem ;

    if exist('strain')==1
      strain_vec{ elem }           = strain' ;
      acum_plas_strain_vec{ elem } = acum_plas_strain ;
    end
  end % if stress

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
    globalReactionForces((timeVar -1)*6 + 1: (timeVar)*6) = Faero(dofsRForces) - Fint(dofsRForces) - Fmas(dofsRForces) - Fvis(dofsRForces) + Fther(dofsRForces) ;
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

% function nodesmat = conv ( conec, coordsElemsMat )
% nodesmat  = [] ;
% nodesread = [] ;
% for i=1:size(conec,1)
%   for j=1:2
%     if length( find( nodesread == conec(i,j) ) ) == 0
%       nodesmat( conec(i,j),:) = coordsElemsMat( i, (j-1)*6+(1:2:5) ) ;
%     end
%   end
% end

% ==============================================================================
%
% function to convert vector of displacements into displacements of element.
%
% ==============================================================================
% _____&&&&&&&&&&&&&& GENERALIZAR PARA RELEASES &&&&&&&&&&&&&&&
  function elemDisps = u2ElemDisps( U, dofselem)
    elemDisps = U( dofselem ) ;
  end

end