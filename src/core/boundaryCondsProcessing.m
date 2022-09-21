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
 
%md This function converts the mesh MEBI information to the data structures used in the numerical simulation

function [ Conec, Nodes, factorLoadsFextCell, loadFactorsFuncCell, diriDofs, neumDofs, KS, userLoadsFilename, userWindVel ] = boundaryCondsProcessing ( mesh, ...
                        materials, ...     % M
                        elements, ...      % E
                        boundaryConds, ... % B
                        analysisSettings  )

Conec  = myCell2Mat( mesh.conecCell ) ;
Nodes  = mesh.nodesCoords ;
nnodes = size( Nodes,1);
KS     = sparse( 6*nnodes, 6*nnodes );


%md Since we want to process the BCs, we keep only the nonzero BCs
%md Computes the number of elements and BCs we have
boundaryTypes  = unique( Conec( :, 3) ) ;
if boundaryTypes(1) == 0,
  boundaryTypes(1)=[];
end

elementTypes   = unique( Conec( :, 2) ) ;
if elementTypes(1)  == 0, % checks if all elements have a type
  error('all elements must be defined');
end

factorLoadsFextCell = {} ;
loadFactorsFuncCell = {} ;
diriDofs            = [] ;

%md loop over the boundary conditions used in the mesh
for indBC = 1:length( boundaryTypes )

  % number of current BC processed
  BCnum = boundaryTypes(indBC) ;
  
  %md loads verification
  %md is loadsCoordSys is not empty, then some load is applied in this BC
  if ~isempty( boundaryConds( indBC ).loadsCoordSys )
    %md The nodal loads vector is computed and assiged to the corresponding BC entry.
    factorLoadsFextCell{ BCnum } = elem2NodalLoads ( Conec, BCnum, elements, boundaryConds( BCnum ),  Nodes ) ;
    % defaul load factor function
    if isempty( boundaryConds(BCnum).loadsTimeFact ),
      boundaryConds(BCnum).loadsTimeFact = @(t) t ;
    end
    loadFactorsFuncCell{ BCnum } = boundaryConds(BCnum).loadsTimeFact ;
  end % if load

  %md displacement verification
  if ~isempty( boundaryConds(BCnum).imposDispDofs ),
    
    %md values and imposed dofs of current BC
    impoDofs = boundaryConds(BCnum).imposDispDofs ;
    impoVals = boundaryConds(BCnum).imposDispVals ;
    
    %md find the elements with the current boundary condition
    elemsWithBC = find( Conec(:,3) == indBC ) ;
    
    %md compute the values of nonHom BC
    [ nonHomDiriVals, bcDiriDofs, nonHomDiriDofs ] = elem2NodalDisps ( Conec, indBC, elemsWithBC, elements, impoDofs, impoVals, Nodes ) ;
    diriDofs = [ diriDofs; bcDiriDofs ] ;

  end % if: disp dofs



  if ~isempty( boundaryConds(BCnum).springDofs ),
    
    %md find the elements with the current boundary condition
    elemsWithBC = find( Conec(:,3) == indBC ) ;

    for inde = 1:length( elemsWithBC )
      %
      if ~strcmp( elements( Conec( elemsWithBC(inde), 2 ) ).elemType, 'node' )    
        error('springDofs can only be assigned to node elements, by the moment!')
      else
        nodes_with_spring = Conec( elemsWithBC(inde), 4+1 )  ;
        for indn = 1:length( nodes_with_spring )
          nodedofs = nodes2dofs( nodes_with_spring(indn), 6 ) ;
          locdofs  = boundaryConds(BCnum).springDofs ;
          KS( nodedofs(locdofs), nodedofs(locdofs) ) = KS( nodedofs(locdofs), nodedofs(locdofs) ) ...
            + spdiags( boundaryConds(BCnum).springVals, 0, length(locdofs), length(locdofs) ) ;
        end
      end
    end

  end % if: springDofs


end % for: elements with boundary condition assigned

%full(KS)

diriDofs = unique( diriDofs) ;

%md remove element if no material is assigned
elemsToRemove = find( Conec( :, 1 ) == 0 ) ;
Conec( elemsToRemove, :  ) = [] ;

%md construction of a vector with the neumann degrees of freedom

%md a zeros filled vector is created
neumDofs = zeros( 6*nnodes, 1 ) ; % maximum possible vector

%md loop for construction of vector of dofs
for elemNum = 1:length( elementTypes )

  %md find the numbers of the elements with the current element type
  elementsNums = find( Conec( :, 2 ) == elementTypes( elemNum ) ) ;

  %md get current element type
  elemType = elements( elementTypes(elemNum) ).elemType ;

  %md if there are any elements with this type
  if length( elementsNums ) > 0

    [numNodes, dofsStep] = elementTypeInfo ( elemType ) ;
    nodes    = Conec( elementsNums, (4+1):(4+numNodes) ) ;
    dofs     = nodes2dofs( nodes, 6)'       ;
    dofs     = dofs(1:dofsStep:end)         ;

    %md remove z displacement for triangles with material assigned (plane movement in x-y)
    if strcmp( elemType, 'triangle' );
      dofs(3:3:end) = [ ] ;
    end

    neumDofs ( dofs ) = dofs ;
  end
end

if length( diriDofs)>0
  neumDofs( diriDofs ) = 0 ;
end
neumDofs = unique( neumDofs ) ;
if neumDofs(1) == 0,
  neumDofs(1)=[] ;
end

% read user load filename
userLoadsFilename = [] ;

if isfield( boundaryConds,'userLoadsFilename' )
  for ind = 1:length(boundaryConds)
    if ~isempty( boundaryConds(ind).userLoadsFilename )
      if ~isempty( userLoadsFilename )
        error('only one user load function can be used!');
      else
        userLoadsFilename = boundaryConds(ind).userLoadsFilename ;
      end
    end
  end
end

% FIx me when userWindVel boundaryConds(wind is not in the first cell)
if isfield( boundaryConds,'userWindVel' )
  userWindVel  = boundaryConds.userWindVel     ;
else
  userWindVel     = [];
end
% ----------------------------------------------------------------------


% ----------------------------------------------------------------------
%md Loop for computing of gravity external force vector
if analysisSettings.booleanSelfWeight == true
  % initialize the gravity load vector
  gravityFactorLoads = zeros( 6*nnodes, 1 ) ; % maximum possible vector
  g = 9.8067 ;

  %md loop in type of element indexes
  for elemNum = 1:length( elementTypes )
    %md find the numbers of the elements with the current element type
    elementsNums = find( Conec( :, 2 ) == elementTypes( elemNum ) ) ;

    %md if there is any element
    if length( elementsNums ) > 0

      %md get current element type
      elemType = elements( elementTypes(elemNum) ).elemType ;

      %md get the material types of the current element type
      materialElemTypes   = unique( Conec( elementsNums, 1) ) ;

      % keep positive-rho materials
      aux = [ ] ;
      for ir = 1:length( materialElemTypes )
        if materials( materialElemTypes(ir) ).density > 0
          aux = [ aux; materialElemTypes(ir) ] ;
        end
      end
      materialElemTypes = aux ;

      %md  loop in different materials of current element
      for matElem = 1: length(materialElemTypes)

        %md  extract the material of the current element type
        rhoElemTypeMaterial = materials( materialElemTypes(matElem) ).density ;

        %md  find the elements with elements = elemNum and material = matElem
        elemntsSameMat  = elementsNums(find(Conec( elementsNums, 1 ) == materialElemTypes(matElem) ) ) ;

        if strcmp( elemType, 'truss') || strcmp( elemType, 'frame') ;

          %md extract connectivity of element
          conecElemsThisMat = Conec( elemntsSameMat, 5:6 );
          % final minus initial nodes coords matrices
          matrixDiffs = Nodes( conecElemsThisMat(:,2),:) - Nodes( conecElemsThisMat(:,1), : );
          vecSquareDiff = sum( (matrixDiffs.^2)' )';
          lengthElems = sqrt( vecSquareDiff ) ;
          elemCrossSecParams = elements(elemNum).elemCrossSecParams ;
          [areaElem, ~, ~, ~, ~ ] = crossSectionProps ( elemCrossSecParams, rhoElemTypeMaterial ) ;

          %md compute nodal selfweight loads
          Fz = - rhoElemTypeMaterial * lengthElems * areaElem * g * 0.5;

          %md compute the dofs where gravitiy is applied: and add loads
          for elem = 1: size(conecElemsThisMat, 1)
            dofs  = nodes2dofs( conecElemsThisMat(elem,:), 6)';
            dofs  = dofs(5:6:end);
            gravityFactorLoads(dofs,1) = gravityFactorLoads(dofs,1) + Fz(elem);
          end

        elseif strcmp( elemType, 'tetrahedron')

          %md extract connectivity of element
          conecElemsThisMat = Conec( elemntsSameMat, 5:8 );

          %md compute the dofs where gravitiy is applied: and add loads
          for elem = 1: size(conecElemsThisMat, 1)
            auxnodes = conecElemsThisMat(elem,:) ;
            dofs     = nodes2dofs( auxnodes, 6)' ;
            dofs     = dofs(5:6:end) ;

            elemCoords = reshape( Nodes( auxnodes, : )', 1, 12 ) ;

            [ ~, ~, vol ] = computeFuncDerivVolTetraSolid( elemCoords ) ;

            Fzval = - vol *.25 * rhoElemTypeMaterial * g ;

            gravityFactorLoads(dofs,1) = gravityFactorLoads(dofs,1) + Fzval ;
          end

        else
          error('this elemType is not implemented using self weight boolean yet')

        end % element type case

      end % for materials

    end % if there are elements

  end % for elementTypes

  %md the number of BC that represent the self weight condition is
  numberOfBCSelfWeight = length( factorLoadsFextCell )  + 1         ;
  factorLoadsFextCell{numberOfBCSelfWeight} = gravityFactorLoads    ;
  loadFactorsFuncCell{numberOfBCSelfWeight} = @(t) 1                ;

end % if booleanSelfWeight
