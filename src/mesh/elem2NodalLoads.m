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
 
%md function that constructs the assembled Fext vector for one given BC

function fext = elem2NodalLoads ( Conec, indBC, elements, boundaryCond, Nodes )

  % declare output fext
  nnodes      = size( Nodes, 1)      ;
  fext        = zeros( 6*nnodes, 1 ) ;

  % get element indexes with current BC
  elemsWithBC = find( Conec(:,3) == indBC ) ;

  % extract BC load base vals and coord system
  loadCoordSys = boundaryCond.loadsCoordSys ;
  loadvals     = boundaryCond.loadsBaseVals ;

  loadedNodes = []                   ;

  % loop in elements with current BC
  for elemInd = 1:length( elemsWithBC )

    elem      = elemsWithBC( elemInd )       ;
    elemInd   = Conec( elem, 2 )             ;
    elemType  = elements( elemInd ).elemType ;

    %md nodal loads
    if strcmp( elemType, 'node') % node

      if strcmp( loadCoordSys, 'global' )
        nodes     = Conec( elem, 4+1 ) ;
      else
        error(' only global flag in load by now.');
      end
      elemNodeLoadsMatrix = loadvals ;

    %md truss
    elseif strcmp( elemType , 'truss') ; %
      error(' not yet.');

    %md frame
    elseif strcmp( elemType , 'frame') ; %
      error(' not yet.');


    %md edge
    elseif strcmp( elemType , 'edge') ; %
      nodes          = Conec( elem, 4+(1:2) ) ;
      
      % vector from node 1 to node 2
      directionVector = Nodes( nodes(2),:) - Nodes( nodes(1),:) ;
      % length of the edge
      lengthElem = norm(directionVector) ;
      % thickness of the edge element
      thickness  = elements( elemInd ).elemCrossSecParams ;
      
      % check oriented vector is defined into x-y plane
      assert( abs( directionVector(3) ) < eps * 10 ,'edge must be defined into x-y plane' )

      % factor of the load for each node
      factor = lengthElem * thickness * 0.5 ;

      if strcmp( loadCoordSys, 'global' )
        Fx = loadvals( 1 ) * factor ;
        Fy = loadvals( 2 ) * factor ;
        Fz = 0 ;
      elseif strcmp( loadCoordSys, 'local' )
        % consider a 90 degrees rotation of the oriented vector of the line element
        % tanget unitary vector
        tangUniVec = directionVector / lengthElem ;
        % normal unitary vector
        normalUniVec = cross( [ 0 0 1 ] , tangUniVec ) ;
        % tension vector
        tensionVec = ( loadvals( 1 ) * tangUniVec  + loadvals( 2 ) * normalUniVec ) * factor; 
        % nodal forces in global coordinates
        Fx = tensionVec(1) ;
        Fy = tensionVec(2) ;
        Fz = 0 ;
      end % if global/local system

      elemNodeLoadsMatrix = ones( length(nodes), 1 )*[Fx 0 Fy 0 Fz 0] ;

      assert( size( elemNodeLoadsMatrix, 2)==6,'error, maybe missing thickness')

    %md triangle tension
    elseif strcmp( elemType , 'triangle') ; %

      nodes = Conec( elem, 4+(1:3) ) ;

      areaElem = 0.5 * norm( cross( ...
        Nodes( nodes(2),:) - Nodes( nodes(1),:) , ...
        Nodes( nodes(3),:) - Nodes( nodes(1),:) ...
        ) ) ;

      if strcmp( loadCoordSys, 'global' )

        Fx = loadvals( 1 ) * areaElem / 3 ;
        Fy = loadvals( 3 ) * areaElem / 3 ;
        Fz = loadvals( 5 ) * areaElem / 3 ;

        assert( sum( loadvals( [ 2 4 6 ] ) == 0 ) == 3, ...
          'error only pressure loads, not moments, create an issue!' );

      elseif strcmp( loadCoordSys, 'local' ) % local coordinates load

        dofsaux = nodes2dofs( nodes , 6 ) ;
        dofs    = dofsaux(1:2:length(dofsaux)) ;
        % compute the normal vector of the element
        normal  = cross( ...
          Nodes( nodes(2),:) - Nodes( nodes(1),:) , ...
          Nodes( nodes(3),:) - Nodes( nodes(1),: ) ) ;
        % and normalize it
        n = normal / norm( normal ) ;

        Fx = n(1) * loadvals( 5 ) * areaElem / 3 ;
        Fy = n(2) * loadvals( 5 ) * areaElem / 3 ;
        Fz = n(3) * loadvals( 5 ) * areaElem / 3 ;

        assert( sum( loadvals( [ 1 2 3 4 6 ] ) == 0 ) == 5, ...
          'error only normal pressure loads in local coords, create an issue!' );

      else
        loadCoordSys
        error(' loadsCoordSys field must be local or global.');
      end % if global/local system

      elemNodeLoadsMatrix = ones( length(nodes), 1 ) * [ Fx 0 Fy 0 Fz 0 ] ;

    end %if elemTypes
    %mdadd loads to matrix of loaded nodes
    loadedNodes = [ loadedNodes ; ...
                    nodes'  elemNodeLoadsMatrix ] ;
  end % for elements

  %md convert to assembled fext vetor
  if exist( 'loadedNodes' ) ~= 0
    for i=1:size( loadedNodes ,1)
      aux = nodes2dofs ( loadedNodes(i,1), 6 ) ;
      fext( aux ) = fext( aux ) + loadedNodes(i,2:7)' ;
    end
  end
