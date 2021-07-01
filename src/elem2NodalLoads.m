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

%md function that constructs the assembled Fext vector for one given BCtype

function fext = elem2NodalLoads ( Conec, indBC, elements, boundaryConds, Nodes )

  elemsWithBC = find( Conec(:,3) == indBC ) ;

  loadedNodes = []                   ;
  nnodes      = size( Nodes, 1)      ;
  fext        = zeros( 6*nnodes, 1 ) ;

  for elemInd = 1:length( elemsWithBC )

    elem            = elemsWithBC( elemInd )               ;
    elemType        = elements.elemType{ Conec( elem, 2 )} ;
    loadCoordSys    = boundaryConds.loadsCoordSys{ indBC } ;

    %md nodal loads
    if strcmp( elemType, 'node') % node

      if strcmp( loadCoordSys, 'global' )
        loadvals = boundaryConds.loadsBaseVals{ indBC } ;
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
      nodes = Conec( elem, 4+(1:2) ) ;

      lengthElem = norm( Nodes( nodes(2),:) - Nodes( nodes(1),:) ) ;
      thickness  = elements.elemTypeGeometry{ Conec( elem, 2 ) } ;
      loadvals   = boundaryConds.loadsBaseVals{ indBC } ;

      factor = lengthElem * thickness * 0.5 ;

      assert( sum( loadvals( [ 2 4 5 6 ] )==0)==4,'error loads added to edge' )

      if strcmp( loadCoordSys, 'global' )
        Fx =   loadvals( 1 ) * factor ;
        Fy =   loadvals( 3 ) * factor ;
        Fz = 0 ;
      elseif strcmp( loadCoordSys, 'local' )
        Fx = - loadvals( 3 ) * factor ;
        Fy =   loadvals( 1 ) * factor ;
        Fz = 0 ;
      end % if global/local system

      elemNodeLoadsMatrix = ones( length(nodes), 1 )*[Fx 0 Fy 0 Fz 0] ;

    %md triangle tension
    elseif strcmp( elemType , 'triangle') ; %

      nodes = Conec( elem, 4+(1:3) ) ;

      areaElem = 0.5 * norm( cross( ...
        Nodes( nodes(2),:) - Nodes( nodes(1),:) , ...
        Nodes( nodes(3),:) - Nodes( nodes(1),:) ...
        ) ) ;

      loadvals = boundaryConds.loadsBaseVals{ indBC } ;

      if strcmp( loadCoordSys, 'global' )

        Fx = loadvals( 1 ) * areaElem / 3 ;
        Fy = loadvals( 3 ) * areaElem / 3 ;
        Fz = loadvals( 5 ) * areaElem / 3 ;

        %~ elseif loadvals(1) == 0 % local coordinates load

        %~ dofsaux = nodes2dofs( nodestrng , 6 ) ;
        %~ dofs    = dofsaux(1:2:length(dofsaux)) ;
        %~ nmod    = norm( cross( ...
        %~ Nodes( nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
        %~ Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) ) ;

        %~ n = cross( ...
        %~ Nodes(nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
        %~ Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) / nmod ;

        %~ Fx = n(1) * loadvals(loadNum, 2+5 ) * area/3 ;
        %~ Fy = n(2) * loadvals(loadNum, 2+5 ) * area/3 ;
        %~ Fz = n(3) * loadvals(loadNum, 2+5 ) * area/3 ;

        %~ if loadvals(2+1) ~= 0 || loadvals(2+3) ~= 0
        %~ error('only local pressure implemented. create an issue!')
        %~ end

        %~ else
        %~ error(' local/global load param must be 1 or 0')

      end % if global/local system

      elemNodeLoadsMatrix = ones( length(nodes), 1 )*[Fx 0 Fy 0 Fz 0] ;

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
