% function that constructs the assembled Fext vector for one BCtype 

function fext = elem2NodalLoads ( Conec, indBC, elements, boundaryConds, Nodes )
  
  elemsWithBC = find( Conec(:,3) == indBC ) ;
  
  loadedNodes = []                ;
  nnodes      = size(Nodes,1)     ;
  fext        = zeros(6*nnodes,1) ;
  
  for elemInd = 1:length( elemsWithBC )
  
    elem            = elemsWithBC( elemInd )               ;
    elemType        = elements.elemType{ Conec(elem,2 )}   ;
    loadCoordSys    = boundaryConds.loadsCoordSys{ indBC } ;

    % nodal loads
    % -----------
    if strcmp( elemType, 'node') % node
  
      if strcmp( loadCoordSys, 'global' )
    
        loadvals = boundaryConds.loadsBaseVals{ indBC } ;
        nodes     = Conec( elem, 4+1 ) ;

      else
        error(' only global flag in load by now.');
      end

      % add loads to matrix of loaded nodes
      loadedNodes = [ loadedNodes ; ...
                      nodes loadvals ] ;


    elseif strcmp( elemType , 'truss') ; %
      error(' not yet.');

    elseif strcmp( elemType , 'frame') ; %
      error(' not yet.');  

    % triangle tension
    % ----------------
    elseif strcmp( elemType , 'triangle') ; %
  
      nodes = Conec( elem, 4+(1:3) ) ;

      area = 0.5 * norm( cross( ...
        Nodes( nodes(2),:) - Nodes( nodes(1),:) , ...
        Nodes( nodes(3),:) - Nodes( nodes(1),:) ...
        ) ) ;

      loadvals = boundaryConds.loadsBaseVals{ indBC } ;

      if strcmp( loadCoordSys, 'global' )
  
        Fx = loadvals( 1 ) * area / 3 ;
        Fy = loadvals( 3 ) * area / 3 ;
        Fz = loadvals( 5 ) * area / 3 ;
      
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

      % add loads to matrix of loaded nodes
      loadedNodes = [ loadedNodes ; ...
                      nodes' ones(length(nodes),1)*[Fx 0 Fy 0 Fz 0] ] ;

    
  
    end %if elemTypes
    


  end % for elements
  

  % convert to assembled fext vetor
  % -------------------------------
  if exist( 'loadedNodes' ) ~= 0
    for i=1:size( loadedNodes ,1)
      aux = nodes2dofs ( loadedNodes(i,1), 6 ) ;
      fext( aux ) = fext( aux ) + loadedNodes(i,2:7)' ;
    end
  end
