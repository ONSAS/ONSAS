function [externalLoadsCases, loadTimeFactors ] = elem2NodalLoads ( externalLoadsCases, loadTimeFactors, boundaryConds, elements, Conec, elem )

  elemTypeNum = Conec( elem, 2 ) ;
  bcNum       = Conec( elem, 3) ;
  
  % nodal loads
  % -----------
  if strcmp( elements.elemType{ elemTypeNum } , 'node') ; % node

    if strcmp( boundaryConds.loadCoordSys{ bcNum },'global' )
    
      loadvals = boundaryConds.loadBaseVals{ bcNum } ;
      node     = Conec( elem, 4+1 ) ; 

    else
      error(' only global flag in load by now.');
    end

  
  % nodal loads
  % -----------
  if strcmp( elements.elemType{ elemTypeNum } , 'triangle') ; % node

  
  
  
  
  
  
  end


    nodalVariableLoads = [ nodalVariableLoads ; ...
                           %~ node loadvals(3:end) ] ;
  
    %~ elseif loadvals(2) == 0 % constant load
      %~ nodalConstantLoads = [ nodalConstantLoads ; ...
                           %~ node loadvals(3:end) ] ;
    
    %~ else
      %~ error(' constant/variable load param must be 1 or 0')
    %~ end
    %~ % ------------------------------------------------------------------    
  
  
  %~ elseif elementsParams{ elemNum } == 5 ; % triangle
  
    %~ nodestrng = Conec( i, 1:3 ) ; 
  
    %~ area = 0.5 * norm( cross( ...
      %~ Nodes( nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
      %~ Nodes( nodestrng(3),:) - Nodes( nodestrng(1),:) ...
      %~ ) ) ;
    
    %~ loadvals = loadsParams{loadNum} ;
  
    %~ if loadvals(1) == 1 % global coordinates load
  
      %~ Fx = loadvals(loadNum, 2+1) * area / 3 ;
      %~ Fy = loadvals(loadNum, 2+3) * area / 3 ;
      %~ Fz = loadvals(loadNum, 2+5) * area / 3 ;
      
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
    %~ end
    
    %~ if loadvals(2) == 1 % variable load
      %~ nodalVariableLoads = [ nodalVariableLoads ; ...
                           %~ nodestrng' ones(3,1)*[Fx 0 Fy 0 Fz 0] ] ;
  
    %~ elseif loadvals(2) == 0 % constant load
      %~ nodalConstantLoads = [ nodalConstantLoads ; ...
                           %~ nodestrng' ones(3,1)*[Fx 0 Fy 0 Fz 0] ] ;
    
    %~ else
      %~ error(' constant/variable load param must be 1 or 0')
    %~ end
         
  %~ end % if type elem
  
  
