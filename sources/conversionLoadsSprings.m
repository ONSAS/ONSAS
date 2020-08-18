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

function [ Conec, nodalVariableLoads, nodalConstantLoads, nodalSprings ] = conversionLoadsSprings ( Nodes, Conec, materialsParams, ...
                        elementsParams, ...
                        loadsParams, ...
                        crossSecsParams, ...
                        springsParams, ...
                        BooleanSelfWheight ...
                        )

% auxiliar elements separation
indsElemsLoad      = find( Conec(:,4+3) ~= 0 ) ;
indsElemsSpri      = find( Conec(:,4+5) ~= 0 ) ;

nodalConstantLoads = [] ;
nodalVariableLoads = [] ;
nodalSprings       = [] ;

indsLoop = unique( [ indsElemsLoad; indsElemsSpri ] ) ;

% looop in auxiliar elements ( with material = 0 ) 
for ind = 1:length(indsLoop)
  
  i = indsLoop ( ind ) ;

  elemNum  = Conec( i, 4 + 2 ) ;
  loadNum  = Conec( i, 4 + 3 ) ;
  crosNum  = Conec( i, 4 + 4 ) ;
  spriNum  = Conec( i, 4 + 5 ) ;

  % --------------------------------------------------------------------
  % --- loads ---
  if loadNum > 0,
  
    if elementsParams{ elemNum } == 1 ; % node
      
      loadvals = loadsParams{loadNum} ;
      node     = Conec( i, 1 ) ; 

      % ------------------------------------------------------------------    
      if loadvals(1) ~= 1 
        error(' wrong global/local flag in load');
      end

      if loadvals(2) == 1 % variable load
        nodalVariableLoads = [ nodalVariableLoads ; ...
                             node loadvals(3:end) ] ;

      elseif loadvals(2) == 0 % constant load
        nodalConstantLoads = [ nodalConstantLoads ; ...
                             node loadvals(3:end) ] ;
      
      else
        error(' constant/variable load param must be 1 or 0')
      end
      % ------------------------------------------------------------------    


    elseif elementsParams{ elemNum } == 5 ; % triangle

      nodestrng = Conec( i, 1:3 ) ; 

      area = 0.5 * norm( cross( ...
        Nodes( nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
        Nodes( nodestrng(3),:) - Nodes( nodestrng(1),:) ...
        ) ) ;
      
      loadvals = loadsParams{loadNum} ;

      if loadvals(1) == 1 % global coordinates load

        Fx = loadvals(loadNum, 2+1) * area / 3 ;
        Fy = loadvals(loadNum, 2+3) * area / 3 ;
        Fz = loadvals(loadNum, 2+5) * area / 3 ;
        
      elseif loadvals(1) == 0 % local coordinates load

        dofsaux = nodes2dofs( nodestrng , 6 ) ;
        dofs    = dofsaux(1:2:length(dofsaux)) ;
        nmod    = norm( cross( ...
          Nodes( nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
          Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) ) ;

        n = cross( ...
          Nodes(nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
          Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) / nmod ;

        Fx = n(1) * loadvals(loadNum, 2+5 ) * area/3 ;
        Fy = n(2) * loadvals(loadNum, 2+5 ) * area/3 ;
        Fz = n(3) * loadvals(loadNum, 2+5 ) * area/3 ;
        
        if loadvals(2+1) ~= 0 || loadvals(2+3) ~= 0
          error('only local pressure implemented. create an issue!')
        end

      else
        error(' local/global load param must be 1 or 0')
      end
      
      if loadvals(2) == 1 % variable load
        nodalVariableLoads = [ nodalVariableLoads ; ...
                             nodestrng' ones(3,1)*[Fx 0 Fy 0 Fz 0] ] ;

      elseif loadvals(2) == 0 % constant load
        nodalConstantLoads = [ nodalConstantLoads ; ...
                             nodestrng' ones(3,1)*[Fx 0 Fy 0 Fz 0] ] ;
      
      else
        error(' constant/variable load param must be 1 or 0')
      end


%~ % add nodal loads
%~ for i = 1:nnodes

  %~ loa = nodesMat(i,4) ;
  %~ if loa > 0
		%~ if loa == 1
			%~ if loadsMat(loa, 1) == 1
				%~ nodalConstantLoads = [ nodalConstantLoads ; ...
															 %~ i loadsMat(loa, 2:7) ] ;
			%~ else
				%~ ms = msgbox('local node considered in node.') ;
			%~ end
    %~ elseif loa == 2
			%~ if loadsMat(loa, 1) == 1
				%~ nodalVariableLoads = [ nodalVariableLoads ; ...
															 %~ i loadsMat(loa, 2:7) ] ;
			%~ else
				%~ ms = msgbox('local node considered in node.') ;
			%~ end
    %~ end
  %~ end     

  %~ sup = nodesMat(i,5) ;
  %~ if sup > 0,
    %~ nodalSprings = [ nodalSprings ; ...
                     %~ i suppsMat(sup, :) ] ;
  %~ end     
%~ end

    
      %~ length = norm( Nodes( nodesElem(2),:) - Nodes( nodesElem(1),:) ) ;

      %~ if loadsMat(loa, 1) == 1  % global coordinates load

        %~ Fx = loadsMat(loa, 2) * length / 2 ;
        %~ Fy = loadsMat(loa, 4) * length / 2 ;
        %~ Fz = loadsMat(loa, 6) * length / 2 ;

        %~ nodalVariableLoads = [ nodalVariableLoads ; ...
                             %~ nodesElem' ones(2,1)*[Fx 0 Fy 0 Fz 0] ] ;
      %~ else
        %~ error('option not implemented, please create an issue.')
      %~ end
     
    end
  end

  % --------------------------------------------------------------------

 % --- SefWheight ---
 
     
     % --------------------------------------------------------------------
  % --- springs ---
  if spriNum > 0,

    nodesElem    = nonzeros( Conec (i, 1:4 ) ) ;
    nodalSprings = [ nodalSprings ; ...
                     nodesElem ones(size(nodesElem))*springsParams{ spriNum } ] ;
    
    if ( elementsParams{elemNum} == 3 )
      if   ( sum( springsParams{ spriNum } < inf ) > 0 ) ...
        && ( sum( springsParams{ spriNum } > 0   ) > 0 )
        
        warning('nodals springs in triangles included as absolut values. create an issue!\n')
      end % area spring considered
    end % triangle if
  end % spring if
  % --------------------------------------------------------------------

end

 if BooleanSelfWheight==1
     indsElems      = find( Conec(:,2) ~= 0 ) ; 
     NumElems       = length(indsElems);
    
    for i = indsElems(1):indsElems(end)
  
      matNum   = Conec( i, 4 + 1 ) ;
      elemNum  = Conec( i, 4 + 2 ) ;
      crosNum  = Conec( i, 4 + 4 ) ;
 
  
     if elementsParams{ elemNum }(1) == 2 ; % truss

            nodestrng = Conec( i, 1:2 ) ;       %element Nodes
            xelem     = Nodes(nodestrng,:);     
            Lelem     = norm( xelem(1,:)-xelem(2,:));%element length
            crossSecsParamsElem =crossSecsParams {crosNum}; %element cross sec params
            if crossSecsParamsElem (1)==1                               %DUDA DE COMO CORREGIR SI NO ES UNO, DONDE SE GUARDA EL AREA EN LA NUEVA NOMENCLATURA?
                Areaelem = pi*crossSecsParamsElem(2)^2/4;
             elseif crossSecsParamsElem (1)==2
                Areaelem = crossSecsParamsElem(2)*crossSecsParamsElem(3);
             else
                Areaelem = crossSecsParamsElem (2)
            end
            Matelem   = materialsParams {matNum} ;   %element material vec
            rhoelem   = Matelem(1);                  %densidad del elemento

            Fz = rhoelem*Lelem*Areaelem*9.8/2;
        
            nodalConstantLoads = [ nodalConstantLoads ; ...
                nodestrng', ones(2,1)*[0 0 0 0 -Fz 0]]; 
            %No habria que sumarlas a nodal ConstaLoads?
       end     
    end 
  end

indsElemsAux             = find( Conec(:,4+1) == 0 ) ;
conecAuxElems            = Conec(indsElemsAux, : ) ;
Conec( indsElemsAux, : ) = [] ;







      



    
			%~ if loa > 0
				%~ elem = i ;
				%~ axi = loadsMat(i,1) ;
				%~ qx = loadsMat(i,2) ;
				%~ qy = loadsMat(i,4) ;
				%~ qz = loadsMat(i,6) ;
				%~ unifLoad = [unifLoad ; elem axi qx qy qz] ;
			%~ end
			
  
    %~ if trng_elem(i,1) == 201
     
      %~ sum = nodalVariableLoads(nodestrng,2:end) + [ ones(3,1)*Fx zeros(3,1) ones(3,1)*Fy zeros(3,1) ones(3,1)*Fz zeros(3,1) ] ;
      %~ nodalVariableLoads(nodestrng,:) = [ nodestrng' sum ] ;
     
    %~ elseif trng_elem(i,1) == 202
    
      %~ sum = nodalVariableLoads(nodestrng,2:end) + [ ones(3,1)*Fx zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) ] ;
      %~ nodalVariableLoads(nodestrng,:) = [ nodestrng' sum ] ;
     
    %~ elseif trng_elem(i,1) == 203
      
      %~ sum = nodalVariableLoads(nodestrng,2:end) + [ zeros(3,1) zeros(3,1) ones(3,1)*Fy zeros(3,1) zeros(3,1) zeros(3,1) ] ;
      %~ nodalVariableLoads(nodestrng,:) = [ nodestrng' sum ] ;
     
    %~ elseif trng_elem(i,1) == 204
    
      %~ sum = nodalVariableLoads(nodestrng,2:end) + [ zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) ones(3,1)*Fz zeros(3,1) ] ;
      %~ nodalVariableLoads = [ nodestrng' sum ] ;  
    %~ end

  %~ end

%~ end  

%~ null = find(nodalSprings(:,1)==0) ;
%~ nodalSprings(null,:) = [] ;
%~ null = find(nodalVariableLoads(:,1)==0) ;
%~ nodalVariableLoads(null,:) = [] ;

%~ Conec = [] ;

%~ % Conec matrix for tetrahedron elements

%~ Conec = [ tet_elem(:,3:end) ones(size(tet_elem,1),1) zeros(size(tet_elem,1),1) ones(size(tet_elem,1),1)*3 ] ; 



