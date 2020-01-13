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


% Function that converts from the gui format to the .m input format

function [Nodes, Conec, nodalVariableLoads, nodalConstantLoads, unifDisLoadL, unifDisLoadG, nodalSprings ] = inputFormatConversion ( nodesMat, conecMat, loadsMat, suppsMat )
  
ndofpnode = 6 ;  

nnodes = size(nodesMat,1) ;
nelems = size(conecMat,1) ;


Nodes = nodesMat(:,1:3) ;

Conec = [];


nodalConstantLoads = [] ;
nodalVariableLoads = [] ;

unifDisLoadL       = [] ;
unifDisLoadG       = [] ;

nodalSprings       = [] ;


% add nodal loads
for i = 1:nnodes

  loa = nodesMat(i,4) ;
  if loa > 0
		if loa == 1
			if loadsMat(loa, 1) == 1
				nodalConstantLoads = [ nodalConstantLoads ; ...
															 i loadsMat(loa, 2:7) ] ;
			else
				ms = msgbox('local node considered in node.') ;
			end
    elseif loa == 2
			if loadsMat(loa, 1) == 1
				nodalVariableLoads = [ nodalVariableLoads ; ...
															 i loadsMat(loa, 2:7) ] ;
			else
				ms = msgbox('local node considered in node.') ;
			end
    end
  end     

  sup = nodesMat(i,5) ;
  if sup > 0,
    nodalSprings = [ nodalSprings ; ...
                     i suppsMat(sup, :) ] ;
  end     
end




for i = 1:nelems

  type = conecMat(i,1) ;
  sec  = conecMat(i,6) ;
  mat  = conecMat(i,7) ;
  loa  = conecMat(i,8) ;
  sup  = conecMat(i,9) ;
  
  switch type

  case 5
    nodestrng = conecMat(i,2:4) ;

    if sup > 0,
      nodalSprings = [ nodalSprings ; ...
                       nodestrng' ones(3,1)*suppsMat(sup, :) ] ;
    end

    if loa > 0,

      area = 0.5 * norm( cross( ...
        Nodes( nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
        Nodes( nodestrng(3),:) - Nodes( nodestrng(1),:) ...
        ) ) ;

      if loadsMat(loa, 1) == 1

        Fx = loadsMat(loa, 2)*area/3 ;
        Fy = loadsMat(loa, 4)*area/3 ;
        Fz = loadsMat(loa, 6)*area/3 ;
        
      elseif loadsMat(loa, 1) == 0

        dofsaux = nodes2dofs( nodestrng , ndofpnode ) ;
        dofs    = dofsaux(1:2:length(dofsaux)) ;
        nmod    = norm( cross( ...
          Nodes( nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
          Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) ) ;
        
        n = cross( ...
          Nodes(nodestrng(2),:) - Nodes( nodestrng(1),:) , ...
          Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) / nmod ;

        Fx = n(1)*loadsMat(loa,6)*area/3 ;
        Fy = n(2)*loadsMat(loa,6)*area/3 ;
        Fz = n(3)*loadsMat(loa,6)*area/3 ;
      
      end

      nodalConstantLoads = [ nodalConstantLoads ; ...
                           nodestrng' ones(3,1)*[Fx 0 Fy 0 Fz 0] ] ;

    end

    case 3 

      Conec = [Conec ; conecMat( i, 2:5 ) conecMat(i,7) conecMat(i,6) type ] ;
     
    case 2 

			Conec = [Conec ; conecMat( i, 2:5 ) conecMat(i,7) conecMat(i,6) type ] ;
			
		case 1	

			Conec = [Conec ; conecMat( i, 2:5 ) conecMat(i,7) conecMat(i,6) type ] ;
			
			%~ if loa > 0
				%~ elem = i ;
				%~ axi = loadsMat(i,1) ;
				%~ qx = loadsMat(i,2) ;
				%~ qy = loadsMat(i,4) ;
				%~ qz = loadsMat(i,6) ;
				%~ unifLoad = [unifLoad ; elem axi qx qy qz] ;
			%~ end
			
  end

end

  
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



%~ % collapse springs

%~ % sum loads
