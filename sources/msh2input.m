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


% Function that process msh_read

function [Nodes, Conec, nodalVariableLoads, nodalSprings] = msh2input( mshfile, p )
  
[ Nodes , node_elem , line_elem , trng_elem, tet_elem ] = solid3DmshRead( mshfile ) ;

ndofpnode = 6 ;  
nnodes = size(Nodes,1) ;


nodalVariableLoads = zeros(nnodes,7) ;
nodalSprings = zeros(nnodes,7) ;
for i = 1:size(trng_elem,1)
  
  nodestrng = trng_elem(i,2:end) ;

  if trng_elem(i,1)<200
    
    if trng_elem(i,1) == 101
      nodalSprings(nodestrng,:) = [nodestrng' ones(3,1)*inf zeros(3,1) ones(3,1)*inf zeros(3,1) ones(3,1)*inf zeros(3,1)] ;
    elseif trng_elem(i,1) == 102
      nodalSprings(nodestrng,[1 2]) = [nodestrng' ones(3,1)*inf ] ;      
    elseif trng_elem(i,1) == 103  
      nodalSprings(nodestrng,[1 4]) = [nodestrng' ones(3,1)*inf ] ;  
    elseif trng_elem(i,1) == 104
      nodalSprings(nodestrng,[1 6]) = [nodestrng' ones(3,1)*inf ] ;  
    end
  
  else  
    
  area = 0.5*norm( cross( Nodes(nodestrng(2),:) - Nodes( nodestrng(1),:) , Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) ) ;
  dofsaux = nodes2dofs( nodestrng , ndofpnode ) ;
  dofs = dofsaux(1:2:length(dofsaux)) ;
  nmod = norm( cross( Nodes(nodestrng(2),:) - Nodes( nodestrng(1),:) , Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) ) ;
  n = cross( Nodes(nodestrng(2),:) - Nodes( nodestrng(1),:) , Nodes( nodestrng(3),:) - Nodes( nodestrng(1),: ) ) / nmod ;
  
  Fx = n(1)*p*area/3 ;
  Fy = n(2)*p*area/3 ;
  Fz = n(3)*p*area/3 ;
  
    if trng_elem(i,1) == 201
     
      sum = nodalVariableLoads(nodestrng,2:end) + [ ones(3,1)*Fx zeros(3,1) ones(3,1)*Fy zeros(3,1) ones(3,1)*Fz zeros(3,1) ] ;
      nodalVariableLoads(nodestrng,:) = [ nodestrng' sum ] ;
     
    elseif trng_elem(i,1) == 202
    
      sum = nodalVariableLoads(nodestrng,2:end) + [ ones(3,1)*Fx zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) ] ;
      nodalVariableLoads(nodestrng,:) = [ nodestrng' sum ] ;
     
    elseif trng_elem(i,1) == 203
      
      sum = nodalVariableLoads(nodestrng,2:end) + [ zeros(3,1) zeros(3,1) ones(3,1)*Fy zeros(3,1) zeros(3,1) zeros(3,1) ] ;
      nodalVariableLoads(nodestrng,:) = [ nodestrng' sum ] ;
     
    elseif trng_elem(i,1) == 204
    
      sum = nodalVariableLoads(nodestrng,2:end) + [ zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) ones(3,1)*Fz zeros(3,1) ] ;
      nodalVariableLoads = [ nodestrng' sum ] ;  
    end
  end
end  

null = find(nodalSprings(:,1)==0) ;
nodalSprings(null,:) = [] ;
null = find(nodalVariableLoads(:,1)==0) ;
nodalVariableLoads(null,:) = [] ;

Conec = [] ;

% Conec matrix for tetrahedron elements

Conec = [ tet_elem(:,3:end) ones(size(tet_elem,1),1) zeros(size(tet_elem,1),1) ones(size(tet_elem,1),1)*3 ] ; 

end  
