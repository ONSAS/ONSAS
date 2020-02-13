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



% === function for reading f and x values ===
% first author: J. Perez Zerpa
% revision: april 2010
% lectura de archivos
% interpretacion de condiciones geometricas y fisicas
%
% %$MeshFormat
% version-number file-type data-size
% %$EndMeshFormat
% %$Nodes
% number-of-nodes
% node-number x-coord y-coord z-coord
% ...
% %$EndNodes
% %$Elements
% number-of-elements
% elm-number elm-type number-of-tags < tag > ... node-number-list 
%
% number-of-tags
% gives the number of integer tags that follow for the n-th element. By default, the
% first tag is the number of the physical entity to which the element belongs; the
% second is the number of the elementary geometrical entity to which the element
% belongs; the third is the number of a mesh partition to which the element
% belongs. All tags must be postive integers, or zero. A zero tag is equivalent to
% no tag.
%
% ...
% %$EndElements
% %$PhysicalNames
% number-of-names
% physical-dimension physical-number "physical-name"
% ...
% %$EndPhysicalNames
%
%
%
% geometrical properties
% point
% 1  todo libre    no considerado
% 2  retengo el x   
% 3 retengo el y
% 4 retengo x e y
%
% line
% 1   libre
% 2  retengo x
% 3 retengo y
% 4 retengo x e y
%
% surface
% 

function [ Nodes , node_elem , line_elem , trng_elem, tet_elem ] = solid3DmshRead( mshfile )

fid = fopen( mshfile ,'r') ;

%f = fscanf(fid , '%c',[ 5 , 20  ] ) ;
X = fgets(fid);   % It has two rows now.x
X = fgets(fid);   % It has two rows now.
X = fgets(fid);   % It has two rows now.
X = fgets(fid);   % It has two rows now.
% ----------------------------------------


% --- number of nodes ---
nnodes = fscanf(fid,'%g',[1 ]) ;  % It has two rows now.

% --- nodes coordinates ---
Nodes = fscanf(fid,'%g %g %g %g',[4 nnodes])' ;   % lee nodos.
%nodes = aux' ;%  descomentar para eliminar numero nodo '(:,2:end) ; % matriz de nodos

%nodes
size(Nodes,1) ;


X = fgets(fid) ; X = fgets(fid) ; X = fgets(fid) ;

nelem = fscanf(fid,'%g',[1 ]) ;

a=fgets(fid,30);

n_node = 0 ;
n_line = 0 ;
n_trng = 0 ;
n_tet  = 0 ;

node_elem = [] ;
line_elem = [] ;
trng_elem = [] ;
tet_elem = [] ;


% $ELM
% numero elemento elem |  dimensao (1 segmento 2 poligono) |  fisical lines or surfaces  |   lineas or surfaces geometricas  |  numero de noos   |  nodo1  | nodo 2 | nodo 3  |
% $ENDELM

for i=1:nelem
   %x = fscanf(fid,'%g %g %g %g %g %g %g %g ',[8 1])
   a = str2num(fgets(fid,200)) ;
      
   type_elem = a(2) ;
   ntags     = a(3) ;

   %~ if abs(rem(i/1000,1))>0
     %~ i
     %~ end
   if ntags == 3
      phy_ent = a(4) ;
      geo_ent = a(5) ;
      partitn = a(6) ;
      
      if  type_elem == 15
         n_node = n_node + 1 ;
         node_elem (n_node,:) = [ n_node  phy_ent  a(7) ] ;
   
      elseif type_elem == 1
         n_line = n_line + 1 ;
         line_elem (n_line , :) = [ n_line phy_ent a(7:8) ]  ;
         
      elseif type_elem == 2
         n_trng = n_trng + 1 ;
         trng_elem (n_trng , :) = [ n_trng phy_ent a(7:9) ] ;
         
      end

   elseif ntags == 2 
      phy_ent = a(4) ;
      geo_ent = a(5) ;
      
      if  type_elem == 15
         n_node = n_node + 1 ;
         node_elem (n_node,:) = [ n_node  phy_ent  a(6) ] ;
   
      elseif type_elem == 1
         n_line = n_line + 1 ;
         line_elem (n_line , :) = [ n_line phy_ent a(6:7) ]  ;
         
      elseif type_elem == 2
         n_trng = n_trng + 1 ;
         trng_elem (n_trng , :) = [ n_trng phy_ent a(6:8) ] ;

      elseif type_elem == 4
         n_tet = n_tet + 1 ;
         tet_elem (n_tet , :) = [ n_tet phy_ent a(6:9) ] ;
         
      end
	  
   else
	  error('insert more tags cases')
   end

end


Nodes      (:,1) = [] ;
%~ node_elem  (:,1) = [] ;
%~ line_elem  (:,1) = [] ;
trng_elem  (:,1) = [] ;
%~ tet_elem   (:,1) = [] ;


fclose(fid);


