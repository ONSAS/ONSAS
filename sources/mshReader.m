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

function [ nodesMat, conecMat ] = mshReader( mshfile )

fid = fopen( mshfile ,'r') ;

% ---- header reading --------------------------
X = fgets(fid);
X = fgets(fid);
if strncmp( X, '2.2',3)
   X = fgets(fid);  X = fgets(fid);
else
  error('wrong format of msh mesh! Gmsh legacy 2 format expected. \n');
end
% ------------------------------


% --- reads physical names if they are defined -------
if strncmp( X, '$Physic',5)
  nPhysicalNames = fscanf(fid,'%g',[1 ]) ;
  numsPhysProp = zeros( nPhysicalNames, 1 ) ;
  strsPhysProp = cell ( nPhysicalNames, 1 ) ; 

  for i=1:nPhysicalNames
    numsPhysProp(i) = fscanf(fid,'%g %g',[2 1])(2) ;  
    strsPhysProp(i) = fscanf(fid,'%s ',[1 1])      ;   
  end

  X = fgets(fid) ;
  X = fgets(fid) ;
end
% ----------------------------------------------------


% ------- read nodes coordinates -------
nnodes = fscanf(fid,'%g',[1 ]) ;
nodesMat = fscanf(fid,'%g %g %g %g',[4 nnodes])' ;
nodesMat(:,1) = [] ;
% --------------------------------------
nodesMat(:,4:5) = 0 ;

X = fgets(fid) ; X = fgets(fid) ; X = fgets(fid) ;

nelem = fscanf(fid,'%g',[1 ]) ;

a = fgets(fid,30) ;

conecMat = zeros(nelem, 9 );

auxPhysNodes = [] ;

for i=1:nelem
  a = str2num(fgets(fid,200)) ;

  type_elem = a(2) ;
  ntags     = a(3) ;
  phy_ent   = a(4) ;

  [bool, ind] = ismember( phy_ent , numsPhysProp ) ;

  if  type_elem == 15
    if bool
      auxstr = strsPhysProp{ind} ;
      auxnodes = a(3+ntags+1) ;
      
      nodesMat(auxnodes,4) = str2num(auxstr(14:15)) ;
      nodesMat(auxnodes,5) = str2num(auxstr(11:12)) ;
      
      auxPhysNodes = [ auxPhysNodes	; auxnodes] ;
    end
     
  else

    if bool
      auxstr = strsPhysProp{ind} ;

      if type_elem == 1
        auxnodes  = [ a(3+ntags+(1:2)) ] ;
        auxnodes2 = [ auxnodes 0 0 ] ;
      
      elseif type_elem == 2
        auxnodes  = [ a(3+ntags+(1:3))   ] ;
        auxnodes2 = [ auxnodes 0 ] ;
      
      elseif type_elem == 4
        auxnodes = [ a(3+ntags+(1:4)) ] ;
        auxnodes2 = auxnodes ;
      else
        error('missing type')
      end
            
      conecMat( i, 1:5) = [ str2num(auxstr(2:3)) auxnodes2 ] ;
      conecMat( i, 6:9) = [ str2num(auxstr( 5: 6)) str2num(auxstr( 8: 9)) ...
                            str2num(auxstr(14:15)) str2num(auxstr(11:12)) ] ;
			for j = 1:length(auxnodes)
				if ~ismember(auxnodes(j), auxPhysNodes)
					nodesMat(auxnodes(j),4:5) = ones(length(auxnodes(j)),2) * diag( [str2num(auxstr(14:15)) str2num(auxstr(11:12)) ] ) ;
        end
      end
                            
    end
        
  end

end

fclose(fid);


