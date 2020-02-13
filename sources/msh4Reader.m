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

% msh4Reader: function for reading gmsh's msh version 4.1 files.
% http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
%
% Input:
%   - mshFilename
%
% Output:
%  - nodesMat: matrix with 4 columns: [x y z physicalTag]
%  - conecMat: matrix with 5 columns: [ n1 n2 n3 n4 physicalTag ]
%              for elements with less than four nodes 0 is used as node.
%  - physicalNames: cell with strings of physical names.
%
% Assumptions:
%  - physical names are saved as strings
%  - maximum of one physical tag per entity 
%  - maximum number of nodes per element: 4 (linear tetrahedron)
%

function [ nodesMat, conecMat, physicalNames ] = msh4Reader( mshFilename )

fid = fopen( mshFilename ,'r') ;

maxLengthLine = 200 ;

% ---- header reading --------------------------
X = fgets(fid);
X = fgets(fid);
if strncmp( X, '4.1',3)
   X = fgets(fid);
   X = fgets(fid);  % read next header
else
  error('wrong format of msh mesh! Gmsh legacy 2 format expected. \n');
end
% ----------------------------------------------


% --- reads physical names if they are defined -------
if strncmp( X, '$Physic',5)
  nPhysicalNames = fscanf(fid,'%g',[1 ])      ;
  numsPhysProp   = zeros( nPhysicalNames, 1 ) ;
  physicalNames  = cell ( nPhysicalNames, 1 ) ; 

  for i=1:nPhysicalNames
    numsPhysProp(i) = fscanf(fid,'%g %g',[2 1])(2) ;  
    physicalNames{i} = fscanf(fid,'%s ',[1 1])     ;   
  end
  X = fgets(fid) ; % read end physical
  X = fgets(fid) ; % read next header

elseif nargout > 2
  error('no physical names header found in msh file.');
end
% ----------------------------------------------------


% --- reads entities if they are defined -------
if strncmp( X, '$Entiti',5)

  matsPhysicalPropsPerEntity = cell(4) ;

  entNumsPerDim = fscanf(fid,'%g %g %g %g',[4  1]) ; fgetl(fid);

  for i=1:4
    vecsPhysicalPropsPerEntity{i} = zeros( entNumsPerDim(i) , 1 ) ;
  end      

  for indDim = 1:4
    colNumTags  = 1+3+3*(indDim>1)+1 ;
    colTags     = 1+3+3*(indDim>1)+2 ;
    if entNumsPerDim(indDim) > 0
      for i=1:entNumsPerDim(indDim)
        aux = str2num( fgets(fid, maxLengthLine ) ) ; 
        if aux(colNumTags) > 0
          vecsPhysicalPropsPerEntity{indDim}(i) = aux( colTags ) ;
        end
      end
    end
  end

  X = fgets(fid) ; % read end entity header
  X = fgets(fid) ; % read next header
end
% ----------------------------------------------------


% --- reads entities if they are defined -------
if strncmp( X, '$Nodes',5)
  aux = fscanf(fid,'%g %g %g %g',[4  1]) ; fgetl(fid);
  
  numEntBlocks = aux(1) ;
  numNodes     = aux(2) ;

  nodesMat = zeros( numNodes, 4 ) ;
  
  for block = 1: numEntBlocks
    aux = fscanf(fid,'%g %g %g %g',[4  1]) ; fgetl(fid);
    if aux(4) > 0 % if there are nodes in the block
      nodesTags = fscanf(fid,'%g \n',[aux(4) 1]);
      matCoords = fscanf(fid,'%g \n',[3, aux(4)])';
      % only saves physical tags for nodes defined as nodes (no inheritance)
      auxphy = vecsPhysicalPropsPerEntity{aux(1)+1}(aux(2)) * ( aux(1) == 0) ;
      nodesMat( nodesTags,:) = [ matCoords ones(aux(4),1)*auxphy ] ;
    end
  end

  X = fgets(fid) ; % read end nodes header
  X = fgets(fid) ; % read next header
end
% ----------------------------------------------------



% --- reads entities if they are defined -------
if strncmp( X, '$Elements',5)
  aux = fscanf(fid,'%g %g %g %g',[4  1]) ; fgetl(fid);
  
  numEntBlocks = aux(1) ;
  numElems     = aux(2) ;

  conecMat = zeros( numNodes, 5 ) ;
  
  for block = 1: numEntBlocks
    aux = fscanf(fid,'%g %g %g %g\n',[4  1]) ; %fgetl(fid);
    
    if aux(4) > 0 % if there are elements in the block
      if aux(1)==2
        auxMatCon = fscanf(fid,'%g %g %g %g \n',[1+aux(1)+1, aux(4)] )' ;
      elseif aux(1)==3
        auxMatCon = fscanf(fid,'%g %g %g %g %g \n',[ 1+aux(1)+1 aux(4)] )' ;
      else
        error('dimension not implemented yet. Please create an issue.')
      end
      
      elemInds = auxMatCon(:,1) ;
      nodes    = auxMatCon(:,2:end) ; 
      auxphy = vecsPhysicalPropsPerEntity{aux(1)+1}(aux(2)) ;

      conecMat( elemInds,1:(aux(1)+1)) = nodes ;
      conecMat( elemInds,5           ) = auxphy ;
    end
    
  end

  X = fgets(fid) ; % read end elements header
  
end

fclose(fid);
