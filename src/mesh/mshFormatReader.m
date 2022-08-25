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
 
%md## mshFormatReader
%md function for reading gmsh's msh version 4.1 files.
%md http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
%md
%md Input:
%md   - mshFilename: the name of the msh file
%md Output:
%md  - nodesMat: matrix with 4 columns: [x y z physicalTag]
%md  - conecMat: matrix with 5 columns: [ n1 n2 n3 n4 physicalTag ]
%md              for elements with less than four nodes 0 is used as node.
%md  - physicalNames: cell with strings of physical names.
%md
%md Assumptions:
%md  - physical names are saved as strings
%md  - maximum of one physical tag per entity
%md  - maximum number of nodes per element: 4 (linear tetrahedron)
%md
function [ nodesMat, conecMat, physicalNames ] = mshFormatReader( mshFilename )

%md
fid = fopen( mshFilename ,'r') ;
maxLengthLine = 200 ;
%md### header reading
%md verifies if the msh version is 4.1
X = fgets( fid ) ;% read start header
assert( strncmp( fgets(fid), '4.1', 3 ), 'wrong format of msh mesh! Gmsh legacy 2 format expected. \n' )
X = fgets(fid);  % read end header
% ----------------------------------------------

%md### reads physical names
%md if they are defined
X = fgets( fid ); % read start next header
if strncmp( X, '$Physic',5)
  nPhysicalNames = fscanf(fid,'%g',[1 ])      ;
  numsPhysProp   = zeros( nPhysicalNames, 1 ) ;
  physicalNames  = cell ( nPhysicalNames, 1 ) ;

  for i=1:nPhysicalNames
    aux                              = fscanf(fid,'%g %g',[2 1]) ;
    numsPhysProp(i)                  = aux(2) ;
    physicalNames{ numsPhysProp(i) } = fscanf(fid,'%s ',[1 1])     ;
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
    vecsPhysicalPropsPerEntity{i} = zeros( entNumsPerDim(i) , 2 ) ;
  end

  %~ vecsPhysicalPropsPerEntity
%~ entNumsPerDim
  %~ stop
  for indDim = 1:4
    colNumTags  = 1+3+3*(indDim>1)+1 ;
    colTags     = 1+3+3*(indDim>1)+2 ;
    if entNumsPerDim(indDim) > 0
      for i=1:entNumsPerDim(indDim)
        aux = str2num( fgets(fid, maxLengthLine ) ) ;
        if aux(colNumTags) > 0
          vecsPhysicalPropsPerEntity{indDim}(i,:) = [ aux(1) aux( colTags ) ] ;
        else
          vecsPhysicalPropsPerEntity{indDim}(i,:) = [ aux(1) 0              ] ;
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
    aux = fscanf(fid,'%g %g %g %g',[4  1]) ;
    fgetl(fid);
    if aux(4) > 0 % if there are nodes in the block
      nodesTags = fscanf(fid,'%g \n',[aux(4) 1]) ;
      matCoords = fscanf(fid,'%g \n',[3, aux(4)])' ;

 % backup
      %~ % only saves physical tags for nodes defined as nodes (no inheritance)
      %~ auxphy = ( aux(1) == 0 ) * vecsPhysicalPropsPerEntity{aux(1)+1}(aux(2))  ;
      %~ nodesMat( nodesTags,:) = [ matCoords ones(aux(4),1)*auxphy ] ;
  % backup

  %~ block
  %~ aux
  %~ vecsPhysicalPropsPerEntity
  %~ stop
      % only saves physical tags for nodes defined as nodes (no inheritance)
      if aux(1) == 0
        vecsPhysicalPropsPerEntity{aux(1)+1} ;
        indEnt = find( vecsPhysicalPropsPerEntity{aux(1)+1}(:,1) == aux(2) ) ;
        auxphy = vecsPhysicalPropsPerEntity{aux(1)+1}( indEnt ,2 ) ;
      else
        auxphy = 0 ;
      end
      nodesMat( nodesTags,:) = [ matCoords ones(aux(4),1)*auxphy ] ;

    end
  end

  X = fgets(fid) ; % read end nodes header
  X = fgets(fid) ; % read next header
end
% ----------------------------------------------------



% --- reads elements if they are defined -------
if strncmp( X, '$Elements',5)
  aux = fscanf(fid,'%g %g %g %g',[4  1]) ; fgetl(fid);

  numEntBlocks = aux(1) ;
  numElems     = aux(2) ;

  conecMat = zeros( numNodes, 5 ) ;

  for block = 1: numEntBlocks
    aux = fscanf(fid,'%g %g %g %g\n',[4  1]) ; %fgetl(fid);

    if aux(4) > 0 % if there are elements in the block
      if     aux(1) == 0
        auxMatCon = fscanf(fid,'%g %g \n'       ,[ 1+aux(1)+1 aux(4)] )' ;
      elseif aux(1) == 1
        auxMatCon = fscanf(fid,'%g %g %g \n'       ,[ 1+aux(1)+1 aux(4)] )' ;
      elseif aux(1) == 2
        auxMatCon = fscanf(fid,'%g %g %g %g \n'    ,[ 1+aux(1)+1 aux(4)] )' ;
      elseif aux(1) == 3
        auxMatCon = fscanf(fid,'%g %g %g %g %g \n' ,[ 1+aux(1)+1 aux(4)] )' ;
      else
        aux
        error('dimension not implemented yet. Please create an issue.')
      end

      elemInds = auxMatCon(:,1)     ;
      nodes    = auxMatCon(:,2:end) ;
      indEnt = find( vecsPhysicalPropsPerEntity{aux(1)+1}(:,1) == aux(2) ) ;
      auxphy = vecsPhysicalPropsPerEntity{aux(1)+1}( indEnt ,2 ) ;
      %~ auxphy = vecsPhysicalPropsPerEntity{aux(1)+1}(aux(2)) ;

      conecMat( elemInds,1:(aux(1)+1)) = nodes ;
      conecMat( elemInds,5           ) = auxphy ;
    end

  end

  X = fgets(fid) ; % read end elements header

end

fclose(fid);
