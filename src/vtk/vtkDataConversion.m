% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
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

%md This function creates the point and cell data matrices for writing
%md vtk files of the solids or structures.

function [ vtkNodes, vtkConec, vtkPointDataCell, vtkCellDataCell ] = vtkDataConversion( modS, modP )

  nelems       = size( modP.Conec, 1 ) ;

  % get the indexes of elements used in the model
  elemTypeInds = unique( modP.Conec( :, 2 ) ) ;

  %md create empty matrices for node coordinates and cell connectivities
  vtkConec = [] ;             vtkNodes = [] ;
  %md create cells for point and cell data (displacemets)
  vtkPointDataCell = {} ;     vtkCellDataCell  = {} ;
  vtkNodalDisps         = [] ;  % after filling the matrix it is assigned to the third column of the cell

  %md loop in element types and add nodes and cells considering the specific
  %md structure and connectivity for each type of element/cell
  %md
  for indType = 1:length( elemTypeInds )

    elemTypeString = modP.elements( elemTypeInds(indType) ).elemType         ;
    elemTypeGeom   = modP.elements( elemTypeInds(indType) ).elemTypeGeometry ;

    % gets all the element numbers corresponding to the current elemType
    elemIndsElemType = find( modP.Conec(:,2)==elemTypeInds(indType) ) ;

    if strcmp( elemTypeString, 'truss' )

      [ currVtkNodes, currVtkConec, currVtkNodalDisps, vtkNormalForces ] ...
        = trussVtkData( modP.Nodes, modP.Conec( elemIndsElemType, 5:end ), ...
        elemTypeGeom, modS.U ) ;

    elseif strcmp( elemTypeString, 'frame' )

      [ currVtkNodes, currVtkConec, currVtkNodalDisps, vtkNormalForces ] ...
        = frameVtkData( modP.Nodes, modP.Conec( elemIndsElemType, 5:end ), ...
        elemTypeGeom, modS.U ) ;

    elseif strcmp( elemTypeString, 'tetrahedron' )


    end % if: type

    % add entries from current element type
    vtkNodes      = [ vtkNodes ;  currVtkNodes ]           ;
    vtkConec      = [ vtkConec ;  currVtkConec ]           ;
    vtkNodalDisps = [ vtkNodalDisps ;  currVtkNodalDisps ] ;

  end % for: elemTypeInds

  if length( vtkNodalDisps ) > 0
    vtkPointDataCell{1,1} = 'VECTORS'       ;
    vtkPointDataCell{1,2} = 'Displacements' ;
    vtkPointDataCell{1,3} = vtkNodalDisps ;
  end
