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

%md this function creates the point and cell data matrices for writing
%md vtk files of the solids or structures.

function [ vtkNodes, vtkConec, vtkPointDataCell, vtkCellDataCell ] = vtkDataConversion( modS, modP )

fprintf('\n\n00000000000\n\n\n')

  elemTypeInds = unique( modP.Conec( :, 2 ) )

  nelems = size( modP.Conec, 1 )

  %md create empty matrices for node coordinates and cell connectivities
  vtkConec = [] ;
  vtkNodes = [] ;

  %md create empty matrices which are filled (if corresponding data is present)
  vtkNodalDisps  = [] ;
  vtkNormalForce = [] ;

  %md loop in element types and add nodes and cells considering the specific
  %md structure and connectivity for each type of element/cell
  %md

  for indType = 1:length( elemTypeInds )
    elemTypeString = modP.elements( elemTypeInds(indType) ).elemType
    elemTypeGeom   = modP.elements( elemTypeInds(indType) ).elemTypeGeometry

    elemIndsElemType = find( modP.Conec(:,2)==elemTypeInds(indType) )

    if strcmp( elemTypeString, 'truss' )

       [ trussVtkNodes, trussVtkConec] ...%, trussVtkNodalDisps, trussVtkNormalForces ] ...
         = trussVtkData( modP.Nodes, modP.Conec( elemIndsElemType, 5:end ), elemTypeGeom, modS.U );
    end
  end % for: elemTypeInds


  vtkNodes = trussVtkNodes ;
  vtkConec = trussVtkConec ;



vtkPointDataCell = {};
vtkCellDataCell = {} ;




% cell data : usar normalForces




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convierte las vigas en sólidos para mostrar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ vtkNodes, vtkConec] ...%, vtkDisps, vtkNormalForces ]
   = trussVtkData( Nodes, Conec, elemTypeGeom, U )

  vtkNodes = [] ;
  vtkConec = [] ;

  nelem = size(Conec,1)

  for i=1:nelem
    i
    nodesElem = Conec(i,1:2)
    dofsElem  = nodes2dofs( nodesElem, 6 )
    coordSubElem = Nodes(nodesElem(:),:)

    [ iniNodes, midNodes, endNodes, sectPar ] = getVtkConnecForCrossSec( elemTypeGeom )

    [ NodesCell, ConecCell ] = vtkBeam2SolidConverter ( coordSubElem, sectPar, U(dofsElem), eye(3) )

    lastNode = size( vtkNodes, 1 ) ;
    ConecCell(:,2:end) = ConecCell(:,2:end) + lastNode ;

    vtkNodes = [ vtkNodes ; NodesCell ] ;
    vtkConec = [ vtkConec ; ConecCell ] ;
  end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ iniNodes, midNodes, endNodes, secc ] = getVtkConnecForCrossSec( elemCrossSecParams )

  if elemCrossSecParams(1) == 1  || elemCrossSecParams(1) == 2 % general or rectangular section
    if elemCrossSecParams(1) == 1
      % equivalent square section using A = wy * wz
      auxh = sqrt( elemCrossSecParams(2) ) ;   auxb = auxh ;
      secc = [ 12 auxb auxh ] ;
    else
      secc = [ 12 elemCrossSecParams(2) elemCrossSecParams(3) ] ;
    end

		iniNodes = [ 1 2 3 4  ] ;
		midNodes = [          ] ;
    endNodes = [ 5 6 7 8  ] ;

  elseif elemCrossSecParams(1) == 3 % circular section

    secc = [ 25 elemCrossSecParams(2) ] ;
		iniNodes = [ 1 2 3 4 9 10 11 12  ] ;
		midNodes = [ 17 18 19 20         ] ;
    endNodes = [ 5 6 7 8 13 14 15 16 ] ;
  end
