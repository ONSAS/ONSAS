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
	totalNodes = 0 ;

  for indType = 1:length( elemTypeInds )

		elemTypeInds(indType) ;
    elemTypeString = modP.elements( elemTypeInds(indType) ).elemType       ;
    elemCrossSecParams   = modP.elements( elemTypeInds(indType) ).elemCrossSecParams ;

    % gets all the element numbers corresponding to the current elemType
    elemIndsElemType = find( modP.Conec(:,2)==elemTypeInds(indType) ) ;

    if strcmp( elemTypeString, 'node' )

      currVtkNodes = [] ;
      currVtkConec = [] ;
      currVtkNodalDisps = [] ;

    elseif strcmp( elemTypeString, 'truss' )

      [ currVtkNodes, currVtkConec, currVtkNodalDisps, vtkNormalForces ] ...
        = trussVtkData( modP.Nodes, modP.Conec( elemIndsElemType, 5:end ), ...
        elemCrossSecParams, modS.U ) ;

    elseif strcmp( elemTypeString, 'frame' )

      [ currVtkNodes, currVtkConec, currVtkNodalDisps, vtkNormalForces ] ...
        = frameVtkData( modP.Nodes, modP.Conec( elemIndsElemType, 5:end ), ...
        elemCrossSecParams, modS.U ) ;

    elseif strcmp( elemTypeString, 'triangle' )

      vtkNormalForces = [] ;
      % reshape the displacements vector
      currVtkNodalDisps = reshape( modS.U(1:2:end)', [3, size(modP.Nodes,1) ])' ;
      % and add it to the nodes matrix
      currVtkNodes   = modP.Nodes + currVtkNodalDisps ;
      vtkStress  = {} ;
      % vtkStress{3,1} = stressMat ;

      % if nargout > numminout
        % add the tetrahedron vtk cell type to the first column
        currVtkConec = [ 5*ones( nelems, 1 )     modP.Conec(:, 5:7 )-1 ] ;
        elem2VTKCellMap = (1:nelems)' ; % default: i-to-i . Columns should be changed.
      % end

    elseif strcmp( elemTypeString, 'tetrahedron' )

      vtkNormalForces = [] ;
      % reshape the displacements vector
      currVtkNodalDisps = reshape( modS.U(1:2:end)', [3, size(modP.Nodes,1) ])' ;
      % and add it to the nodes matrix
      currVtkNodes = modP.Nodes + currVtkNodalDisps ;

      vtkStress = {} ;
      % vtkStress{3,1} = stressMat ;

            % if nargout > numminout
              % add the tetrahedron vtk cell type to the first column
      currVtkConec    = [ 10*ones( nelems, 1 )     modP.Conec(:, 5:8 )-1 ] ;
      elem2VTKCellMap = (1:nelems)' ; % default: i-to-i . Columns should be changed.

    end % if: type

    % add entries from current element type
    vtkNodes      = [ vtkNodes ;  currVtkNodes ]           ;

    if size( vtkConec, 1 ) > 0
      vtkConec( ...
        (size(vtkConec,1)+1):(size(vtkConec,1)+size(currVtkConec,1)), 1:(size(currVtkConec,2)) ) ...
        = [currVtkConec(:,1) currVtkConec(:,2:end)+totalNodes] ;
    elseif size( currVtkConec, 1 ) > 0
      vtkConec = [currVtkConec(:,1) currVtkConec(:,2:end)+totalNodes] ;
    end

    vtkNodalDisps = [ vtkNodalDisps ;  currVtkNodalDisps ] ;

		totalNodes = totalNodes + size(currVtkNodes, 1) ;


  end % for: elemTypeInds



  if length( vtkNodalDisps ) > 0
    vtkPointDataCell{1,1} = 'VECTORS'       ;
    vtkPointDataCell{1,2} = 'Displacements' ;
    vtkPointDataCell{1,3} = vtkNodalDisps ;
  end

  stressMat =   modS.Stress ;


%   % Scalars vals
%   svm             = [] ;
%   vecSigI         = [] ;
%   vecSigII        = [] ;
%   vecSigIII       = [] ;
%   normalForceMat  = [] ;
%
%   % Vectors vals
%   vecI    = [] ;
%   vecII   = [] ;
%   vecIII  = [] ;
%
%   cellPointData   = cell(0) ;
%   cellCellData    = cell(0) ;
%   cellTensorData  = cell(0) ;
%
%   % the point data is assumed to be the displacements field
%   cellPointData = cell(1,3) ;
%   cellPointData{1,1} = 'VECTORS' ; cellPointData{1,2} = 'Displacements' ; cellPointData{1,3} = vtkDispMat ;
%
% return
%   if size( vtkNormalForce, 1 ) > 0
%     cellCellData{1,1} = 'SCALARS' ; cellCellData{1,2} = 'Normal_Force' ; cellCellData{1,3} = vtkNormalForce ;
%   end
%
   if size( stressMat, 1 ) > 0
     for indType = 1:length( elemTypeInds )
       elemTypeString = modP.elements( elemTypeInds(indType) ).elemType       ;
       if strcmp( elemTypeString, 'tetrahedron' )

         nelems    = size( stressMat, 1 ) ;

         sxx = stressMat(:, 1 ) ;
         syy = stressMat(:, 2 ) ;
         szz = stressMat(:, 3 ) ;
         tyz = stressMat(:, 4 ) ;
         txz = stressMat(:, 5 ) ;
         txy = stressMat(:, 6 ) ;

         svm = sqrt( 0.5*( (sxx-syy).^2 + (syy-szz).^2 + (szz-sxx).^2  + 6*(tyz.^2 + txz.^2 + txy.^2) ) ) ;
%
%     for i = 1:nelems
%       tensor = [ sxx(i) txy(i) txz(i) ; txy(i) syy(i) tyz(i) ; txz(i) tyz(i) szz(i) ] ;
%
%       [vec, val] = eig(tensor) ;
%
%       [ sigI  , indI   ] = max( diag(val) ) ;
%       [ sigIII, indIII ] = min( diag(val) ) ;
%
%       if indI ~= indIII
%         indII = setdiff( [1 2 3], [indI indIII]) ;
%       else
%         indII = indI ;
%       end
%
%       aux = diag(val) ;
%       sigII = aux(indII) ;
%
%       % Sigma I
%       vecSigI = [vecSigI ; sigI] ;
%       vecI = [vecI ; vec(:,indI)'/(norm(vec(:,indI)'))*abs(sigI)] ;
%
%       % Sigma II
%       vecSigII = [ vecSigII ; sigII ] ;
%       vecII = [ vecII ; vec(:,indII)'/(norm(vec(:,indII)'))*abs(sigII)] ;
%
%       % Sigma III
%       vecSigIII = [vecSigIII ; sigIII] ;
%       vecIII = [vecIII ; vec(:,indIII)'/(norm(vec(:,indIII)'))*abs(sigIII)] ;
%
%     end
         vtkCellDataCell{1,1} = 'SCALARS' ; vtkCellDataCell{1,2} = 'Von_Mises'   ; vtkCellDataCell{1,3} = svm ;
         vtkCellDataCell{2,1} = 'SCALARS' ; vtkCellDataCell{2,2} = 'Sxx'   ; vtkCellDataCell{2,3} = sxx ;
%     cellCellData{2,1} = 'VECTORS' ; cellCellData{2,2} = 'vI'          ; cellCellData{2,3} = vecI ;
%     cellCellData{3,1} = 'VECTORS' ; cellCellData{3,2} = 'vII'         ; cellCellData{3,3} = vecII ;
%     cellCellData{4,1} = 'VECTORS' ; cellCellData{4,2} = 'vIII'        ; cellCellData{4,3} = vecIII ;
%     cellCellData{5,1} = 'SCALARS' ; cellCellData{5,2} = 'SigI'        ; cellCellData{5,3} = vecSigI ;
%     cellCellData{6,1} = 'SCALARS' ; cellCellData{6,2} = 'SigII'       ; cellCellData{6,3} = vecSigII ;
%     cellCellData{7,1} = 'SCALARS' ; cellCellData{7,2} = 'SigIII'      ; cellCellData{7,3} = vecSigIII ;
%
%   elseif size( cellStress, 1 ) > 0 && size( cellStress{3}, 1) > 0
%     stressMat = cellStress{3} ;
%     nelems    = size( stressMat, 1 ) ;
%
%     sxx = stressMat(:, 1 ) ;
%
%     cellCellData{1,1} = 'SCALARS' ; cellCellData{1,2} = 'SigXLoc'   ; cellCellData{1,3} = sxx ;
%
       end
     end
   end
