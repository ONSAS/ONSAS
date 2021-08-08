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

%md this function creates the conectivities and coordinates data structures for producing
%md vtk files of the solids or structures.

function [ vtkNodes, vtkConec, cellPointData, cellCellData ] = vtkDataFormater( modelCurrSol, modelProperties)

numminout = 4 ;

Conec = modelProperties.Conec ;
Nodes = modelProperties.Nodes ;
U     = modelCurrSol.U        ;

elemTypeInds = unique( Conec(:,2) ) ;

nelems = size( Conec, 1 ) ;

vtkDispMat = [] ;

for indType = 1:length( elemTypeInds )

  % get the string of the current type of element
  elemTypeString = modelProperties.elements( elemTypeInds(indType) ).elemType ;

  if strcmp( elemTypeString, 'tetrahedron' )

    vtkNormalForces = [] ;
    % reshape the displacements vector
    vtkDispMat = reshape( U(1:2:end)', [3, size(Nodes,1) ])' ;
    % and add it to the nodes matrix
    vtkNodes = Nodes + vtkDispMat ;

    vtkStress = {} ;
    % vtkStress{3,1} = stressMat ;

    % if nargout > numminout
      % add the tetrahedron vtk cell type to the first column
      vtkConec = [ 10*ones( nelems, 1 )     Conec(:, 5:8 )-1 ] ;
      elem2VTKCellMap = (1:nelems)' ; % default: i-to-i . Columns should be changed.
    % end

  elseif strcmp( elemTypeString, 'triangle' )

    vtkNormalForces = [] ;
    % reshape the displacements vector
    vtkDispMat = reshape( U(1:2:end)', [3, size(Nodes,1) ])' ;
    % and add it to the nodes matrix
    vtkNodes = Nodes + vtkDispMat ;

    vtkStress = {} ;
    % vtkStress{3,1} = stressMat ;

    % if nargout > numminout
      % add the tetrahedron vtk cell type to the first column
      vtkConec = [ 5*ones( nelems, 1 )     Conec(:, 5:7 )-1 ] ;
      elem2VTKCellMap = (1:nelems)' ; % default: i-to-i . Columns should be changed.
    % end


  end %if type of element

  % if strcmp( elemTypeString,'truss')
  %
  % Nodesvtk = [] ;
  % vtkNormalForces = [] ;
  % vtkStressMat    = [] ;
  %
  % if nargout>1
  %
  %   counterNodesVtk = 0 ;
  %   counterCellsVtk = 0 ;
  %   Conecvtk        = [] ;
  %
  %   % mapping from the original connectivity to the vtk cells
  %   elem2VTKCellMap = (1:nelems)' ; % default: i-to-i . Colums should be changed.
  % end
  %
  % % loop in elements
  % for i = 1:nelems
  %
  %   elemMat  = Conec(i, 4+1 ) ;
  %   elemType = elementsParamsMat( Conec(i, 4+2 ), 1 ) ;
  %   elemCrossSecParams = crossSecsParamsMat( Conec(i, 4+4 ) , : ) ;
  %
  %   nodeselem  = Conec(i,1:2)' ;
  %   dofselem   = nodes2dofs( nodeselem , 6 ) ;
  %   dispsElem  = U( dofselem ) ;
  %
  %   coordsElem         = coordsElemsMat( i, 1 : 12 )'  ;
  %   %~ coordsElem(1:2:6 ) = Nodes( nodeselem(1), : )' ;
  %   %~ coordsElem(7:2:12) = Nodes( nodeselem(2), : )' ;
  %
  %   %~ dispsElem     = zeros(12,1) ;
  %   %~ dispsElem( 6) =  .5 ;
  %   %~ dispsElem(12) = -.5 ;
  %
  %   [ xdef, ydef, zdef, conecElem, titax, titay, titaz, Rr ] = outputFrameElementPlot( coordsElem, dispsElem, elemType ) ;
  %
  %   % nodes sub-division
  %   ndivNodes = conecElem( end, 2 ) ; % last node in the sub-division
  %
  %   % Local x coords in the reference config
  %   xloc = linspace( coordsElem(1), coordsElem( 7), ndivNodes )' ;
  %   yloc = linspace( coordsElem(3), coordsElem( 9), ndivNodes )' ;
  %   zloc = linspace( coordsElem(5), coordsElem(11), ndivNodes )' ;
  %
  %   for j = 1 : (ndivNodes-1)
  %
  %     dispsSubElem       = zeros(12,1) ;
  %     dispsSubElem(1:6 ) = [ xdef(j  )-xloc(j  ) titax(j  ) ydef(j  )-yloc(j  ) titay(j  ) zdef(j  )-zloc(j  ) titaz(j  ) ]' ;
  %     dispsSubElem(7:12) = [ xdef(j+1)-xloc(j+1) titax(j+1) ydef(j+1)-yloc(j+1) titay(j+1) zdef(j+1)-zloc(j+1) titaz(j+1) ]' ;
  %
  %     coordSubElem             = [ xloc(j:(j+1))  yloc(j:(j+1))  zloc(j:(j+1)) ] ;
  %
  %     if elemCrossSecParams(1) == 1  || elemCrossSecParams(1) == 2 % general or rectangular section
  %       if elemCrossSecParams(1) == 1
  %         % equivalent square section using A = wy * wz
  %         auxh = sqrt( elemCrossSecParams(2) ) ;   auxb = auxh ;
  %         secc = [ 12 auxb auxh ] ;
  %       else
  %         secc = [ 12 elemCrossSecParams(2) elemCrossSecParams(3) ] ;
  %       end
  %
	% 			iniNodes = [ 1 2 3 4  ] ;
	% 			midNodes = [          ] ;
  %       endNodes = [ 5 6 7 8  ] ;
  %
  %     elseif elemCrossSecParams(1) == 3 % circular section
  %
  %       secc = [ 25 elemCrossSecParams(2) ] ;
  %
	% 			iniNodes = [ 1 2 3 4 9 10 11 12  ] ;
	% 			midNodes = [ 17 18 19 20         ] ;
  %       endNodes = [ 5 6 7 8 13 14 15 16 ] ;
  %
  %     end
  %
  %     nNodesPerExt = length( endNodes ) ;
  %     nNodesMidEle = length( midNodes ) ;
  %     totalNodes   = nNodesPerExt*2 + nNodesMidEle ;
  %
  %     if length( iniNodes ) ~= nNodesPerExt, error('create issue.'), end
  %
  %     [ NodesCell, ConecCell ] = vtkBeam2SolidConverter ( coordSubElem, secc, dispsSubElem, Rr ) ;
  %
  %     % disps
  %     DispsNodesCell = zeros( size( NodesCell )) ;
  %     DispsNodesCell(iniNodes,:) =  ones( nNodesPerExt, 1 )*(dispsSubElem(1:2:5 )') ;
  %     DispsNodesCell(endNodes,:) =  ones( nNodesPerExt, 1 )*(dispsSubElem(7:2:11)') ;
  %     if nNodesMidEle > 0
  %       avgdisps   = (dispsSubElem(1:2:5)' + dispsSubElem(7:2:11)')*.5 ;
  %       DispsNodesCell(midNodes,:) =ones( nNodesMidEle, 1 )*(avgdisps) ;
  %     end
  %
  %     if nargout > numminout % Conectivity is required as output
  %       ConecCell(:,2:end) = ConecCell(:,2:end)+counterNodesVtk ;
  %       Conecvtk = [ Conecvtk ; ConecCell ]   ;
  %
  %       elem2VTKCellMap( i, 1:9 ) = (1:(ndivNodes-1)) + counterCellsVtk ;
  %       counterCellsVtk = counterCellsVtk + ( ndivNodes - 1 ) ;
  %     end
  %
  %     % add nodes of the initial section and the displacements
  %     Nodesvtk   = [ Nodesvtk   ; NodesCell ] ; % reorders
  %     vtkDispMat = [ vtkDispMat ; DispsNodesCell ] ;
  %
  %     counterNodesVtk = counterNodesVtk + totalNodes ;
  %
  %   end % if divisions solids
  %
  %   vtkNormalForces = [ vtkNormalForces ; ones( (ndivNodes-1),1)*normalForces(i) ] ;
  %   vtkStressMat    = [ vtkStressMat    ; ones( (ndivNodes-1),1)*stressMat(i)    ] ;
  %
  % end
  %
  % vtkStress = {}                ;
  % vtkStress{3,1} = vtkStressMat ;


end % for types elem


[ cellPointData, cellCellData ] = vtkData( modelProperties.outputDir, modelProperties.problemName, vtkDispMat, {} ) ;


function [ cellPointData, cellCellData ] = vtkData( outputdir, problemName, vtkDispMat )

  % Scalars vals
  svm             = [] ;
  vecSigI         = [] ;
	vecSigII        = [] ;
	vecSigIII       = [] ;
  normalForceMat  = [] ;

  % Vectors vals
  vecI    = [] ;
	vecII   = [] ;
	vecIII  = [] ;

  cellPointData   = cell(0) ;
	cellCellData    = cell(0) ;
  cellTensorData  = cell(0) ;

  % the point data is assumed to be the displacements field
  cellPointData = cell(1,3) ;
  cellPointData{1,1} = 'VECTORS' ; cellPointData{1,2} = 'Displacements' ; cellPointData{1,3} = vtkDispMat ;

return
	if size( vtkNormalForce, 1 ) > 0
    cellCellData{1,1} = 'SCALARS' ; cellCellData{1,2} = 'Normal_Force' ; cellCellData{1,3} = vtkNormalForce ;
	end

  if size( cellStress, 1 ) > 0 && size( cellStress{3}, 1) > 0 && size( cellStress{3}, 2) > 1
    stressMat = cellStress{3} ;
    nelems    = size( stressMat, 1 ) ;

    sxx = stressMat(:, 1 ) ;
    syy = stressMat(:, 2 ) ;
    szz = stressMat(:, 3 ) ;
    tyz = stressMat(:, 4 ) ;
    txz = stressMat(:, 5 ) ;
    txy = stressMat(:, 6 ) ;

    svm = sqrt( 0.5*( (sxx-syy).^2 + (syy-szz).^2 + (szz-sxx).^2  + 6*(tyz.^2 + txz.^2 + txy.^2) ) ) ;

    for i = 1:nelems
      tensor = [ sxx(i) txy(i) txz(i) ; txy(i) syy(i) tyz(i) ; txz(i) tyz(i) szz(i) ] ;

      [vec, val] = eig(tensor) ;

      [ sigI  , indI   ] = max( diag(val) ) ;
      [ sigIII, indIII ] = min( diag(val) ) ;

      if indI ~= indIII
        indII = setdiff( [1 2 3], [indI indIII]) ;
      else
        indII = indI ;
      end

      aux = diag(val) ;
      sigII = aux(indII) ;

      % Sigma I
      vecSigI = [vecSigI ; sigI] ;
      vecI = [vecI ; vec(:,indI)'/(norm(vec(:,indI)'))*abs(sigI)] ;

      % Sigma II
      vecSigII = [ vecSigII ; sigII ] ;
      vecII = [ vecII ; vec(:,indII)'/(norm(vec(:,indII)'))*abs(sigII)] ;

      % Sigma III
      vecSigIII = [vecSigIII ; sigIII] ;
      vecIII = [vecIII ; vec(:,indIII)'/(norm(vec(:,indIII)'))*abs(sigIII)] ;

    end
    cellCellData{1,1} = 'SCALARS' ; cellCellData{1,2} = 'Von_Mises'   ; cellCellData{1,3} = svm ;
    cellCellData{2,1} = 'VECTORS' ; cellCellData{2,2} = 'vI'          ; cellCellData{2,3} = vecI ;
    cellCellData{3,1} = 'VECTORS' ; cellCellData{3,2} = 'vII'         ; cellCellData{3,3} = vecII ;
    cellCellData{4,1} = 'VECTORS' ; cellCellData{4,2} = 'vIII'        ; cellCellData{4,3} = vecIII ;
    cellCellData{5,1} = 'SCALARS' ; cellCellData{5,2} = 'SigI'        ; cellCellData{5,3} = vecSigI ;
    cellCellData{6,1} = 'SCALARS' ; cellCellData{6,2} = 'SigII'       ; cellCellData{6,3} = vecSigII ;
    cellCellData{7,1} = 'SCALARS' ; cellCellData{7,2} = 'SigIII'      ; cellCellData{7,3} = vecSigIII ;

  elseif size( cellStress, 1 ) > 0 && size( cellStress{3}, 1) > 0
    stressMat = cellStress{3} ;
    nelems    = size( stressMat, 1 ) ;

    sxx = stressMat(:, 1 ) ;

    cellCellData{1,1} = 'SCALARS' ; cellCellData{1,2} = 'SigXLoc'   ; cellCellData{1,3} = sxx ;

  end
