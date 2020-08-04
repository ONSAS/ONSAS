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

function [ cellPointData, cellCellData, filename ] = vtkData( outputdir, problemName, indplot, vtkNormalForce, cellStress, vtkDispMat )

  filename = [ outputdir problemName '_' sprintf('%04i',indplot-1) '.vtk'] ;

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
  
	if size(vtkNormalForce,1) > 0
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
    
    cellCellData{1,1} = 'SCALARS' ; cellCellData{1,2} = 'sig x loc'   ; cellCellData{1,3} = sxx ;
  
  end

  cellPointData = cell(1,3) ;
  cellPointData{1,1} = 'VECTORS' ; cellPointData{1,2} = 'Displacements' ; cellPointData{1,3} = vtkDispMat ;


