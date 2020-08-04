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


% function for writing vtk files of deformed configurations of structures.
% Creates the file filename with the nodes coordinates given in nodes,
% the conectivity given in conect and with the point and element data
% given in cellPointData and cellCellData, respectively.

function vtkWriter( filename, nodes, conect, cellPointData, cellCellData )


format_A = '%18.6e' ;
format_B = '%20.8e' ;

% sizes
Nnodes = size(nodes,1) ;
nelem = size(conect,1) ;

fid = fopen( filename ,'w+') ;

% header
fprintf(fid,'# vtk DataFile Version 2.0') ;
fprintf(fid,'\n')  ;
%
fprintf(fid, filename ) ;
fprintf(fid,'\n') ;
%
fprintf(fid,'ASCII') ;
fprintf(fid,'\n') ;


% node coordinates matrix
fprintf(fid,'DATASET UNSTRUCTURED_GRID') ;
fprintf(fid,'\n') ;
fprintf(fid,'\n') ;
%
fprintf(fid,'POINTS %g float', Nnodes ) ;
fprintf(fid,'\n') ;
%
for i=1:Nnodes,
  fprintf(fid,[ ' ' format_A ' ' format_A ' ' format_A ' \n'], nodes( i,: )  ) ;
end
%
fprintf(fid,'\n') ;

nTotNumbers = 0;
for i=1:nelem
  if (conect(i,1) == 1) || (conect(i,1) == 2)
    nTotNumbers = nTotNumbers + 3 ;
  elseif (conect(i,1) == 3) || (conect(i,1) == 4)
    nTotNumbers = nTotNumbers + 5 ;
  elseif (conect(i,1) == 12)
    nTotNumbers = nTotNumbers + 8+1 ;
  elseif (conect(i,1) == 25) 
    nTotNumbers = nTotNumbers + 20+1 ;
  end
end

fprintf(fid,'CELLS %g %g', nelem, nTotNumbers ) ;
fprintf(fid,'\n') ;
for i=1:nelem,
  if (conect(i,1) == 1) || (conect(i,1) == 2)
    fprintf(fid,'2 %g %g \n', conect(i,1+(1:2))-1 ) ;
  
  elseif (conect(i,1) == 3) || (conect(i,1) == 4)
    fprintf(fid,'4 %g %g %g %g \n', conect(i,1+(1:4))-1 ) ;
  
  elseif (conect(i,1) == 12)
    fprintf(fid,'8 %g %g %g %g %g %g %g %g \n', conect(i,1+(1:8))-1 ) ;
  
  elseif (conect(i,1) == 25)
    fprintf(fid, '20 %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n', conect(i,1+(1:20))-1 ) ;
  end
end
%
fprintf(fid,'\n') ;


fprintf(fid,'CELL_TYPES %g\n', nelem ) ;
for i=1:nelem,
  if (conect(i,1) == 2) || (conect(i,1) == 3)
    fprintf(fid,'3\n') ;
  elseif (conect(i,1) == 4) 
    fprintf(fid,'10\n') ;
  elseif (conect(i,1) == 5) 
    fprintf(fid,'9\n') ;
  elseif (conect(i,1) == 12) 
    fprintf(fid,'12\n') ;
  elseif (conect(i,1) == 25)
    fprintf(fid,'25\n') ;
  end
end  
fprintf(fid,'\n') ;

%
% Cell PointData 
%
fprintf(fid, 'POINT_DATA  %8i \n' , Nnodes ) ;
for k = 1:size(cellPointData,1)
  auxtype = cellPointData(k,1) ;
  auxstr  = cellPointData(k,2) ;
  auxdata = cell2mat(cellPointData(k,3)) ;
  if strcmp(auxtype, 'VECTORS')
    fprintf(fid, ['VECTORS ' '%s' ' float\n'], [auxstr{:}] ) ;
    for i = 1 :Nnodes
      fprintf(fid,[' ' format_B ' ' format_B ' ' format_B ' \n'], auxdata(i,:) ) ;
    end
    fprintf(fid,'\n') ;
  elseif strcmp(auxtype, 'SCALARS')
    fprintf(fid,['LOOKUP_TABLE default\n']) ;
    % Data ...
  end  
end

%
% Cell CellData
%

if size( cellCellData, 1 ) > 0 
  fprintf(fid, 'CELL_DATA  %8i \n' , nelem ) ;
  for k = 1:size(cellCellData,1)
    auxtype = cellCellData(k,1) ;
    auxstr  = cellCellData(k,2) ;
    auxdata = cell2mat(cellCellData(k,3)) ;

    if strcmp(auxtype, 'SCALARS')
      fprintf(fid, ['SCALARS ' '%s' ' float 1\n'], [auxstr{:}]) ;
      fprintf(fid,['LOOKUP_TABLE default\n']) ;
      for i = 1:nelem
        fprintf(fid,[' ' format_B ' \n'], auxdata(i) ) ;
      end
      fprintf(fid,'\n') ;
    elseif strcmp(auxtype, 'VECTORS')
      fprintf(fid, ['VECTORS ' '%s' ' float\n'], [auxstr{:}] ) ;
      for i = 1:nelem
        fprintf(fid,[' ' format_B ' ' format_B ' ' format_B ' \n'], auxdata(i,:) ) ;
      end
      fprintf(fid,'\n') ;
    end
  end 
end  


  %~ fprintf(fid,['VECTORS test float\n']) ;
  %~ for i = 1 :nelem
      %~ tensor = [ sxx(i) txy(i) txz(i) ; txy(i) syy(i) tyz(i) ; txz(i) tyz(i) szz(i) ] ;
      %~ [vec,val] = eig(tensor);
      %~ [mval,ind] = max(abs(diag(val))) ;
      
    %~ fprintf(fid,[' ' format_B ' ' format_B ' ' format_B ' \n'], (vec(:,ind)/norm(vec(:,ind))*mval )' ) ;
  %~ end
  %~ fprintf(fid,'\n') ;


fclose(fid);

