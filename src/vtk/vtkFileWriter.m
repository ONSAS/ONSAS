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
 

%md function for writing vtk files of deformed configurations of structures.
%md Creates the file filename with the nodes coordinates given in nodes,
%md the conectivity given in conect and with the point and element data
%md given in cellPointData and cellCellData, respectively.

function vtkFileWriter( filename, nodes, conect, cellPointData, cellCellData )

format_A = '%18.6e' ;
format_B = '%20.8e' ;

% sizes
Nnodes = size( nodes , 1 ) ;
nelem  = size( conect, 1 ) ;

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

% count the number of numbers to print in the CELLS block
nTotNumbers = 0;
cellTypes = unique( conect(:,1) ) ;

for i = 1:length(cellTypes)
  % count elements with the current cell type
  nelemCurrType = sum( conect( :, 1 )== cellTypes(i) ) ;

  if ( cellTypes(i) == 5 ) % if it is triangle add (3+1)*nelemCurrType
    nTotNumbers = nTotNumbers + (3+1)*nelemCurrType ;

  elseif ( cellTypes(i) == 10 ) % if it is tetrahedron add (4+1)*nelemCurrType
    nTotNumbers = nTotNumbers + (4+1)*nelemCurrType ;

  elseif ( cellTypes(i) == 12 ) % if it is tetrahedron add (4+1)*nelemCurrType
    nTotNumbers = nTotNumbers + (8+1)*nelemCurrType ;

  elseif ( cellTypes(i) == 25 ) % if it is tetrahedron add (4+1)*nelemCurrType
    nTotNumbers = nTotNumbers + (20+1)*nelemCurrType ;

  end % if celltype
end % for celltypes


fprintf(fid,'CELLS %g %g', nelem, nTotNumbers ) ;
fprintf(fid,'\n') ;
for i=1:nelem,
  % if (conect(i,1) == 1) || (conect(i,1) == 2)
  %   fprintf(fid,'2 %g %g \n', conect(i,1+(1:2))-1 ) ;

  if (conect(i,1) == 5) %triangle
    fprintf(fid,'3 %g %g %g \n', conect(i,1+(1:3)) ) ;

  elseif (conect(i,1) == 10) % tetrahedron
    fprintf(fid,'4 %g %g %g %g \n', conect(i,1+(1:4)) ) ;

  elseif (conect(i,1) == 12)
    fprintf(fid,'8 %g %g %g %g %g %g %g %g \n', conect(i,1+(1:8)) ) ;

  elseif (conect(i,1) == 25)
    fprintf(fid, '20 %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n', conect(i,1+(1:20)) ) ;
  end
end
%
fprintf(fid,'\n') ;


%
% Cell type numbers
%
fprintf(fid,'CELL_TYPES %g\n', nelem ) ;
for i=1:nelem,
  fprintf(fid,'%3i\n',conect(i,1) ) ;
end
fprintf(fid,'\n') ;

%
% Cell PointData
%
if size( cellPointData, 1 ) > 0
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
      % Data
      fprintf(fid, ['SCALARS ' '%s' ' float 1\n'], [auxstr{:}] ) ;
      fprintf(fid,['LOOKUP_TABLE default\n']) ;
      for i = 1 :Nnodes
        fprintf(fid,[' ' format_B '\n'], auxdata(i,:) ) ;
      end
    end
  end
end

%
% Cell CellData
%

if size( cellCellData, 1 ) > 0
  fprintf(fid, 'CELL_DATA  %8i \n' , nelem ) ;
  for k = 1:size( cellCellData, 1 )
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
