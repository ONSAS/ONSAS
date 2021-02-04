% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
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


%function mesh2vtk( ...
%  filename , nodes , conect , cell1, cell2 ) %, ...
  %~ Nvector_point_data , vector_point_data_names , vector_point_data_values , ...
  %~ Nfield_data        ,             field_names ,              field_values )
  
  format_A = '%16.4e' ;
  format_B = '%20.8e' ;
  filename=[ outputdir problemName  '.vtk'] ;
  cell1 = cellDisp ;
  cell2 = cellStress ;
  Conec = Conec(:,1:4) ;
  % tamanios
  Nnodes = size(Nodes,1) ;
  nelem = size(Conec,1) ;

  fid = fopen( filename ,'w+') ;


  % escribe encabezado
  fprintf(fid,'# vtk DataFile Version 2.0') ;
  fprintf(fid,'\n')  ;
  %
  fprintf(fid, filename ) ;
  fprintf(fid,'\n') ;
  %
  fprintf(fid,'ASCII') ;
  fprintf(fid,'\n') ;


  % escribe matrices de coord d nodos
  fprintf(fid,'DATASET UNSTRUCTURED_GRID') ;
  fprintf(fid,'\n') ;
  fprintf(fid,'\n') ;
  %
  fprintf(fid,'POINTS %g float', Nnodes ) ;
  fprintf(fid,'\n') ;
  %
  for i=1:Nnodes,
    %
    fprintf(fid,[ ' ' format_A ' ' format_A ' ' format_A ' \n'], Nodes(i,:) ) ;
    %
  end
  %
  fprintf(fid,'\n') ;



  fprintf( fid , 'CELLS %g %g', nelem , nelem*(4+1) ) ;
  fprintf(fid,'\n') ;
    
  %
  for i=1:nelem,
    fprintf(fid,'%g %g %g %g %g \n', 4, Conec(i,:)-1 ) ;
  end
  %
  fprintf(fid,'\n') ;
  
  
   fprintf (fid,'CELL_TYPES %g \n', nelem)
      for i=1:nelem
        fprintf(fid,'10\n')
      end


  fprintf(fid, 'POINT_DATA  %8i \n' , Nnodes ) ;
      %~ fprintf(fid,['SCALARS  name  float 1 \n']) ;
      %~ fprintf(fid,'LOOKUP_TABLE default\n') ;

      %~ for i=1:Nnodes
        %~ % esto es para graficar las flechas de desp tambien en la mlla indeformada
        %~ fprintf(fid,['  ' format_B '  \n'], auxvalores(i) ) ;
      %~ end
  
      %~ fprintf(fid,'\n') ;


      fprintf(fid,['VECTORS ' cell1{2} ' float \n']) ;

      for i=1:Nnodes
        % esto es para graficar las flechas de desp tambien en la mlla indeformada
        fprintf(fid,['  ' format_B '  ' format_B '  ' format_B '  \n'], cell1{3}((3*i-2):(3*i)) ) ;
      end
      %lo agregue yo Joaquin pero esta para cambiar 
      fprintf(fid,['VECTORS ' 'Cargas' ' float \n']) ;
      for i = 1:Nnodes
        if ismember(i,nodalVariableLoads(:,1))
          a=find(nodalVariableLoads(:,1)==i);
          fprintf(fid,['  ' format_B '  ' format_B '  ' format_B '  \n'], [nodalVariableLoads(a,1) nodalVariableLoads(a,3) nodalVariableLoads(a,5)] ) ;
        else
          fprintf(fid,['  ' format_B '  ' format_B '  ' format_B '  \n'], [ 0 0 0 ] ) ;
        end  
      end  
  
      fprintf(fid,'\n') ;


   
  
  %~ % escribe cell data
  %~ if Nfield_data > 0
    fprintf(fid, 'CELL_DATA %8i \n' , nelem ) ;
    fprintf(fid, ['TENSORS ' cell2{2} ' float \n'] ) ;
    matdata = cell2{3};
    for i=1:size( matdata,1)
      vecdata = matdata(i,:) ;

   s11 = vecdata(1);
   s22 = vecdata(2);
   s33 = vecdata(3);

   s23 = vecdata(4);
   s32 = s23;

   s13 = vecdata(5);
   s31 = s13;
 
   s12 = vecdata(6);
   s21 = s12;
 svm = sqrt(  1/2 * ( (s11-s22)^2 + (s22-s33)^2 + (s33-s11)^2 + 6*( s23^2 + s31^2 + s12^2 ) )  ) ;   
  


      fprintf(fid,[' ' format_B ' ' format_B ' ' format_B ' ' ...
       format_B ' ' format_B ' ' format_B ' ' ...
       format_B ' ' format_B ' ' format_B ' ' ...
       ' \n'], ...
      vecdata(1), vecdata(6), vecdata(5), ...
      vecdata(6), vecdata(2), vecdata(4), ...
      vecdata(5), vecdata(4), svm  ) ;
    end
      fprintf(fid,'\n') ;

    %~ for j=1:Nfield_data
      
      %~ auxvalores = cell2mat(field_values(j))  ;
      %~ auxcomps   = size(auxvalores,2) ;
  
      %~ fprintf(fid, '%s %3i %8i float \n' , cell2mat(field_names(j)) , auxcomps , nelem ) ;
  
      %~ if size(auxvalores,1)~=nelem, size(auxvalores,1), nelem, error('tamanio de field'), end
      
      %~ for indelem =1:nelem
        %~ for k=1:auxcomps
          %~ fprintf( fid, [' ' format_B ' '],auxvalores(indelem,k) ) ;
        %~ end
        %~ fprintf( fid, '\n');
      %~ end
    %~ end
    %~ fprintf(fid, 'POINT_DATA %8i \n' , Nnodes ) ;
  %~ end
  
  fclose(fid);
