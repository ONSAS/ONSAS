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


function [nodesMat conecMat] = txt2ONSAS(fname)
  fid = fopen(fname, 'r') ;

  % Nodes title
  title = fgets(fid) ;
  nodesVec = fscanf(fid, '%g'); 
  % Conec title
  title = fgets(fid) ;
  conecVec = fscanf(fid, '%g');
  fclose(fid) ;
  
  aux = length(nodesVec);
  nnodes = aux/5 ;
  nodesMat = [] ;
  for i = 1:nnodes
    nodesMat = [nodesMat ; nodesVec((i-1)*5+1:i*5)'] ;
  end
  
  aux = length(conecVec) ;
  nelems = aux/7 ;
  conecMat=[];
  for i = 1:nelems
    conecMat = [conecMat ; conecVec((i-1)*9+1:i*9)'] ;
  end

end
