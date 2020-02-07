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


function [ nodesMat, conecMat ] = esmacParser( nodesMatinp, conecMatinp, physicalNames )

matInds = zeros( length( physicalNames), 5 ) ;
for i=1:size( matInds, 1)
  for j=1:5
    matInds(i,j) = str2num( physicalNames{i}( 1+(j-1)*3+(1:2)) ) ;
  end
end

nodesMat = zeros( size(nodesMatinp,1), 8 ) ;
nodesMat(:,1:3) = nodesMatinp(:, 1:3) ;
indsNZ = find( nodesMatinp(:,4) ) ;
nodesMat(indsNZ,4:end) = matInds( nodesMatinp(indsNZ,4), :) ;

conecMat = zeros( size(conecMatinp,1), 9 ) ;
conecMat(:,1:4) = conecMatinp(:, 1:4) ;


indsNZ = find( conecMatinp(:,5) ) ;
conecMat(indsNZ,5:end) = matInds( conecMatinp(indsNZ, 5), :) ;
