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
 
function [ nodesMat, conecMat ] = meshFileReader( fileName )

fileExtension = fileName( (end-2):end ) ;

if strcmp( fileExtension , 'msh' )
  %md reads data from msh file
  [ nodesMatinp, conecMatinp, physicalNames ] = mshFormatReader( fileName ) ;

  %md converts strings to integers indexes matrix
  matInds = zeros( length( physicalNames), 4 ) ;
  for i=1:size( matInds, 1)
    for j=1:4
      matInds(i, j) = str2num( physicalNames{i}( 1+(j-1)*3+(1:2)) ) ;
    end
  end

  nodesMat = zeros( size( nodesMatinp, 1 ), 3 ) ;
  nodesMat(:,1:3) = nodesMatinp(:, 1:3) ;

  % adds melcs parameters to nodes with params defined
  %~ indsNZ = find( nodesMatinp(:,4) ) ;
  %~ nodesMat(indsNZ,4:end) = matInds( nodesMatinp(indsNZ,4), :) ;

  conecMat = zeros( size( conecMatinp, 1 ), 4+4 ) ;
  conecMat( :, 4+(1:4)) = conecMatinp(:, 1:4) ;

  %md adds MEBI parameters to elements with params defined
  indsNZ = find( conecMatinp(:,5) ) ;
  conecMat( indsNZ, 1:4 ) = matInds( conecMatinp( indsNZ, 5), :) ;

elseif strcmp( fileExtension , 'dxf' )
  [ nodesMat, conecMat ] = dxfReader( fileName ) ;
%
else
  error('extension not implemented yet. Please report an issue.')
end
