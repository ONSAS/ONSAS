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

function [numNodes, dofsStep] = elementTypeInfo ( elemType )

if strcmp( elemType, 'node');
  numNodes = 1 ;
  dofsStep = 2 ;

elseif strcmp( elemType, 'truss') || strcmp( elemType, 'edge')
  numNodes = 2 ;
  dofsStep = 2 ;

elseif strcmp( elemType, 'frame')
  numNodes = 2 ;
  dofsStep = 1 ;

elseif strcmp( elemType, 'tetrahedron')
  numNodes = 4 ;
  dofsStep = 2 ;

elseif strcmp( elemType, 'triangle')
  numNodes = 3 ;
  dofsStep = 2 ;

end
