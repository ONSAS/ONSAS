%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.


% Octave function for nodes-to-dof conversion
function [dofs] = nodes2dofs( nodes , degreespernode )
n    = length(nodes);
dofs = zeros( n*degreespernode , 1 ) ;
for i=1:n
  for j=1:degreespernode
    dofs( (i-1)*degreespernode + j ) = degreespernode * (nodes(i)-1) + j  ;
  end
end
