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
 
% Function that converts a vector of nodes indices to a vector of the
% corresponding degrees of freedom.

function [dofs] = nodes2dofs( nodes , degreespernode )
nodes= nodes(:) ;
n    = length(nodes);
dofs = zeros( n*degreespernode , 1 ) ;
for i=1:n
  dofs( (i-1)*degreespernode + (1:degreespernode) ) = [ (degreespernode*(nodes(i)-1))+(1:degreespernode) ]' ;
end
