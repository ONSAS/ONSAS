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
 
% --------------------------------------------------------------------------------------------------

function R=RotationMatrix(ndofpnode,locglomat)

R = zeros ( ndofpnode*2, ndofpnode*2 ) ;

for k=1:2
    aux1 = [(1+(k-1)*6) : 2 : (1+(k-1)*6+4)] ; % displacement dofs
    aux2 = [(2+(k-1)*6) : 2 : (2+(k-1)*6+4)] ; % rotation     dofs
    R ( aux1, aux1 ) = locglomat ;
    R ( aux2, aux2 ) = locglomat ;
end
