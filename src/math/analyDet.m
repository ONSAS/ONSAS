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
 
function detVal = analyDet(A)

detVal =   A(1,1)*A(2,2)*A(3,3)  ...
          - A(1,1)*A(2,3)*A(3,2)  ...
          - A(1,2)*A(2,1)*A(3,3)  ...
          + A(1,2)*A(2,3)*A(3,1)  ...
          + A(1,3)*A(2,1)*A(3,2)  ...
          - A(1,3)*A(2,2)*A(3,1) ;
