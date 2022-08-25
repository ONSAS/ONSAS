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

function [x21, d21, l, l0] = corotLenCoords(xs ,dg)

  x21 = xs(4:6) - xs(1:3) ;
  d21 = dg(7:9) - dg(1:3) ;

  l0 = sqrt( ( x21       )' * ( x21       ) ) ; % reference length
  l  = sqrt( ( x21 + d21 )' * ( x21 + d21 ) ) ; % deformed length

end