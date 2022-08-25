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
 
function [ Ro, lengthElem ] = beamRefConfRotMat( x ) ;

  assert( iscolumn(x), 'coordinates must be in a column vector.')

  lengthElem = norm(x) ;

  exL = x / lengthElem ;

  eyG = [0 1 0]' ;    ezG = [0 0 1]' ;

  % Vector normal to beam in reference configuration
  if ( abs( exL(1) ) > 1e-8*lengthElem ) || ( abs( exL(2) ) > 1e-8*lengthElem ) ; % if exL it is not ezG
    aux = cross( ezG, exL ) ;
    eyL = aux / norm( aux ) ;
  else
    eyL = eyG ;
  end
  ezL = cross( exL, eyL ) ;

  Ro  = [ exL eyL ezL ] ;

end
