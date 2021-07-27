% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
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

function Ro = beamRefConfRotMat( x ) ;

  exL = x / norm(x) ;

  eyG = [0 1 0]' ;
  ezG = [0 0 1]' ;

  % Vector normal to beam in reference configuration
  aux = cross ( ezG, exL ) ;

  if norm( aux ) > 1e-15
    eyL = aux / norm( aux );
  else
    eyL = eyG ;
  end

  ezL = cross( exL, eyL ) ;

  Ro = [ exL eyL ezL ] ;
end
